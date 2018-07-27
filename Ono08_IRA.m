%% Ono_IRA 2008 code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 30; 


% Maximum/minimum bound on input: 

ulim = 1; 

% System matrices: 
    % Sampling time of the discrete system:
        delT = 0.25;
    % A matrix of the 2D double integrator:
        A = [ 1 delT; 
              0 1   ;];
          
    % B matrix of the 2D double integrator:
        B = [delT^2/2; delT]; 
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = ones(size(B,1)*T,1);
       w = zeros(size(x,1)*T,1);
       
   % Define the matrices of the discrete time system: 
       Ad = [];
       Bd = zeros(size(A,2)*(T+1),size(u,1)); 
       Gd = zeros(size(A,2)*T,size(w,1));
       
    for i = 1:(T+1)
           
           Ad = [Ad; A^(i-1);];
           
    end 

       for i = 1:(T+1)
           if i == 1
               
                Ab{i}= zeros(size(A,2),size(B,2)); 
               
           else
               
                Ab{i}=A^(i-2)*B; 
                
           end
    
       end

       No = numel(Ab);
       X = tril(true(No));
       Y = hankel(ones(No,1));
       [R,C] = find(Y);
       U = accumarray(R,find(X),[],@(n){n});
       Z = cell(No);
       for k = 1:No
           Z(U{k}) = Ab(k);
       end
       Z(~X) = {zeros(size(Ab{1}))};
       Bd = cell2mat(Z);
       Bd = Bd(:,1:end-1);
       
% Randomly generate the disturbance vector from the standard normal.
    cov_mat_diag = 0.001*diag([1 0;]); 
    cov_mat = kron(eye(T+1),cov_mat_diag); 
    
% Initial conditions: 
    x0 = [0.4;0];
    xtarget = linspace(-0.495,0,T)';
% Generate nominal x (Note this is a code snippet taken from SReachTools):
%     mean_concat_disturb = kron(ones(time_horizon,1), ...
%                                sys.dist.parameters.mean);
%     cov_concat_disturb  = kron(eye(time_horizon), ...
%                                sys.dist.parameters.covariance);
%     mean_X_sans_input = Ad * initial_state.parameters.mean + G *...
%         mean_concat_disturb;
%     cov_X_sans_input = Z * initial_state.parameters.covariance * Z' + ...
%         G * cov_concat_disturb * G';
    mean_concat_disturb = kron(ones(T+1,1),[0;0]);
    cov_concat_disturb  = kron(eye(T+1), cov_mat_diag);
    cov_concat_disturb(1:2,1:2) = zeros(2);
    mean_X_sans_input = Ad * x0 + mean_concat_disturb;
    cov_X_sans_input = cov_concat_disturb;

% Generate bounds: 
    h = [-1 0; 1 0;];
    hbig = kron(eye(T),h);
%     g = 0.5*[1;1];
%     gbig = repmat(g,T,1);
    g = linspace(0.5,0.1, T);
    gbig = kron(g,[1,1])';
    
disp('================================')
disp('Finding feasible deltas:');
disp('================================')
   
    
    
%% Ono's converge alpha value
alpha_on_iter = @(n) 0.7 * (0.98)^n;

%% House-keeping for the event where the bisection fails to take off
lb_stoch_reach_avoid = -1;
optimal_input_vector = nan(size(u,1) * T,1);

%% Define when a slack variable will be set to zero
myeps = 1e-8;
%% Desired accuracy
desired_accuracy = 0.001; 

%% Compute M --- the number of polytopic halfspaces to worry about
no_linear_constraints = size(hbig,1);

Delta = 0.05;

%% Prepration for iterated risk allocation
% Variables (re)initialized
iter_count = 0;                 % No. of iterations done in Ono's risk
                                % allocation algorithm
% Given Delta, construct delta_i as Delta/M
delta_vec = Delta/no_linear_constraints * ones(no_linear_constraints,1);


%% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
sigma_vector = sqrt(diag(hbig*cov_X_sans_input(3:end,3:end)*hbig'));
N_active = no_linear_constraints;
% Converge to one more decimal precision
while N_active > 0 && iter_count < 50
    %% Store the previous optimal value

    %% Solve the feasibility problem
    % Construct the back-off (Hessem's term) in the constraints
    scaled_norminv=sigma_vector.*...
                      norminv(ones(no_linear_constraints,1)- delta_vec);
    tstart = tic;
    cvx_begin quiet
        variable U_vector(size(B,2)*T,1);
        variable mean_X(size(mean_X_sans_input,1), 1);
        variable slack_variables(no_linear_constraints, 1) nonnegative;
        minimize norm(slack_variables,1);
        subject to
            mean_X == mean_X_sans_input + Bd*U_vector; 
            abs(U_vector) <= ulim; 
            hbig*mean_X(3:end)<=gbig - scaled_norminv + slack_variables;
    t1 = toc(tstart);
    cvx_end;
    t2 = toc(tstart);
    time_to_solve1(iter_count+1) = t2 - t1;

    %% Number of active/infeasible constraints via complementary
    %% slackness --- non-zero slack variables imply infeasible \delta_i
    N_active=nnz(slack_variables >= myeps);

    %% Break if N_active is zero or all are active
    if N_active == 0 || N_active == no_linear_constraints
        break;
    end

    %% Compute \alpha value --- decides on how quickly inactive \delta_i
    %% are tightened: 
    alpha = alpha_on_iter(iter_count); 

    %% For all the inactive constraints, \delta_i^+ \gets \alpha
    %% \delta_i + (1-\alpha) (1-normcdf(g-hx/sigma_vec(i)))
    inactive_indx = find(slack_variables < myeps);
    % Compute relevant g-hx^\ast
    correction_deltas = gbig(inactive_indx) - ...
        hbig(inactive_indx,:) * mean_X(3:end);
    % Update inactive delta
    delta_vec(inactive_indx) = alpha * delta_vec(inactive_indx) +...
        (1 - alpha) * (1 - normcdf(...
            correction_deltas./sigma_vector(inactive_indx)));

    %% Collect the headroom gained: \delta_residual \gets \Delta -
    %% \sum_i \delta_i
    if Delta - sum(delta_vec) < 0
        disp('Breaking early since we got a negative residual');
        break
    end
    delta_residual =  Delta - sum(delta_vec);

    %% For all the active constraints, distribute the remaining headroom
    active_indx = find(slack_variables >= myeps);
    % \delta_i \gets \delta_i + \delta_residual/N_active 
    delta_vec(active_indx) = delta_vec(active_indx) +...
        delta_residual / N_active;
%     [iter_count N_active cvx_optval delta_vec(2:2:end)']
    %% Update the iteration count
    iter_count = iter_count + 1;
end
% min_delta_true_with_motion = 1-mvncdf(-repmat(g(1),1,T),...
%     repmat(g(2),1,T),mean_X(3:2:end)', cov_X_sans_input(3:2:end,3:2:end));    

fprintf('Done with Delta: %1.4f, N_active: %2d\n',...
        Delta,...
        N_active);
if N_active == 0 
    % If no constraint is active, then a feasible input policy has been
    % found --- Dream higher (decrease delta)
    disp(' ')
    disp('Found a feasible initial delta!');

else
    error('Given Delta requirement is not feasible.');
end
% figure()
% plot(1:(T+1),mean_X(1:2:end))    
% hold on
% plot(1:(T+1), -g(1)*ones(1,T+1),'r')
% plot(1:(T+1), g(2)*ones(1,T+1),'r')
% plot(1:(T+1), xtarget*ones(1,T+1),'g')
% xlabel('time')
% ylabel('x')
% title('Init delta');
disp('================================')
disp('Now working on control objective:');
disp('================================')

%% Prepration for iterated risk allocation
% Variables (re)initialized
iter_count = 1;                 % No. of iterations done in Ono's risk
                                % allocation algorithm
opt_value_prev = 10000;         % Previous optimal value --- exit cond.
opt_value = 0;                  % Optimal value --- |slack variables|_1
% Given delta_vec that is feasible we will now solve for a control
% objective. 


%% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
sigma_vector = sqrt(diag(hbig*cov_X_sans_input(3:end,3:end)*hbig'));
input_state_ratio = 0.0001;
opt_value_array(1) = (input_state_ratio*sum(abs(U_vector))/(ulim*T) + sum(abs(mean_X(3:2:end)-xtarget))/(2*g(1)*T));
% Converge to one more decimal precision
while abs(opt_value_prev - opt_value) >  1e-4
    %% Store the previous optimal value
    opt_value_prev = opt_value;

    %% Solve the feasibility problem
    % Construct the back-off (Hessem's term) in the constraints
    scaled_norminv=sigma_vector.*...
                      norminv(ones(no_linear_constraints,1)- delta_vec);
    cvx_precision BEST
    tstart = tic;
    cvx_begin quiet
        variable U_vector(size(B,2)*T,1);
        variable mean_X(size(mean_X_sans_input,1), 1);
        minimize (input_state_ratio*sum(abs(U_vector))/(ulim*T) + sum(abs(mean_X(3:2:end)-xtarget))/(2*g(1)*T));
        subject to
            mean_X == mean_X_sans_input + Bd*U_vector; 
            abs(U_vector) <= ulim; 
            hbig*mean_X(3:end)<=gbig - scaled_norminv;
    t1 = toc(tstart);
    cvx_end;
    t2 = toc(tstart);
    time_to_solve2(iter_count) = t2 - t1;
    if ~strcmp(cvx_status, 'Solved')
        keyboard
    end
    opt_value = cvx_optval;
    opt_value_array(iter_count+1) = opt_value;
    
    %% Number of active/infeasible constraints via complementary
    %% slackness --- non-zero slack variables imply infeasible \delta_i
    N_active=nnz( abs(hbig*mean_X(3:end) - gbig + scaled_norminv)<myeps);

    %% Break if N_active is zero or all are active
    if N_active == 0
        disp('No issues with the current delta allocation');
        break;
    elseif N_active == no_linear_constraints
        error('Not possible to solve');
    end

    %% Compute \alpha value --- decides on how quickly inactive \delta_i
    %% are tightened: 
    alpha = alpha_on_iter(iter_count); 

    %% For all the inactive constraints, \delta_i^+ \gets \alpha
    %% \delta_i + (1-\alpha) (1-normcdf(g-hx/sigma_vec(i)))
    inactive_indx = find(hbig*mean_X(3:end) < gbig - scaled_norminv);
    % Compute relevant g-hx^\ast
    correction_deltas = gbig(inactive_indx) - ...
        hbig(inactive_indx,:) * mean_X(3:end);
    % Update inactive delta
    delta_vec(inactive_indx) = alpha * delta_vec(inactive_indx) +...
        (1 - alpha) * (1 - normcdf(...
            correction_deltas./sigma_vector(inactive_indx)));
    %% Collect the headroom gained: \delta_residual \gets \Delta -
    %% \sum_i \delta_i
    if Delta - sum(delta_vec) < 0
        disp('Breaking early since we got a negative residual');
        break;
    end
    delta_residual =  Delta - sum(delta_vec);

    %% For all the active constraints, distribute the remaining headroom
    active_indx = find(abs(hbig*mean_X(3:end) - gbig + scaled_norminv)<myeps);
    % \delta_i \gets \delta_i + \delta_residual/N_active 
    delta_vec(active_indx) = delta_vec(active_indx) +...
        delta_residual / N_active;
%     [iter_count N_active opt_value]
    %% Update the iteration count
    iter_count = iter_count + 1;
end
% min_delta_true_with_motion = 1-mvncdf(-repmat(g(1),1,T),repmat(g(2),1,T),mean_X(3:2:end)', cov_X_sans_input(3:2:end,3:2:end));    

fprintf('Done with Delta: %1.4f, N_active: %2d\n',...
        Delta,...
        N_active);
%%
figure(1)
hold on
plot(1:(T+1), [-gbig(1);-gbig(1:2:end)],'r')
plot(1:(T+1), [gbig(2);gbig(2:2:end)],'r')
plot(2:(T+1), xtarget,'go')
plot(1:(T+1),mean_X(1:2:end))    
xlabel('time')
ylabel('x')
title('Trajectory')

figure(2);
hold on
plot(opt_value_array)
xlabel('iteration count')
ylabel('cost')
title('Final J cost');

fprintf('Total Solve Time: %1.4f seconds\n',sum(time_to_solve1)+sum(time_to_solve2))