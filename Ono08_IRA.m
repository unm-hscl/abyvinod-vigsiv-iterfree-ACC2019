%% Ono_IRA 2008 code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 10; 


% Maximum/minimum bound on input: 

ulim = 0.2; 

% max risk: 

Delta = 0.05;

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
    x0 = [0.1;0];

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
    mean_X_sans_input = Ad * x0 + mean_concat_disturb;
    cov_X_sans_input = cov_concat_disturb;

% Generate bounds: 
    h = [-1 0; 1 0;];
    hbig = kron(eye(T+1),h);
    g = [0.1; 0.1];
    gbig = repmat(g,T+1,1);
    
    
%% Ono's converge alpha value
    alpha_on_iter = @(n) 0.7 * (0.98)^n;

    %% House-keeping for the event where the bisection fails to take off
    lb_stoch_reach_avoid = -1;
    optimal_input_vector = nan(size(u,1) * T,1);
    
    %% Define when a slack variable will be set to zero
    myeps = 1e-6;
    %% Desired accuracy
    desired_accuracy = 0.1; 

    %% Compute M --- the number of polytopic halfspaces to worry about
    no_linear_constraints = size(hbig,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = sqrt(diag(hbig*cov_X_sans_input*hbig'));

    %% Bisection over Delta between 0 and 0.5    
    Delta_lb = 0;
    Delta_ub = 0.5;
while (Delta_ub - Delta_lb) > desired_accuracy
        Delta = (Delta_ub + Delta_lb)/2;

        %% Prepration for iterated risk allocation
        % Variables (re)initialized
        iter_count = 0;                 % No. of iterations done in Ono's risk
                                        % allocation algorithm
        opt_value_prev = 10000;         % Previous optimal value --- exit cond.
        opt_value = 0;                  % Optimal value --- |slack variables|_1
        % Given Delta, construct delta_i as Delta/M
        delta_vec = Delta/no_linear_constraints * ones(no_linear_constraints,1);

        % Converge to one more decimal precision
        while abs(opt_value - opt_value_prev) >= desired_accuracy/10
            %% Store the previous optimal value
            opt_value_prev = opt_value;

            %% Solve the feasibility problem
            % Construct the back-off (Hessem's term) in the constraints
            scaled_norminv=sigma_vector.*...
                              norminv(ones(no_linear_constraints,1)- delta_vec);
            cvx_begin 
                variable U_vector(size(B,2)*T,1);
                variable mean_X(size(mean_X_sans_input,1), 1);
                variable slack_variables(no_linear_constraints, 1);
                minimize norm(slack_variables,1);
                subject to
                    mean_X == mean_X_sans_input + Bd*U_vector; 
                    abs(U_vector) <= ulim; 
                    abs(mean_X)<=0.01;
                    hbig*mean_X<=gbig-sigma_vector.*norminv(1-delta_vec) + slack_variables; 
                    slack_variables >= 0;
            cvx_end
            opt_value = norm(slack_variables,1);

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
                hbig(inactive_indx,:) * mean_X;
            % Update inactive delta
            delta_vec(inactive_indx) = alpha * delta_vec(inactive_indx) +...
                (1 - alpha) * (1 - normcdf(...
                    correction_deltas./sigma_vector(inactive_indx)));

            %% Collect the headroom gained: \delta_residual \gets \Delta -
            %% \sum_i \delta_i
            delta_residual = Delta - sum(delta_vec);

            %% For all the active constraints, distribute the remaining headroom
            active_indx = find(slack_variables >= myeps);
            % \delta_i \gets \delta_i + \delta_residual/N_active 
            delta_vec(active_indx) = delta_vec(active_indx) +...
                delta_residual / N_active;

            %% Update the iteration count
            iter_count = iter_count + 1;
        end
        fprintf('Done with Delta: %1.4f, N_active: %2d ',...
                Delta,...
                N_active);
        if N_active == 0
            % If no constraint is active, then a feasible input policy has been
            % found --- Dream higher (decrease delta)
            Delta_ub = Delta;
            lb_stoch_reach_avoid = 1-Delta;
            optimal_input_vector = U_vector;
            disp(' Aim for more probability');
        else
            % All constraints are active, then no solution could be found ---
            % Dream lower (increase delta)
            Delta_lb = Delta;
            disp(' Aim for less probability');
        end
end   
    
% % Generate intial delta and do some housekeeping: 
%     del = Delta/((T+1)*size(h,1)); 
%     delta = repmat(del,size(h,1)*(T+1),1);
%     Jp =0.1;
%     epsil = 1e-6;
%     sigma_vector = sqrt(diag(hbig*cov_X_sans_input*hbig'));
% while abs(norm(slack_variables,1)-Jp)<epsil
%     if Jp ~= 0
%         Jp = norm(slack_variables,1);
%     end
%     cvx_clear
%     cvx_begin
%         variable u(size(B,2)*T)
%         variable mean_X_w_input(size(A,2)*(T+1))
%         variable slack_variables(size(hbig,1)) nonnegative
%         minimize norm(slack_variables,1)%((sum(abs(u))) + (sum(abs(mean_X_w_input))))
% 
%         subject to
% 
%             mean_X_w_input == mean_X_sans_input + Bd*u; 
%             abs(u) <= ulim; 
% %             abs(mean_X_w_input)<=0.01;
%             hbig*mean_X_w_input<=gbig-sigma_vector.*norminv(1-delta) + slack_variables; 
%     cvx_end
%  end
%     