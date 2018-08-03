%% Ono_IRA 2008 code
% Coder: Abraham Vinod and Vignesh Sivaramakrishnan

disp('---------OnoIRA2008-----------')
disp(' ')

ono_opt_mean_X = nan(size(Ad,1),1);
ono_opt_val = nan;
ono_opt_input_vector = nan(size(Bd,2),1); 
ono_opt_value_array = nan;

if Delta>0.5
    warning('Skipping Ono''s formulation since Delta is not <0.5');
else
    
% Generate bounds: 
    hbig = kron(eye(T),h);
    gbig = kron(g,[1,1])';

disp('================================')
disp('Finding feasible deltas:');
disp('================================')



% Ono's converge alpha value
    alpha_on_iter = @(n) 0.7 * (0.98)^n;

% House-keeping for the event where the bisection fails to take off
    lb_stoch_reach_avoid = -1;
    optimal_input_vector = nan(size(Bd,1),1);

% Define when a slack variable will be set to zero
    myeps = 1e-6;
% Desired accuracy
    desired_accuracy = 0.001; 

% Compute M --- the number of polytopic halfspaces to worry about
    no_linear_constraints = size(hbig,1);


%% Prepration for iterated risk allocation to find feasible deltas
% Variables (re)initialized
    iter_count = 0;                 % No. of iterations done in Ono's risk
                                    % allocation algorithm
    % Given Delta, construct delta_i as Delta/M
        delta_vec = Delta/no_linear_constraints *...
            ones(no_linear_constraints,1);


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
        variable U_vector(size(Bd,2),1);
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
    tot_time1(iter_count+1) = cvx_cputime;

    % Number of active/infeasible constraints via complementary
    % slackness --- non-zero slack variables imply infeasible \delta_i
        N_active=nnz(slack_variables >= myeps);

    % Break if N_active is zero or all are active
        if N_active == 0 || N_active == no_linear_constraints
            break;
        end

    % Compute \alpha value --- decides on how quickly inactive \delta_i
    % are tightened: 
        alpha = alpha_on_iter(iter_count); 

    % For all the inactive constraints, \delta_i^+ \gets \alpha
    % \delta_i + (1-\alpha) (1-normcdf(g-hx/sigma_vec(i)))
        inactive_indx = find(slack_variables < myeps);
    % Compute relevant g-hx^\ast
        correction_deltas = gbig(inactive_indx) - ...
            hbig(inactive_indx,:) * mean_X(3:end);
    % Update inactive delta
        delta_vec(inactive_indx) = alpha * delta_vec(inactive_indx) +...
            (1 - alpha) * (1 - normcdf(...
                correction_deltas./sigma_vector(inactive_indx)));

    % Collect the headroom gained: \delta_residual \gets \Delta -
    % \sum_i \delta_i
        if Delta - sum(delta_vec) < 0
            disp('Breaking early since we got a negative residual');
            error('Given Delta requirement is not feasible.');
        end
        delta_residual =  Delta - sum(delta_vec);

    % For all the active constraints, distribute the remaining headroom
        active_indx = find(slack_variables >= myeps);
    % \delta_i \gets \delta_i + \delta_residual/N_active 
        delta_vec(active_indx) = delta_vec(active_indx) +...
            delta_residual / N_active;

    % Update the iteration count
        iter_count = iter_count + 1;
end


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


disp('================================')
disp('Now working on control objective:');
disp('================================')

%% Prepration for main iterated risk allocation problem
% Variables (re)initialized
    iter_count = 1;                 % No. of iterations done in Ono's risk
                                    % allocation algorithm
    opt_value_prev = 10000;         % Previous optimal value --- exit cond.
    opt_value = 0;                  % Optimal value --- |slack variables|_1
    % Given delta_vec that is feasible we will now solve for a control
    % objective. 


%% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
sigma_vector = sqrt(diag(hbig*cov_X_sans_input(3:end,3:end)*hbig'));
ono_opt_value_array(1) = (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
    sum(abs(mean_X(3:2:end)-xtarget))/(2*g(1)*T));

% Converge to one more decimal precision
    while abs(opt_value_prev - opt_value) >  1e-4
        %% Store the previous optimal value
        opt_value_prev = opt_value;

        %% Solve the control problem
        % Construct the back-off (Hessem's term) in the constraints
        scaled_norminv=sigma_vector.*...
                          norminv(ones(no_linear_constraints,1)- delta_vec);
        cvx_precision BEST
        tstart = tic;
        cvx_begin quiet
            variable U_vector(size(Bd,2),1);
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
        tot_time2(iter_count) = cvx_cputime;
        if ~strcmp(cvx_status, 'Solved')
            keyboard
        end
        opt_value = cvx_optval;
        ono_opt_value_array(iter_count+1) = opt_value;

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
            disp('Breaking early since we get a negative residual');
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

ono_opt_mean_X = mean_X;
ono_opt_val = ono_opt_value_array(end);
ono_opt_input_vector = U_vector; 

fprintf('Done with Delta: %1.4f, N_active: %2d\n\n',...
        Delta,...
        N_active);

fprintf('Total CVX Run Time: %1.4f seconds\n',...
    sum(tot_time1)+sum(tot_time2))
disp('------------------------------------')
fprintf('Total CVX Solve Time: %1.4f seconds\n\n',sum(time_to_solve1)+sum(time_to_solve2))
disp('------------------------------------')
fprintf('Total Run Time: %1.4f seconds\n', sum(tot_time1)+sum(tot_time2)+sum(time_to_solve1)+sum(time_to_solve2))
end

