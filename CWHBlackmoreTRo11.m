%% Blackmore TRo 2011 Code with Ono08 Dynamics & Formulation
% Coder: Vignesh Sivaramakrishnan

% System matrices: 

    disp('---------BlackmorePC11-----------')
    disp(' ')
    fprintf('No. of particles: %d\n',N);
    
    large_constant = 5000;
% Generate the cov_matrix for optimization problem given covariance.
    cov_mat = kron(eye(T+1),covariance_disturbance); 
    
% Vectorize the target trajectory for the optimization problem. 
    xtargetbig = repmat(xtarget,1,N);

% Randomly generate the disturbance vector from the standard normal.

    mean_GdTimesW = repmat(mean_disturbance,T,N);
    
    GdTimesW = mvnrnd(mean_GdTimesW', cov_X_sans_input)';
  
    

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
    tstart = tic;
    cvx_clear
        cvx_precision BEST
    cvx_begin quiet
        variable U_vector(size(Bd,2),1);
        variable xBl(size(mean_X_sans_input,1),N);
        variable mean_X(size(mean_X_sans_input,1),1);
        variable d(N) binary;

        %minimize (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
            %sum(sum(abs(xBl(1:2:end,1:end)-xtargetbig)))/(2*g(1)*T)/N);
        minimize (sum(sum((xBl-xtargetbig).^2))/N);

        subject to
          mean_X == Ad*x0+ Bd*U_vector;

          xBl(1:end,1:N) == GdTimesW+repmat(mean_X,1,N);

          abs(U_vector) <= ulim;

          for i = 1:N
              h*xBl(:,i) - gb(:) <= large_constant*(d(i));
          end
          1/N*sum(d)<=Delta;

    t1 = toc(tstart);
    cvx_end;
    t2 = toc(tstart);
    time_to_solve = t2 - t1;
    
    if strcmpi(cvx_status,'Solved')
        blackmore_opt_mean_X = mean_X;
        blackmore_opt_val = cvx_optval;
        blackmore_opt_input_vector = U_vector;            
    else
        blackmore_opt_mean_X = nan(length(mean_X_sans_input),1);
        blackmore_opt_val = nan;
        blackmore_opt_input_vector = nan(size(Bd,2),1);         
    end
    
    fprintf('Total CVX Run Time for %1i particles: %1.4f seconds\n',...
        N,cvx_cputime)
    disp('------------------------------------')
    fprintf('Total CVX Solve Time for %1i particles: %1.4f seconds\n'...
        ,N,time_to_solve)

    d = full(d);
    disp('------------------------------------')
    fprintf('Total Run Time: %1.4f seconds\n', cvx_cputime + time_to_solve)