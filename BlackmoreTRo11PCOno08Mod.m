%% Blackmore TRo 2011 Code with Ono08 Dynamics & Formulation
% Coder: Vignesh Sivaramakrishnan

% System matrices: 

    disp('---------BlackmorePC11-----------')
    disp(' ')
    fprintf('No. of particles: %d\n',N);
    
    large_constant = 5000;
% Generate the cov_matrix for optimization problem given covariance.
    cov_mat = kron(eye(T+1),cov_mat_diag); 
    
% Vectorize the target trajectory for the optimization problem. 
    xtargetbig = repmat(xtarget,1,N);

% Randomly generate the disturbance vector from the standard normal.

    mean_w_vec = repmat(mean_w(:),1,T+1);
    
    for i = 1:N
        w(:,i) = mvnrnd(zeros(1,size(Ad,1)),cov_mat)';
    end

% Generate bounds: 
    htemp = zeros(size(Ad,2)*T);
    
    for k = 1:(size(Ad,2)*T)
        htemp(k,2*(k-1)+1:2*k) = [1 0];   
    end
    htemp = htemp(1:(size(Ad,2)*T)/2,1:(size(Ad,2)*T));
    
    
    input_state_ratio = 0.0001;
    

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
    tstart = tic;
    cvx_clear
        cvx_precision BEST
    cvx_begin quiet
        variable U_vector(size(Bd,2),1);
        variable xBl(size(Ad,2)*T,N);
        variable mean_X(size(Ad,2)*T,1);
        variable d(N) binary;

        minimize (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
            sum(sum(abs(xBl(1:2:end,1:end)-xtargetbig)))/(2*g(1)*T)/N);

        subject to
          mean_X == Ad(3:end,:)*x0+ Bd(3:end,:)*U_vector;

          xBl(1:end,1:N) == Gd(3:end,:)*w(3:end,1:N)+repmat(mean_X,1,N);

          abs(U_vector) <= ulim;

          for i = 1:N
              kron(htemp,h(1))*xBl(:,i) - g(:) <= large_constant*(d(i));
              kron(htemp,h(2))*xBl(:,i) - g(:) <= large_constant*(d(i));
          end
          1/N*sum(d)<=Delta;

    t1 = toc(tstart);
    cvx_end;
    t2 = toc(tstart);
    time_to_solve = t2 - t1;
    
    blackmore_opt_mean_X = [x0;mean_X];
    blackmore_opt_val = cvx_optval;
    blackmore_opt_input_vector = U_vector;            
    
    fprintf('Total CVX Run Time for %1i particles: %1.4f seconds\n',...
        N,cvx_cputime)
    disp('------------------------------------')
    fprintf('Total CVX Solve Time for %1i particles: %1.4f seconds\n'...
        ,N,time_to_solve)

    d = full(d);
    disp('------------------------------------')
    fprintf('Total Run Time: %1.4f seconds\n', cvx_cputime + time_to_solve)