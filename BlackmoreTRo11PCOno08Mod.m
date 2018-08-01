%% Blackmore TRo 2011 Code with Ono08 Dynamics & Formulation
% Coder: Vignesh Sivaramakrishnan

% Time Horizon: 

T = 30; 

% Number of particles: 

N = 100;


% Maximum/minimum bound on input: 

ulim = 1; 

% System matrices: 

    [Ad,Bd] = doubIntModel(T,delT);
       
% Randomly generate the disturbance vector from the standard normal.
    cov_mat_diag = 0.001*diag([1 0;]); 
    cov_mat = kron(eye(T+1),cov_mat_diag); 
    
% Initial conditions: 
    x0 = [0.4;0];
    xtarget = linspace(-0.495,0,T)';
    xtargetbig = repmat(xtarget,1,N);

% Randomly generate the disturbance vector from the standard normal.
    for i = 1:N
        w(:,i) = mvnrnd(zeros(1,size(Ad,1)),cov_mat)';
    end

% Generate bounds: 
    h = [-1 1];
    htemp = zeros(size(Ad,2)*T);
    
    for k = 1:(size(Ad,2)*T)
        htemp(k,2*(k-1)+1:2*k) = [1 0];   
    end
    htemp = htemp(1:(size(Ad,2)*T)/2,1:(size(Ad,2)*T));
    
    hbig(:,:,1) = kron(eye(T),[h(1,:); 0,0;]);
    g = linspace(0.5,0.1, T);
    gbig = kron(g,[1,1])';
    
    
    input_state_ratio = 0.0001;
    

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
    tstart = tic;
    cvx_clear
        cvx_precision BEST
    cvx_begin quiet
        variable U_vector(size(Bd,2),1);
        variable x(size(Ad,2)*T,N);
    %     variable d(T,size(h,2),N) binary
        variable d(N) binary
    %     variable d(size(h,2),N) binary

        minimize (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
            sum(sum(abs(x(1:2:end,1:end)-xtargetbig)))/(2*g(1)*T));


        subject to


          x(1:end,1:N) == w(3:end,1:N)+repmat(Ad(3:end,:)*x0+...
                          Bd(3:end,:)*U_vector,1,N);

          abs(U_vector) <= ulim;

          for i = 1:N
    %          -kron(htemp,h(1))*x(:,i) + g(:) <= 500*(1-d(:,1,i));
    %           kron(htemp,h(1))*x(:,i) - g(:) <= 500*(d(:,1,i));
    %          -kron(htemp,h(2))*x(:,i) + g(:) <= 500*(1-d(:,2,i));
    %           kron(htemp,h(2))*x(:,i) - g(:) <= 500*(d(:,2,i));

              kron(htemp,h(1))*x(:,i) - g(:) <= 500*(d(i));
              kron(htemp,h(2))*x(:,i) - g(:) <= 500*(d(i));


    %            -kron(htemp,h(1))*x(:,i) + g(:) <= 500*(1-d(1,i));
    %             kron(htemp,h(1))*x(:,i) - g(:) <= 500*(d(1,i));
    %            -kron(htemp,h(2))*x(:,i) + g(:) <= 500*(1-d(2,i));
    %             kron(htemp,h(2))*x(:,i) - g(:) <= 500*(d(2,i));

          end

    %           kron(eye(N),kron(htemp,h(1)))*x(:,1:N) - kron(eye(N),g(:)) <= 500*(d(1:N));
    %           kron(kron(htemp,h(2)),N,1)*x(:,1:N) - kron(eye(N),g(:)) <= 500*(d(1:N));

    %       1/N*sum(pos(sum(sum(d,2),1)))<=0.05;
            1/N*sum(d)<=0.05;
    %       1/N*sum(pos(sum(sum(d,2),1)))<=0.05;

    t1 = toc(tstart);
    cvx_end;
    t2 = toc(tstart);
    time_to_solve = t2 - t1;
    disp('---------BlackmorePC11-----------')
    disp(' ')
    fprintf('Total CVX Run Time for %1i particles: %1.4f seconds\n',...
        N,cvx_cputime)
    disp('------------------------------------')
    fprintf('Total CVX Solve Time for %1i particles:: %1.4f seconds\n'...
        ,N,time_to_solve)

    d = full(d);
