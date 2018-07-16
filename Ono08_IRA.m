%% Ono_IRA 2008 code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 10; 

% Number of particles: 

N = 5; 

% Maximum/minimum bound on input: 

ulim = 0.2; 

% max risk: 

Delta = 0.05;

% System matrices: 
    % Sampling time of the discrete system:
        delT = 0.25;
    % A matrix of the 2D double integrator:
        A = [ 1 1; 
              0 1;];
          
    % B matrix of the 2D double integrator:
        B = [0; 0.033]; 
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = ones(size(B,1)*T,1);
       w = zeros(size(x,1),N);
       
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
               
                Ab{i}=A^(i-1)*B; 
                
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
    cov_mat_diag = diag([0.001 0;]); 
    cov_mat = kron(eye(T+1),cov_mat_diag); 
    w = mvnrnd(zeros(size(A,2)*(T+1),1),cov_mat)';
    
% Initial conditions: 
    x0 = [0.01;0];
    
    
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
    cov_X_sans_input = 0 + cov_concat_disturb;

% Generate obstacles: 
    h = [1 0; -1 0;];
    hbig = kron(eye(T+1),h);
    g = [1 0.5];
    gbig = repmat(g',T+1,1);
% Generate intial delta: 
    del = Delta/((T+1)*N); 
    delta = repmat(del,size(h,1)*(T+1),1);
    Jp =0.1;
    epsil = 1e-6;

%  while abs(sum(abs(u))-Jp)<epsil
    Jp = sum(abs(u));
    cvx_clear
    cvx_begin
        variable u(size(B,2)*T)
        variable mean_X_w_input(size(A,2)*(T+1))

        minimize(sum(abs(u)))

        subject to

            mean_X_w_input == mean_X_sans_input + Bd*u; 
            abs(u) <= ulim; 
            hbig*mean_X_w_input<=gbig-sqrt(hbig'*cov_X_sans_input*hbig)*norminv(1-delta);
    cvx_end
%  end
    