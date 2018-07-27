%% Blackmore TRo 2011 Code with Ono08 Dynamics & Formulation
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
    w = mvnrnd(
    
% Initial conditions: 
    x0 = [0.4;0];
    xtarget = linspace(-0.495,0,T)';

% Generate bounds: 
    h = [-1 0; 1 0;];
    hbig = kron(eye(T),h);
    g = linspace(0.5,0.1, T);
    gbig = kron(g,[1,1])';