%% Blackmore TRo 2011 Code for Ono 2008
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
    for i = 1:N
        w(:,i) = mvnrnd(zeros(size(A,2)*(T+1),1),cov_mat)';
    end
    
% Randomly generate the initial conditions: 
    for i = 1:N
        x0(1:2,i) = [0,0];
    end
    
% Randomly generate obstacles: 
    h = [1 0; -1 0;];
    h = kron(eye(T+1),h);
    g = [0.01 -0.02];
    g = repmat(g',T+1,1);
    
       
cvx_clear
cvx_begin
    variable u(size(B,2)*T)
    variable x(size(A,2)*(T+1),N)
    variable z(1,N) integer
    
    minimize(sum(abs(u)))
    
    subject to
        
          for i = 1:N
              
             x(:,i) == Ad*x0(1:2,i)+Bd*u+w(:,i);
              
          end
          
          u <= ulim;
          u >= 0.01;
          
      for i = 1:N
          for j = 1:T
                g - h*x(:,i) <= 10*(1-z(i))-1e-6;
%                 h*x(:,i)-g <= 10*z(i);
          end
      end
      
%          1/N*sum(sum(z)) <= Delta; 
          
          
cvx_end
