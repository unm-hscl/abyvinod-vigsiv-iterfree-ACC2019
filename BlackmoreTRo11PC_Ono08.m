%% Blackmore TRo 2011 Code for Ono 2008
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 10; 

% Number of particles: 

N = 100; 

% Maximum/minimum bound on input: 

ulim = 0.2; 

% max risk: 

Delta = 0.1;

% System matrices: 
    % Sampling time of the discrete system:
        delT = 0.25;
    % A matrix of the 2D double integrator:
        A = [ 1 1; 
              0 1;];
          
    % B matrix of the 2D double integrator:
        B = [0; 0.033]; 
    % Define the matrices of the discrete time system: 
       Ad = [];
       Bd = zeros(size(A,2)*(T+1),size(size(B,1)*T,1)); 
       Gd = zeros(size(A,2)*T,size(A,2)*(T+1));
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = zeros(T,1);
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
       Bp = cell2mat(Z);
       Bd(:,1:size(u,1)) = Bp(:,1:size(u,1));

       % Randomly generate the disturbance vector from the standard normal.
        for i = 1:N
            w(:,i) = mvnrnd(zeros(1,length(A)*(T+1)),eye(length(A)*(T+1)))';
        end
        w(1:2,1:N) = repmat([0;0],1,N); 
