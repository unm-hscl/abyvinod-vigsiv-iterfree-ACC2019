%% Blackmore TRo 2011 Code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all
tic
% Time Horizon: 

    T = 20; 
    
% Number of Particles: 

    N = 10;

% Initial positition: 
    
    x0 = [1; 1; 1; 1];
    
% System matrices: 
    % Sampling time of the discrete system:
        delT = 0.25;
    % A matrix of the 2D double integrator:
        A = [ 1 delT 0    0;
              0    1 0    0;
              0    0 1 delT;
              0    0 0    1;];
    % B matrix of the 2D double integrator:
        B = [0.5*delT^2          0; 
                   delT          0; 
                      0 0.5*delT^2;
                      0       delT;];
        G = eye(4);
        C = eye(4); 
        
    % Define the trajectories to be used: 
       x = zeros(size(A,2)*(T+1),1);
       u = zeros((T+1)*size(B,2),1);
       w = zeros(size(A,1)*(T+1),1);
       y = zeros(size(C,1)*T,1);
       v = zeros(size(C,1)*T,1);
       
   % Define the matrices to the discrete time system: 
       Ad = zeros(size(A,2)*T+1,size(A,2));
       Bd = zeros(size(A,2)*T-1,size(B,2)*(T-1)); 
       Gd = zeros(size(A,2)*T-1,size(A,2)*T);
       
       
       for i = 1:(T+1)
           
           Ad(size(A,2)*(i-1)+1:size(A,2)*i,:) = A^(i-1);
           
       end
       
       for i = 0:T
           if i == 0
               
                Ab{i+1}= zeros(size(A,2),size(B,2)); 

           else
               
                Ab{i+1}=A^(i-1)*B; 
                
           end
            
           
           
       end
       
       N = numel(Ab);
       X = tril(true(N));
       Y = hankel(ones(N,1));
       [R,C] = find(Y);
       U = accumarray(R,find(X),[],@(n){n});
       Z = cell(N);
       for k = 1:N
           Z(U{k}) = Ab(k);
       end
       Z(~X) = {zeros(size(Ab{1}))};
       Bd = cell2mat(Z);
       
       
       
       
       for i = 1:T+1
           if i == 1
               
                Ga{i}= zeros(size(A,2)); 
           
           else
               
               Ga{i}=A^(i-2);
               
           end
           
       end

       N = numel(Ga);
       X = tril(true(N));
       Y = hankel(ones(N,1));
       [R,C] = find(Y);
       U = accumarray(R,find(X),[],@(n){n});
       Z = cell(N);
       for k = 1:N
           Z(U{k}) = Ga(k);
       end
       Z(~X) = {zeros(size(Ga{1}))};
       Gd = cell2mat(Z);
       
   x = Ad*x0+Bd*u+Gd*w; 
       
       toc
       