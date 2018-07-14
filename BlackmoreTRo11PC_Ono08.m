%% Blackmore TRo 2011 Code for Ono 2008
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 10; 

% Number of particles: 

N = 10; 

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

       for i = 0:T
           if i == 0
               
                Ab{i+1}= zeros(size(A,2),size(B,2)); 
               
           else
               
                Ab{i+1}=A^(i-1)*B; 
                
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
        w(1:2,i) = [0; 0];
    end
        
    
% Randomly generate the initial conditions: 
    for i = 1:N
        x0(1:2,i) = [0,0];
    end
    
% Define objects:
    
     ob_a(:,:,1) = [ 1 0;
             -1 0 ];
    ob_b(:,1) =  [0;-0.02];
    
       
cvx_clear
tstart = tic;
cvx_precision best
cvx_begin
    variable u(size(B,2)*T)
    variable x(size(A,2)*(T+1),N)
    variable z(N,1) binary
    variable g(N,size(ob_a,3)) binary
    variable e(T,N,size(ob_a,3)) binary
    variable d(T,size(ob_a,1),N,size(ob_a,3)) binary
    
    minimize(sum(abs(u)))
    
    subject to
        
          for i = 1:N
              
             x(:,i) == Ad*x0(1:2,i)+Bd*u+w(:,i);
              
          end
          
%      abs(u) <= 0.02;
          
    for i = 1:N
         for k = 1:size(ob_a,3)
            for l = 1:size(ob_a,1)  
               for j = 2:T
                  -ob_a(l,:,k)*x((2*(j-1))+1:2*j,i) + ob_b(l,k) <= 500*(1-d(j-1,l,i,k));
                  ob_a(l,:,k)*x((2*(j-1))+1:2*j,i) - ob_b(l,k) <= 500*(d(j-1,l,i,k));

               end
            end
           
         end
        
    end
     
    for i = 1:N
         for k = 1:size(ob_a,3) 
               for j = 1:T
                   
                  -sum(d(j,:,i,k)) + size(ob_a,2) <= 100*(e(j,i,k));
                   sum(d(j,:,i,k)) - size(ob_a,2)+1 <= 100*(1-e(j,i,k));
                   
               end
         end
    end
    
    for i = 1:N
         for k = 1:size(ob_a,3)
               sum(e(:,i,k))  <=  100*(1-g(i,k));
              -sum(e(:,i,k))+1  <=  100*(g(i,k));
         end
    end
    
    
    for i = 1:N
        
            sum(g(i,:)) <= 100*(1-z(i));
           -sum(g(i,:))+1 <= 100*(z(i));
        
    end
      
%              1/N*sum(z)<=Delta;
            
t1 = toc(tstart);
cvx_end;
t2 = toc(tstart);
time_to_solve = t2 - t1;
