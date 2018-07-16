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

Delta = 0.1;

% System matrices: 
    % A matrix:
        A = [ 1 1; 
              0 1;];
          
    % B matrix:
        B = [0; 0.033]; 
        G = eye(4);
    % Define the matrices of the discrete time system: 
       Ad = [];
       Bd = zeros(size(A,2)*(T+1),size(size(B,1)*T,1)); 
       Gd = zeros(size(A,2)*T,size(A,2)*(T+1));
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = zeros(size(B,2)*T,1);
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
cov_mat = [0.001 0; 0 0;]; 
cov_mat_huge = kron(eye((T+1)),cov_mat); 
for i = 1:N
    w(:,i) = mvnrnd(zeros(1,length(A)*(T+1)),cov_mat_huge)';
end
w(1:2,1:N) = repmat([0;0],1,N);

% Generate bounds randomly. 

g(1) = normrnd(-0.01,0.01);
g(2) = -normrnd(0.01,0.01);
gg = kron(ones((T-1),1),g');
g_huge = kron(ones(N*(T-1),1),g'); 

% Generate h matrix for all particles for all time. 

h(1,:) = [1  0];
h(2,:) = [-1 0];
hh = kron(eye((T-1)),h); 
h_huge = kron(eye(N*(T-1)),h); 
% Generate initial conditions randomly: 

x0 = mvnrnd([0.01; 0],cov_mat)';


% Begin optimization problem: 
cvx_clear
cvx_begin

    variable u(size(B,2)*T)
    variable x(size(A,2)*(T+1),N)
    variable z(N,1) binary
    
    minimize(sum(abs(u)))
    
    subject to
    
          for i = 1:N
              
             x(:,i) == Ad*x0+Bd*u+w(:,i);
              
          end
          
          abs(u) <= ulim;
          
       for i = 1:N
          hh*x(size(A,2)+1:end-size(A,2),i) - gg  <=  100*(1-z(i));
       end
          
%          1/(N)*sum(z)>=1-Delta;

cvx_end
x = full(x); 
        
        
