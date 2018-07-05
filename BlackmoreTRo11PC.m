%% Blackmore TRo 2011 Code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

    T = 20; 
    
% Number of Particles: 

    N = 10;
    
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
       u = zeros(size(B,1)*T,size(B,2));
       w = zeros(size(A,1)*T,1);
       y = zeros(size(C,1)*T,1);
       v = zeros(size(C,1)*T,1);
   % Define the matrices to the discrete time system: 
       Ad = zeros(size(A,2)*T,size(A,2));
       Ad(1:4,1:4) = eye(4); 
       Bd = zeros(size(A,2)*T,size(B,2)*T); 
       Gd = zeros(size(A,2)*T,size(A,2)*T);
       
       
%        New Answer! Simply define a cell vector of matrices, which must all be the same size. The rest happens automagically. Some fake data:
% 
% B{1} = [0,1;2,3];
% B{2} = [4,5;6,7];
% B{3} = [8,9;10,11];
% B{4} = [12,13;14,15];
% Create numeric matrix A with the elements of B on the diagonals:
% 
% N = numel(B);
% X = tril(true(N));
% Y = hankel(ones(1,N));
% [R,C] = find(Y);
% U = accumarray(R,find(X),[],@(n){n});
% Z = cell(N);
% for k = 1:N
%     Z(U{k}) = B(k);
% end
% Z(~X) = {zeros(size(B{1}))};
% A = cell2mat(Z);
       
       
       for i = 2:T
           Ad(4*(i-1)+1:4*i,:) = A^(i-1);
           if i ~= T-1
            Ab(4*(i-1)+1:4*(i),:) = A^(i-1);
           else
               Ab(1:4,:) = ones(4); 
           end
       end
        BB = 
        Bd(2*size(A,2)+1:end,1:size(B,2)) = Ab*B;