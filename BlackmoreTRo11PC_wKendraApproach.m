%% Blackmore TRo 2011 Code
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all
tic
% Time Horizon: 

    T = 30; 
    
% Number of Particles: 

    N = 5;

% Initial positition: 
    
    x0 = [0; 0; 0; 0];
    
% Desired position: 

    xref = [200;0;200;0;];
    xrefh = repmat(xref,T+1,N);
    
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
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = zeros(size(B,2)*T,1);
       w = zeros(size(G,1)*T,1);
       
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
       Bp = cell2mat(Z);
       Bd(:,1:size(u,1)) = Bp(:,1:size(u,1));


       
       for i = 1:(T+1)
           if i == 1
               
                Ga{i}= zeros(size(G,2)); 
           
           else
               
                Ga{i}=A^(i-1);
               
           end
           
       end

       No = numel(Ga);
       X = tril(true(No));
       Y = hankel(ones(No,1));
       [R,C] = find(Y);
       U = accumarray(R,find(X),[],@(n){n});
       Z = cell(No);
       for k = 1:No
           Z(U{k}) = Ga(k);
       end
       Z(~X) = {zeros(size(Ga{1}))};
       Gp = cell2mat(Z);
       Gd = Gp(:,1:size(w,1));
       

       
       
% Randomly generate the disturbance vector from the standard normal.
    for i = 1:N
        w(:,i) = mvnrnd(zeros(1,length(A)*T),eye(length(A)*T))';
    end
    
% Define objects:
    
        ob_a = [ 1 0  0 0;
              0 0  1 0;
             -1 0  0 0;
              0 0 -1 0;
              1 0  0 0;
              0 0  1 0;
             -1 0  0 0;
              0 0 -1 0;];
       ob_a = repmat(ob_a,1,T+1); 

       ob_b =  [50;50;-100;-100;10;120;-40;-140];


% Generate desired trajectory for the entire time horizon: 
    Q = [200 0 0 0; 
         0 100 0 0; 
         0 0 200 0; 
         0 0 0 100;];
    R = 0.001*eye(size(B,2)*T);
    
%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
Qhugep=kron(eye(T+1),Q);
% t = ones(T,size(ob_a,1)*N*size(ob_a,3));

cvx_clear
cvx_begin
    variable u(size(B,2)*T)
    variable x(size(A,2)*(T+1),N)
    variable d(T-1,N) binary
    
    minimize(1/N*trace((sum(x-xrefh,2))'*Qhugep*sum((x-xrefh),2)) + u'*R*u)
%     minimize( sum(z) + 1/N*sum(sum(h)))
    subject to
    
          for i = 1:N
              
             x(:,i) == Ad*x0+Bd*u+Gd*w(:,i);
              
          end

          abs(u) <= 15;
               
      for i = 1:N
           for l = 1:size(ob_a(:,1))
               ob_a(l,:)*x(:,i) - ob_b(l) <= 50000*(1-d(i))-1;
%               -ob_a(l,:)*x(:,i) + ob_b(l) <= 50000*d(i);
           end
      end
%             1/N*sum(sum(d,2))<=0.02;
            
cvx_end
    
    x = full(x);
    figure
    P1 = Polyhedron('V', [50, 50; 50, 100; 100, 100; 100, 50;]);
    P2 = Polyhedron('V', [10, 120; 10, 140; 40, 140; 40, 120;]);
    P1.plot()
    hold on
    P2.plot()
    for i=1:N
        plot(x(1:4:T*4,i),x(3:4:T*4,i),'+');
   end
    axis([-100 400 -100 400])

toc
       