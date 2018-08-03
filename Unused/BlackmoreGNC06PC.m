%% Blackmore TRo2011 Particle Approach Code
% Coder: Vignesh Sivaramakrishnan
% Date: 6/29/2018

% NOTE: If other figures and workspace variables need to be present please
% comment out the following line!
clear, clc, close all

%% Initial Model Specifications: 
% Section gives specifics with regard to disturbances, initial conditions,
% and other important model parameters.

% Initial system matricies.

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
    % Disturbance matrix of the 2D double integrator: 
        G = eye(4);

% Number of particles to generate. 
    N = 50; 
    
% Time Horizon.
    T = 30;
    
% Desired target trajectory. 
    xref = [200 0 200 0]';
    xrefh = repmat(xref,T,N);
    
% Randomly generate the disturbance vector from the standard normal.
    for i = 1:N
        W(:,i) = mvnrnd(zeros(1,length(A)*T),eye(length(A)*T))';
    end
%     plot(W(:),'+')

% Consider the initial state is known:
    x0 = [0,0,0,0]';

% Generate future state trajectories now we have sampled the disturbance.
    
    % Preallocate concatnated xP vector (pre-input applied): 
        Xp = zeros(length(A)*T,N);
        Xp(1:4,:) = repmat(x0,1,N);
    
    % Generate state tracjectories: 
        for i = 1:N
            for j = 1:T-1
                
                Xp((4*(j)+1):4*(j+1),i) = A*Xp(4*(j-1)+1:4*j,i)+ G*W(4*(j-1)+1:4*j,i);

            end
        end
        
% Approximate the chance constraints in terms of generated particles. 
    
    % Define the obstacle(s) for which we will evalutate the binary function
    % d:
    ob_a(:,:,1) = [ 1 0  0 0;
              0 0  1 0;
             -1 0  0 0;
              0 0 -1 0;];
    ob_a(:,:,2) = [ 1 0  0 0;
              0 0  1 0;
             -1 0  0 0;
              0 0 -1 0;];
    ob_b(:,1) =  [50;50;-100;-100];
    ob_b(:,2) =  [10;120;-40;-140];
    
     
% Approximate the expected cost function in terms of particles.

    % Generate desired trajectory for the entire time horizon: 
        Q = [200 0 0 0; 
             0 100 0 0; 
             0 0 200 0; 
             0 0 0 100;];
        R = 0.001*eye(size(B,2)*T);

        
        
        
 
    

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
Qhugep=kron(eye(T),Q);
% t = ones(T,size(ob_a,1)*N*size(ob_a,3));

cvx_clear
cvx_precision best
cvx_begin
    variable U(size(B,2)*T)
    variable Xr(size(A,2)*T,N)
    variable z(N,1) binary
    variable g(N,size(ob_a,3)) binary
    variable e(T,N,size(ob_a,3)) binary
    variable d(T,size(ob_a,1),N,size(ob_a,3)) binary
    
    minimize(sum(z)+1/N*trace((sum(Xr-xrefh,2))'*Qhugep*sum((Xr-xrefh),2)) + U'*R*U)
%     minimize( sum(z) + 1/N*sum(sum(h)))
    subject to
    
        Xr(1:4,:) == Xp(1:4,:);
        % Generate state tracjectories: =
        for i = 1:N
            for j = 1:T-1
%                     abs(U(2*(j-1)+1:2*j))<=[0.5;0.5]; 
                    Xr((4*j+1):4*(j+1),i) == B*U(2*(j-1)+1:2*j) + A*Xr(4*(j-1)+1:4*(j),i)+G*W(4*(j-1)+1:4*j,i);
            end
        end
          U <=  20;
          U >= -20;
          
          
      for i = 1:N
         for k = 1:size(ob_a,3)
            for l = 1:size(ob_a,1)  
               for j = 1:T
                   ob_a(l,:,k)*Xr((4*(j-1))+1:4*j,i) - ob_b(l,k) + 500*(1-d(j,l,i,k)) >= 0;
                  ob_a(l,:,k)*Xr((4*(j-1))+1:4*j,i) - ob_b(l,k) - 500*(d(j,l,i,k)) <= 0;

               end
            end
           
         end
        
      end
     
    for i = 1:N
         for k = 1:size(ob_a,3) 
               for j = 1:T
                   
                  -sum(d(j,:,i,k)) + size(ob_a,2) <= 100*(1-e(j,i,k));
                   sum(d(j,:,i,k)) - size(ob_a,2)+1 <= 100*(e(j,i,k));
                   
               end
         end
    end
    
    for i = 1:N
         for k = 1:size(ob_a,3)
               sum(e(:,i,k))-1  <=  100*(g(i,k))-1;
              -sum(e(:,i,k))+1  <=  100*(1-g(i,k));
         end
    end
    
    
    for i = 1:N
        
            sum(g(i,:))-1 <= 100*(z(i))-1;
           -sum(g(i,:))+1 <= 100*(1-z(i));
        
    end
      
            1/N*sum(z)<=0.02;
            
cvx_end
    
    Xr = full(Xr);
    figure
    P1 = Polyhedron('V', [50, 50; 50, 100; 100, 100; 100, 50;]);
    P2 = Polyhedron('V', [10, 120; 10, 140; 40, 140; 40, 120;]);
    P1.plot()
    hold on
    P2.plot()
    for i=1:N
        plot(Xr(1:4:T*4,i),Xr(3:4:T*4,i),'+');
   end
    axis([-100 400 -100 400])
    
    
% cvx_clear
% cvx_begin
%     variable U(size(B,2)*T)
%     variable Xr(size(A,2)*T,1)
%     
%     minimize( trace((Xr-xrefh)'*Qhuge*(Xr-xrefh)))
% %     minimize( sum(z) + 1/N*sum(sum(h)))
%     subject to
%     
%         Xr(1:4) == Xp(1:4,1);
%         % Generate state tracjectories: =
%             for j = 1:T-1 
%                     Xr((4*j+1):4*(j+1)) == B*U(2*(j-1)+1:2*j) + A*Xr(4*(j-1)+1:4*(j));
%             end
%             abs(U) <= 200;
%             
%             
% cvx_end
%     
%     Xr = full(Xr);
% %     for i=1:N
%         plot(Xr(1:4:T*4),Xr(3:4:T*4),'-+');
% %     end
%     axis([-100 250 -100 250])