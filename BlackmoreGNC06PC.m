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
        G = [1 1 1 1];

% Number of particles to generate. 
    N = 100; 
    
% Time Horizon.
    T = 30;
    
% Randomly generate the disturbance vector from the standard normal.
    for l = 1:N
        W(:,:,l) = mvnrnd(zeros(1,length(A)*T),eye(length(A)*T))';
    end
%     plot(W(:),'+')

% Consider the initial state is known:
    x0 = [0,0,25,0]';

% Generate future state trajectories now we have sampled the disturbance.
    
    % Preallocate concatnated xP vector (pre-input applied): 
        Xp = zeros(length(W),1,N);
        Xp(1:4,:,:) = repmat(x0,1,1,N);
    
    % Generate state tracjectories: 
        for i = 1:N
            for j = 1:T
                if j == 1
                    Xp((4*(j)+1):4*(j+1),:,i) = A*Xp(1:4,:,i)+G*W(1:4,:,i);
                elseif j == T
                    break;
                else
                    Xp((4*(j)+1):4*(j+1),:,i) = A*Xp(4*(j-1):4*j-1,:,i)+G*W(4*(j-1):4*j-1,:,i);
                end
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
    
    % Evaluate binary variable d which represents, given an object O1,
    % it indicates whether the linear constraints that represent the object
    % have been crossed by the particles. 
    
    d = zeros(T,size(ob_a,1),N);
    
    for i = 1:N
        for j = 1:T
            for k = 1:size(ob_a,3)
                for l = 1:size(ob_a,1)
                    if j == 1
                            d(1,l,i,k) = ob_a(l,:,k)*Xp(1:4,:,i)>=ob_b(l,k);
                        else
                            d(j,l,i,k) = ob_a(l,:,k)*Xp((4*(j-1))+1:4*j,:,i)>=ob_b(l,k);
                    end
                end
            end
        end
    end
    
    % Sanity check done to determine if d vector is generated properly:
    % NOTE: Please comment out if you don't need to check. 
    
        P1 = Polyhedron('V', [50, 50; 50, 100; 100, 100; 100, 50;]);
        P2 = Polyhedron('V', [10, 120; 10, 140; 40, 140; 40, 120;]);
        P1.plot()
        hold on
        P2.plot()
        plot(Xp(1:4:T*4,:,2),Xp(3:4:T*4,:,2),'-k+');
        axis([-100 250 -100 250])
    
    % Evaluate binary variable e which indicates if an obstacle has been
    % avoided at time t.
     for i = 1:N
        for j = 1:T
            for k = 1:size(ob_a,3)
                if sum(d(j,:,i,k)) == size(d(j,:,i,k),2)
                    e(j,i,k) = 1; 
                else
                    e(j,i,k) = 0;
                end
            end
            
        end
     end
    
     % Evaluate binary variable g which indicates an obstacle is avoided
     % for all time: 
     for i = 1:N
         for k = 1:1:size(ob_a,3)
            if sum(e(:,i,k)) > 0
                g(i,k) = 1; 
            else
                g(i,k) = 0;
            end
         end
     end 
     
     % Evaluate binary variable z which indicates all obstacles are avoided
     % for all time steps by particle i: 
     for i = 1:N
         for k = 1:1:size(ob_a,3)
            
         end
     end  

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T. 
    

