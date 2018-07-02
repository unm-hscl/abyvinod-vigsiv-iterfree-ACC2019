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
    N = 50; 
    
% Time Horizon.
    T = 10;
    
% Desired target trajectory. 
    xref = [200 0 200 0]';
    
% Randomly generate the disturbance vector from the standard normal.
    for l = 1:N
        W(:,:,l) = mvnrnd(zeros(1,length(A)*T),eye(length(A)*T))';
    end
%     plot(W(:),'+')

% Consider the initial state is known:
    x0 = [0,0,25,0]';

% Generate future state trajectories now we have sampled the disturbance.
    
    % Preallocate concatnated xP vector (pre-input applied): 
        Xp = zeros(length(A)*T,1,N);
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
        for i=1:N/16
            plot(Xp(1:4:T*4,:,i),Xp(3:4:T*4,:,i),'-+');
        end
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
     
     g = zeros(N,size(ob_a,3));
     
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
     
     z = zeros(N,1);
     
     for i = 1:N
            if sum(g(i,:)) > 0
                z(i) = 1; 
            else
                z(i) = 0;
            end
     end
     
     % Generate the chance constraint to be used in the optimization
     % problem: 
     
     nu = 1/N*sum(z);
     
% Approximate the expected cost function in terms of particles.

    % Generate desired trajectory for the entire time horizon: 
        xref = repmat(xref,T,1);
        Q = 50*eye(length(A));
        R = 0.001*eye(2*T);
        h = zeros(T,1);

    % The expected cost is modified from Vitus et al. (2012): 
        for i = 1:N
            for j = 1:T
                if j == 1
                   h(j,i) = (Xp(1:4,:,i)-xref(1:4))'*Q*(Xp(1:4,:,i)-xref(1:4));
                else
                   h(j,i) = (Xp((4*(j-1))+1:4*j,:,i)-xref((4*(j-1))+1:4*j))'*Q*(Xp((4*(j-1))+1:4*j,:,i)-xref((4*(j-1))+1:4*j));
                end
                Eh(j,1) = 1/N*sum(h(:,i));
            end
            
        end
        
        
        
 
    

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.

cvx_begin
    variable U(size(B,2)*T)
    variable Xr(size(A,2)*T,1,N)
    variable d(T,size(ob_a,1),N,size(ob_a,1)) binary
    variable e(T,N,size(ob_a,3)) binary
    variable g(N,size(ob_a,3)) binary
    variable z(N) binary
    variable Eh(T,1)
    
    
    minimize( sum( + N*quad_form(U,R)))
    subject to
    
        % Generate state tracjectories: 
        for i = 1:N
            for j = 1:T
                if j == 1
                    Xr((4*(j)+1):4*(j+1),:,i) == A*Xp(1:4,:,i)+G*W(1:4,:,i);
                elseif j == T
                    break;
                else
                    Xr((4*(j)+1):4*(j+1),:,i) == A*Xr(4*(j-1):4*j-1,:,i)+G*W(4*(j-1):4*j-1,:,i)+B*U(2*(j-1):2*j-1,:);
                end
            end
        end
    
    for i = 1:N
        for j = 1:T
            for k = 1:size(ob_a,3)
                for l = 1:size(ob_a,1)
                    if j == 1
                            2000*d(1,l,i,k) + ob_a(l,:,k)*Xr(1:4,:,i)-ob_b(l,k)>=0;
                        else
                            2000*d(j,l,i,k) + ob_a(l,:,k)*Xr((4*(j-1))+1:4*j,:,i)-ob_b(l,k)>=0;
                    end
                end
            end
        end
    end
    
    
    % Evaluate binary variable e which indicates if an obstacle has been
    % avoided at time t.
     for i = 1:N
        for j = 1:T
            for k = 1:size(ob_a,3)
                sum(d(j,:,i,k)) - (size(ob_a,3)-1)<=2000*e(j,i,k);
            end
            
        end
     end
    
     % Evaluate binary variable g which indicates an obstacle is avoided
     % for all time: 
     
     
     for i = 1:N
         for k = 1:1:size(ob_a,3)
            sum(e(:,i,k)) <= 2000*g(i,k);
         end
     end 
     
     % Evaluate binary variable z which indicates all obstacles are avoided
     % for all time steps by particle i: 
     
     for i = 1:N
        sum(g(i,:)) <= 2000*z(i);
     end
         

cvx_end

