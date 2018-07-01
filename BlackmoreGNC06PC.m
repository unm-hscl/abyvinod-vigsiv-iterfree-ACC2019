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
    for k = 1:N
        W(:,:,k) = mvnrnd(zeros(1,length(A)*T),eye(length(A)*T))';
    end
%     plot(W(:),'+')

% Consider the initial state is known:
    x0 = [1,0,2,0]';

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
                    break
                else
                    Xp((4*(j)+1):4*(j+1),:,i) = A*Xp(4*(j-1):4*j-1,:,i)+G*W(4*(j-1):4*j-1,:,i);
                end
            end
        end
        

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T. 
    

