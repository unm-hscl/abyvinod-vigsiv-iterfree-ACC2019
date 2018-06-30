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
        G = [1; 1; 1; 1;];

% Number of particles to generate. 
    N = 100; 
    
% Time Horizon.
    T = 30;
    
% Randomly generate the disturbance vector from the standard normal.
    W = mvnrnd(zeros(1,length(A)*N),eye(length(A)*N))';
    plot(W(:),'+')

% Consider the initial state is known:
    x0 = [1,0,2,0]';

% Generate future state trajectories now we have sampled the disturbance.
    
    % Preallocate concatnated x vector: 
        X = zeros(1,length(W))';
        X(1:length(x0)) = x0;
    % Preallocate A vector: 
        Ac = zeros(length(A)*length(X),length(A));
        Ac(1:length(A),1:length(A)) = eye(length(A)); 
        Ac(length(A)+1:2*length(A),1:length(A)) = A; 
        
    
    % Generate state tracjectories: 
        for i = 1:length(X)
            
            if i == 1
                Ac(1:length(A),1:length(A)) = eye(length(A));
                X(i+4:i+8) = Ac*X(i:i+4)+G
            elseif i == 2
                Ac(length(A)+1:2*length(A),1:length(A)) = A;
            end
            
        end
        

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T. 
    

