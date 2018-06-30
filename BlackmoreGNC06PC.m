%% Blackmore GNC06 Particle Approach Code
% Coder: Vignesh Sivaramakrishnan
% Date: 6/29/2018

% NOTE: If other figures and workspace variables need to be present please
% comment out the following line!
clear, clc, close all

%% Initial Model Specifications: 
% Section gives specifics with regard to disturbances, initial conditions,
% and other important model parameters.

% Initial system matricies.

    % sampling time of the discrete system
    delT = 0.25;
    % A matrix of the 2D double integrator
    A = [ 1 delT 0    0;
          0    1 0    0;
          0    0 1 delT;
          0    0 0    1;];
    % B matrix of the 2D double integrator
    B = [0.5*delT^2          0; 
               delT          0; 
                  0 0.5*delT^2;
                  0       delT;];

% Number of particles to generate. 
    N = 100; 
    
% Time Horizon
    T = 30;
    
% Randomly generate the disturbance vector from the standard normal
    W = mvnrnd(zeros(1,length(A)*N),eye(length(A)*N))';
    plot(W(:),'+')

%% Generate samples
% Section generates samples according to BlackmoreGNC06 to be used in the
% optimization problem. 

% Consider the initial state is known:

    x0 = [1,0,2,0]';

%% Run optimization problem for an optimal control policy
% We run a control policy to determine

