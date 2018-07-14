%% Blackmore TRo 2011 Code for Ono 2008
% Coder: Vignesh Sivaramakrishnan
% 
clc, clear, close all

% Time Horizon: 

T = 10; 

% Number of particles: 

N = 100; 

% Maximum/minimum bound on input: 

ulim = 0.2; 

% max risk: 

Delta = 0.1;

% System matrices: 
    % Sampling time of the discrete system:
        delT = 0.25;
    % A matrix of the 2D double integrator:
        A = [ 1 1; 
              0 1;];
          
    % B matrix of the 2D double integrator:
        B = [0; 0.033]; 
    % Define the matrices of the discrete time system: 
       Ad = [];
       Bd = zeros(size(A,2)*(T+1),size(u,1)); 
       Gd = zeros(size(A,2)*T,size(w,1));
        
    % Preallocate the trajectory vectors to be used: 
       x = zeros(size(A,2)*(T+1),1); 
       u = ones(size(B,1)*T,1);
       w = zeros(size(x,1),N);
       
