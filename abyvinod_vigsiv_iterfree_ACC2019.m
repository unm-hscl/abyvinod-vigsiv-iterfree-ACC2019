%% abyvnod-vigsiv-iterfree-ACC2019
% This code runs comparisons of the 


clc, clear, close all, cvx_clear

% Load parameters: 

% Time Horizon: 

    T = 30; 
    
% Probability of being outside the safe set: 

    Delta = 0.05;

% Initial conditions:   

    x0 = [0.4;0];
    xtarget = linspace(-0.495,0,T)';   
    
% Disturbance parameters: 

    cov_mat_diag = 0.001*diag([1 0;]); 
    mean_w = [0;0];
    
    
% Maximum/minimum bound on input: 

    ulim = 1; 
    
% Sampling time of the discrete system:

    delT = 0.25;

% Number of particles: 

    N = 100;

Ono08_IRA
BlackmoreTRo11PCOno08Mod
%% Plotting
    figure(1)
    hold on
    % h1 = plot(1:(T+1), [-gbig(1);-gbig(1:2:end)],'r');
    % plot(1:(T+1), [gbig(2);gbig(2:2:end)],'r');
    polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
        [gbig(1),gbig(1:2:end)']]';
    P = Polyhedron('V',polyvertex);
    h1 = plot(P,'alpha',0.1);
    % P1 = polytope(hbig,gbig);
    h2 = plot(2:(T+1), xtarget,'go');
    h3 = plot(1:(T+1),mean_X(1:2:end),'-b','LineWidth',1);
    h4 = plot(1:T,x(1:2:end,1),'.');
    plot(1:T,x(1:2:end,2:N),'.')  
    xlabel('time')
    ylabel('x')
    title('Trajectory')
    legend([h1 h2 h3 h4],{'Safe Set','Target Trajectory','Ono2008 IRA Method','Blackmore11 PC Method'});

    figure(2);
    hold on
    plot(opt_value_array)
    xlabel('iteration count')
    ylabel('cost')
    title('Final J cost for OnoIRA2008');