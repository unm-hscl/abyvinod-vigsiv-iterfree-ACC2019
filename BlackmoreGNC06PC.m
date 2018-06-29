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
      
          

% % Distribution of initial state vector. Note: Only the x and y positions
% % have uncertainity for this problem. Their velocities start at zero. 
%     d = 100;
%     mu = [0 0];
%     Sigma = [1 0; 0 1];
%     x = linspace(1,3,d); y = linspace(1,2,d);
%     [X1,X2] = meshgrid(x,y);
%     x0Dist = mvnpdf([X1(:) X2(:)],mu,Sigma);
%     x0Dist2 = reshape(x0Dist,length(y),length(x));
%     surf(x,y,x0Dist2);
%     caxis([min(x0Dist2(:))-.5*range(x0Dist2(:)),max(x0Dist2(:))]);
%     axis([1 3 1 2 0 .2])
%     xlabel('x'); ylabel('y'); zlabel('Probability Density');

% Distribution of the disturbance.
    % Distribution of initial state vector. Note: Only the x and y positions
% have uncertainity for this problem. Their velocities start at zero. 
    d = 100;
    mu = [0 0 0 0];
    Sigma = eye(4);
    xerr = linspace(-3,3,d); xdoterr = linspace(-0.5,0.5,d);
    yerr = linspace(-2,2,d); ydoterr = linspace(-0.5,0.5,d);
    [X1,X2] = meshgrid(xerr,yerr);
    [X3,X4] = meshgrid(xdoterr,ydoterr);
    v0Dist = mvnpdf([X1(:) X2(:) X3(:) X4(:)],mu,Sigma);
    v0Dist2 = reshape(v0Dist,length(yerr),length(xerr));
    surf(xerr,yerr,v0Dist2);
    caxis([min(v0Dist2(:))-.5*range(v0Dist2(:)),max(v0Dist2(:))]);
    axis([-3 3 -2 2 0 .2])
    xlabel('x'); ylabel('y'); zlabel('Probability Density');

% Number of particles to generate. 
    N = 100; 
% Time Horizon
    T = 30;

%% Generate samples
% Section generates samples according to BlackmoreGNC06 to be used in the
% optimization problem. 

% % Generate samples of the initial conditions:
%     xsamp = randsample(x0Dist,N);
%     c = ismember(x0Dist2,xsamp);
%     [j,k] = find(c);
%     x0 = zeros(length(A),N);
%     x0(1,:) = randsample(x(j),N); 
%     x0(3,:) = randsample(y(k),N); 
%     plot(x0(1,:),x0(3,:),'k.')
%     axis([-150 150 -150 150])

% Consider the initial state is known:

      x0 = [1,0,2,0]';

% Generate samples of the disturbance:
    v = zeros(length(A),N,T-1);
    for i = 1:T-1
        vsamp = randsample(v0Dist,N);
        c = ismember(v0Dist2,vsamp);
        [j,k] = find(c);      
        v(1,:,i) = randsample(xerr(j),N);
        v(2,:,i) = randsample(xdoterr(j),N);
        v(3,:,i) = randsample(yerr(k),N); 
        v(4,:,i) = randsample(ydoterr(j),N);
    end



%% Run optimization problem for an optimal control policy
% We run a control policy to determine

