%% abyvnod-vigsiv-iterfree-ACC2019
    % This code runs comparisons of the 


    clear;
    close all;
    clc;
    cvx_clear;
    
    cd('SReachTools-private');
    srtinit
    cd('../');
    % Load parameters: 

    % Time Horizon: 

        T = 30; 

    % Probability of being outside the safe set: 

        Delta = 0.05;
%       Delta = 0.6;

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

    % Number of particles for BlackmorePCApproach: 

        N = 10;
    
    % Run the following scripts (which should be in the same directory) 
    % with parameters above: 

        Ono08_IRA
        disp(' ');
        disp(' ');
        PiecewiseLinearRA
        disp(' ');
        disp(' ');
        BlackmoreTRo11PCOno08Mod


%% Plotting
    plot_markersize = 15;
    figure(1)
    clf
    hold on
    polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
        [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
    P = Polyhedron('V',polyvertex);
    h1 = plot(P,'alpha',0.1);
    h2 = plot(2:(T+1), xtarget,'go','MarkerSize',plot_markersize);
    h3 = plot(1:(T+1),ono_opt_mean_X(1:2:end),'bx','LineWidth',1,'MarkerSize',plot_markersize);
    h4 = plot(1:T,blackmore_opt_mean_X(1:2:end,1),'ks','MarkerSize',plot_markersize);
    h5 = plot(1:(T+1),onopwl_opt_mean_X(1:2:end),'md','LineWidth',1,'MarkerSize',plot_markersize);
%   h4 = plot(1:T,x(1:2:end,1),'.');
%     plot(1:T,x(1:2:end,2:N),'.')  
    xlabel('time')
    ylabel('x')
    title('Trajectory')
    legend([h1 h2 h3 h4 h5],{'Safe Set','Target Trajectory',...
        sprintf('Ono2008 IRA Method (Cost: %1.3f)',ono_opt_val),sprintf('Blackmore11 PC Method (Cost: %1.3f)',blackmore_opt_val),sprintf('Piecewise linear approach (Cost: %1.3f)',onopwl_opt_val)});
    box on;
    set(gca,'FontSize',20)
    
    % Figure 2 shows the cost of Ono's method with each iteration.
    figure(2);  
    clf
    hold on
    plot(ono_opt_value_array)
    xlabel('Iteration Count')
    ylabel('Cost, J')
    title('Final J cost for OnoIRA2008');
    
    figure(1);
    

% %% Monte-Carlo simulation
% n_mcarlo_sims = 1e5;
% collection_of_input_vectors = [blackmore_opt_input_vector, ono_opt_input_vector, onopwl_opt_input_vector];
% % Have a collection of string
% 
% % SReachTools for Monte-Carlo simulation
% sys=getChainOfIntegLtiSystem(2, delT, Polyhedron('lb',-ulim,'lb',ulim),RandomVector('Gaussian',mean_w,cov_mat_diag));    
% for input_vec_indx = 1:3
%     U_vector = collection_of_input_vectors(:,input_vec_indx);
%     % This function returns the concatenated state vector stacked columnwise
%     x = generateMonteCarloSims(...
%             n_mcarlo_sims,...
%             sys,...
%             x0,...
%             T,...
%             U_vector);
%     % all does it column-wise
%     particlewise_result = all(hbig*x <= gbig);
%     prob_estim = sum(particlewise_result)/n_mcarlo_sims;
%     cost_estim = (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
%             sum(sum(abs(x(1:2:end,1:end)-xtargetbig)))/(2*g(1)*T)/n_mcarlo_sims);
%     % Add a %s before
%     fprintf('Monte-Carlo simulation using %1.0e particles | Prob. of Hx<=g = %1.3f | Cost = %1.3f\n',...
%             n_mcarlo_sims,...
%             prob_estim,...
%             cost_estim);
% end