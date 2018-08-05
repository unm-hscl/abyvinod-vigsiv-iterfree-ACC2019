%% abyvnod-vigsiv-iterfree-ACC2019, Figure 2
    % This code runs comparisons of the 


    clear;
    close all;
    clc;
    cvx_clear;
    
    cd('SReachTools-private');
    srtinit
    cd('../');
    % Load parameters: 

    % Probability of being outside the safe set: 

        Delta = 0.05;



        cov_mat_diag = 0.0001*diag([1 0;]); 
        mean_w = [0;0];


    % Maximum/minimum bound on input: 

        ulim = 1; 

    % Sampling time of the discrete system:

        delT = 0.25;

    % Number of particles for BlackmorePCApproach: 

        N = 200;
    % Desired accuracy
    desired_accuracy = 0.001;       
        
    %% Cost ratios b/n input and state --- scalarization term
    input_state_ratio = 0.0001;
    
    T_array = 30; 
    
    for i = 1:length(T_array)
        T = T_array(i);
        % Initial conditions:    
        x0 = [0.4;0];
        xtarget = linspace(-0.2,0,T)'; 
    %% Bounds on the safe set
    h = [-1 0; 1 0;];
%     g = linspace(0.5,0.5, T);
    g = linspace(0.5,0.5, T);
    % Disturbance parameters: 
        
    % Prep
    % System matrices: 

    % Generate a large cov_mat for the optimizaiton problem.
        cov_mat = kron(eye(T+1),cov_mat_diag); 

    % Generate nominal x (Note this is a code snippet taken from SReachTools):
        sys=getChainOfIntegLtiSystem(2, delT, Polyhedron('lb',-ulim,'lb',ulim),RandomVector('Gaussian',mean_w,cov_mat_diag));    
        [Ad, Bd, Gd] = getConcatMats(sys, T);
        [~, mean_X_sans_input, cov_X_sans_input] = getHmatMeanCovForXSansInput(sys, x0, T);      
        
        
    %% Run the following scripts (which should be in the same directory) 
    %% with parameters above: 

        try
            Ono08_IRA
        catch
            disp('Ono''s method failed');
        end
        disp(' ');
        disp(' ');
        PiecewiseLinearRA
        
        
        Ono08TTS(i) = ono_time_to_solve;
        PWLRA(i) = onopwl_time_to_solve;   
    end
    
%% Plotting
%     plot_markersize = 15;
%     plot_fontSize = 20;
    plot_markersize = 10;
    plot_fontSize = 10;
    figure(1)
    clf
    hold on
    polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
        [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
    P = Polyhedron('V',polyvertex);
    h1 = plot(P,'alpha',0.1);
    h11 = scatter(1,x0(1),plot_markersize*10,'filled');
    h2 = plot(2:(T+1),xtarget,'go','MarkerSize',plot_markersize,'LineWidth',2);
    h3 = plot(2:(T+1),ono_opt_mean_X(1:2:end),'bx','LineWidth',1,'MarkerSize',plot_markersize);
    h5 = plot(2:(T+1),onopwl_opt_mean_X(1:2:end),'md','LineWidth',1,'MarkerSize',plot_markersize);    
%   h4 = plot(1:T,x(1:2:end,1),'.');
%     plot(1:T,x(1:2:end,2:N),'.')  
    xlabel('time')
    ylabel('x')
    title('Trajectory')
    legend([h1 h11 h2 h3 h5],{'Safe Set','Initial state','Target Trajectory',...
        sprintf('Ono2008 IRA Method (Cost: %1.3f)',ono_opt_val),sprintf('Piecewise linear approach (Cost: %1.3f)',onopwl_opt_val)});
    box on;
    set(gca,'FontSize',plot_fontSize)