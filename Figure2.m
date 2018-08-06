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
    
    maxlierror = 1E-3;
    lbdelta = 1E-4;
    [onopwl_invcdf_approx_m, onopwl_invcdf_approx_c,...
    onopwl_lb_deltai] = RolleLerpClosedForm(Delta,lbdelta,maxlierror);
    
    T_array = 20:10:80; 
    
    for i = 1:length(T_array)
        T = T_array(i);
        % Initial conditions:    
        x0 = [0.4;0];
%         xtarget = linspace(-0.4,-0.2,T)'; 
xtarget = linspace(-0.5,-0.2,T)';
    %% Bounds on the safe set
    h = [-1 0; 1 0;];
%     g = linspace(0.5,0.5, T);
%     g = linspace(0.5,0.3, T);
 g = linspace(0.5,0.3, T);
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
    plot(T_array,Ono08TTS); hold on; plot(T_array,PWLRA)
    legend('OnoCDC2008','PWL OnoCDC2008')
    xlabel('Time Horizon (seconds)');
    ylabel('Time to compute (CVX) (seconds)');