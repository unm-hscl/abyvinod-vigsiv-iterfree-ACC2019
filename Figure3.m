%% abyvnod-vigsiv-iterfree-ACC2019: Figure 3
    % This code runs comparisons of OnoCDC08, PWLRealizationOno, and
    % BlackmoreTRo11 for a single tracjetory for a satellite rendevous
    % with a window closing in time and compares how the methods fair.
    %
    % NOTE: The CWH dynamics and portions of the problem setup are from the
    % SReachTools toolbox. 
    % computationally. 
    %
    % REQUIRED DEPENDENCIES: - SReachTools 
    %                        - MATLAB symbolic toolbox
    %                        - MATLAB global optimization toolbox
    
    %% Housekeeping: 
    
        clear;
        close all;
        clc;
        cvx_clear;
        cd('SReachTools-private');
        srtinit
        cd('../');
    
    %% CWH system params
        ulim=0.1;
        mean_disturbance = zeros(4,1);
        covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
        % Define the CWH (planar) dynamics of the deputy spacecraft relative to the chief spacecraft as a LtiSystem object
        sys = getCwhLtiSystem(4,...
                              Polyhedron('lb', -ulim*ones(2,1),...
                                         'ub',  ulim*ones(2,1)),...
                              StochasticDisturbance('Gaussian',...
                                                    mean_disturbance,...
                                                    covariance_disturbance));
                                                
                                                
        T=5; % Stay within a line of sight cone for 4 time steps and 
                        % reach the target at t=5% Safe Set --- LoS cone
                        
        N = 5;
        
    %% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and |vy|<=vymax
        ymax=2;
        vxmax=0.5;
        vymax=0.5;
        hp = [1, 1, 0, 0;           
                     -1, 1, 0, 0; 
                      0, -1, 0, 0;
                      0, 0, 1,0;
                      0, 0,-1,0;
                      0, 0, 0,1;
                      0, 0, 0,-1];
        gbflag = 0;
        gbp = [0;
              0;
              ymax;
              vxmax;
              vxmax;
              vymax;
              vymax];
        safe_set = Polyhedron(hp, gbp);
        %% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
        target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                                'ub', [0.1; 0; 0.01; 0.01]);
        target_tube = TargetTube('reach-avoid',safe_set, target_set, T);                    

        x0 = [-1;         % Initial x relative position
              -1;         % Initial y relative position
               0;            % Initial x relative velocity
               0];           % Initial y relative velocity
        xtarget = repmat([0; 0; 0; 0;],T,1);
       [h, gb] = target_tube.concat();
        h = h(5:end,5:end);
        gb = gb(5:end,1);
        slice_at_vx_vy = x0(3:4);             

        %% Generate matrices for optimal mean trajectory generation
            [~, mean_X_sans_input, cov_X_sans_input] =...
                getHmatMeanCovForXSansInput(sys, x0, T);
            [Ad, Bd, Gd] = getConcatMats(sys, T);
        
        % Probability of being outside the safe set: 

            Delta = 0.2;
            
        % Generate the PW realization of the distribution for PWLRA: 
            maxlierror=1e-2;
            g = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
            fun_monotone = 'mono-inc';
            lower_bound = 1E-5;
            upper_bound = Delta; 
            function_handle = @(z) -g(z);

            [PWA_overapprox_m, PWA_overapprox_c]...
                = getPWAOverAndUnderApprox(lower_bound,upper_bound,...
                maxlierror,function_handle,fun_monotone);
            
    %% Run the following scripts (which should be in the same directory) 
    %% with parameters above: 

        try
            CWHOno08_IRA
        catch
            disp('Ono''s method failed');
        end
        disp(' ');
        disp(' ');
        CWHPiecewiseLinearRA
        disp(' ');
        disp(' ');
        CWHBlackmoreTRo11


    %% Plotting trajectories of each method: 

        plot_markersize = 9;
        plot_fontSize = 9;
        fig3 = figure(1);
        clf
        hold on
        h1 = plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        h2 = plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        h3 = scatter(x0(1),x0(2),10*plot_markersize,'go','filled');
        h4 = scatter([x0(1); ono_opt_mean_X(1:4:end)],...
                    [x0(2); ono_opt_mean_X(2:4:end)],...
                    30, 'bx');
        h5 = scatter([x0(1); onopwl_opt_mean_X(1:4:end)],...
                    [x0(2); onopwl_opt_mean_X(2:4:end)],...
                    30, 'md');
        h6 = scatter([x0(1); blackmore_opt_mean_X(1:4:end)],...
                    [x0(2); blackmore_opt_mean_X(2:4:end)],...
                    30, 'ks');

        
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        
        xlabel('\textbf{Position (y)}')
        ylabel('\textbf{Position (x)}')
        title('\textbf{Trajectory}')
        set(gca,'FontSize',plot_fontSize)
        legend([h1 h2 h3 h4 h5 h6],{'Safe Set',...
            'Target Set','Initial Position',...
            sprintf('Ono2008 IRA Method (Cost: %1.3f)',...
            ono_opt_val), sprintf('Piecewise linear approach (Cost: %1.3f)',...
            onopwl_opt_val),...
            sprintf('Blackmore11 PC Method, %i Particles (Cost: %1.3f)',...
            N,blackmore_opt_val)},'Location','SouthOutside','FontSize',9);
        box on;
        set(fig3,'PaperUnits','centimeters');
        set(fig3,'PaperPosition',[0 0 8.8 8.8]);
        grid on
        axis tight
        fig3 = tightfig(fig3);
        hgexport(fig3,'Figure3',hgexport('factorystyle'),'Format', 'png')
    

    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5;  
        xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
        collection_of_input_vectors = [blackmore_opt_input_vector, ono_opt_input_vector, onopwl_opt_input_vector];
        collection_of_costs = [blackmore_opt_val, ono_opt_val, onopwl_opt_val];
        % Have a collection of string
        collection_of_method_strings = {'BlackmoreTRO11',...
                                        'OnoCDC2008    ',...
                                        'PWL OnoCDC2008'};
    % SReachTools for Monte-Carlo simulation
        max_rel_abserror = 0.1;
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',Delta,max_rel_abserror);
        for input_vec_indx = 1:3
            U_vector = collection_of_input_vectors(:,input_vec_indx);
            method_str = collection_of_method_strings{input_vec_indx};
            true_cost = collection_of_costs(input_vec_indx);
            % This function returns the concatenated state vector stacked columnwise
            X_mcarlo_sans_init_state = generateMonteCarloSims(...
                    n_mcarlo_sims,...
                    sys,...
                    x0,...
                    T,...
                    U_vector);
        % all does it column-wise
        particlewise_result = all(h*X_mcarlo_sans_init_state <= gb);
        prob_estim = sum(particlewise_result)/n_mcarlo_sims;
        cost_estim = mean(sum((X_mcarlo_sans_init_state-xtarget_mcarlo).^2));    
        relative_abserror_in_cost = abs(cost_estim - true_cost)/true_cost;
        if prob_estim >= 1-Delta && relative_abserror_in_cost <= max_rel_abserror
            fprintf('PASSD: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                    method_str,... 
                    n_mcarlo_sims,...
                    prob_estim,...
                    relative_abserror_in_cost);
        else
            fprintf('ERROR: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                    method_str,... 
                    n_mcarlo_sims,...
                    prob_estim,...
                    relative_abserror_in_cost);
        end
    end
