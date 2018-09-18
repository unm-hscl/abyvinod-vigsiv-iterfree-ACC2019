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
        mean_w = zeros(4,1);
        covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
        % Define the CWH (planar) dynamics of the deputy spacecraft relative to the chief spacecraft as a LtiSystem object
        sys = getCwhLtiSystem(4,...
                              Polyhedron('lb', -ulim*ones(2,1),...
                                         'ub',  ulim*ones(2,1)),...
                              StochasticDisturbance('Gaussian',...
                                                    mean_w,...
                                                    covariance_disturbance));
                                                
                                                
        T=5; % Stay within a line of sight cone for 4 time steps and 
                        % reach the target at t=5% Safe Set --- LoS cone
                        
        N = 100;
        Q = kron(eye(T),diag([10 1 10 1])); 
        R = kron(eye(T),0.001*eye(2));
        
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

        x0 = [-1.15;         % Initial x relative position
              -1.15;         % Initial y relative position
               0;            % Initial x relative velocity
               0];           % Initial y relative velocity
        xtarget = repmat([0; 0; 0; 0;],T,1);
       [h, gb] = target_tube.concat();
        hbig = h(8:end,5:end);
        gbig = gb(8:end);
        state_offset = 1;
        slice_at_vx_vy = x0(3:4);             

        %% Generate matrices for optimal mean trajectory generation
            [~, mean_X_sans_input, cov_X_sans_input] =...
                getHmatMeanCovForXSansInput(sys, x0, T);
            [Ad, Bd, Gd] = getConcatMats(sys, T);
        
        % Probability of being outside the safe set: 

            Delta = 0.6;
            
        if Delta <= 0.5
            % Generate the PW realization of the distribution for PWLRA: 
            maxlierror_phiinv=1e-2;
            phiinv = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
            fun_monotone_phiinv = 'mono-inc';
            lower_bound_phiinv = 1E-5;
            upper_bound_phiinv = Delta; 
            function_handle = @(z) -phiinv(z);
            [~,~,PWA_negphiinv_underapprox_m,...
                PWA_negphiinv_underapprox_c] =...
                getPWAOverAndUnderApprox(lower_bound_phiinv,...
                upper_bound_phiinv,maxlierror_phiinv,...
                function_handle,fun_monotone_phiinv);
            PWA_phiinv_overapprox_m = - PWA_negphiinv_underapprox_m;
            PWA_phiinv_overapprox_c = - PWA_negphiinv_underapprox_c;
        else
            PWA_phiinv_overapprox_m = zeros(T);
            PWA_phiinv_overapprox_c = zeros(T);
            lower_bound_phiinv = 0;
            upper_bound_phiinv = 0;
        end
        
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        maxlierror_logphi=5e-4;
        K = 5;
        fun_monotone_logphi = 'mono-inc';
        lower_bound_logphi = -K;
        upper_bound_logphi = norminv(exp(-maxlierror_logphi)); 
        function_handle = logphi;
        [~, ~, PWA_logphi_underapprox_m, PWA_logphi_underapprox_c] =...
            getPWAOverAndUnderApprox(lower_bound_logphi,...
            upper_bound_logphi,maxlierror_logphi,function_handle,...
            fun_monotone_logphi);

        % Compute underapproximation for log(1-x)
        log1minusx = @(z) log(1-z);
        maxlierror_log1minusx=5e-4;
        fun_monotone_log1minusx = 'mono-dec';
        lower_bound_log1minusx = log(1-Delta);
        upper_bound_log1minusx = log(normcdf(K)); 
        function_handle = log1minusx;
        [PWA_log1minusx_overapprox_m, PWA_log1minusx_overapprox_c,~,~] =...
            getPWAOverAndUnderApprox(lower_bound_log1minusx,...
            upper_bound_log1minusx,maxlierror_log1minusx,...
            function_handle,fun_monotone_log1minusx);
            
        %% Run the following scripts (which should be in the same directory) 
        %% with parameters above: 

        [ono_time_to_solve,ono_total_time,ono_opt_input_vector,...
             ono_opt_mean_X,ono_opt_val] = Ono08_IRA...
             (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,...
             mean_X_sans_input,cov_X_sans_input,state_offset,Q,R);


        [onopwl_time_to_solve,onopwl_total_time,onopwl_opt_input_vector,...
            onopwl_opt_mean_X,onopwl_opt_val] = PiecewiseAffineWithDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,cov_X_sans_input,...
            PWA_phiinv_overapprox_m,PWA_phiinv_overapprox_c,lower_bound_phiinv,upper_bound_phiinv,state_offset,Q,R);

        [pwa_time_to_solve,pwa_total_time,pwa_opt_input_vector,...
            pwa_opt_mean_X,pwa_opt_val] = PiecewiseAffineNoDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
            cov_X_sans_input,PWA_logphi_underapprox_m, PWA_logphi_underapprox_c,...
            maxlierror_logphi,lower_bound_logphi,PWA_log1minusx_overapprox_m,...
            PWA_log1minusx_overapprox_c,maxlierror_log1minusx,lower_bound_log1minusx,state_offset,Q,R);
        
        stoc_x0_flag = 0; 
        
        for run_indx = 1:3
            fprintf('Run : %d\n',run_indx);
            [blackmore_time_to_solve(run_indx),blackmore_total_time(run_indx),blackmore_opt_input_vector(:,run_indx),...
                blackmore_opt_mean_X(:,run_indx),blackmore_opt_val(run_indx)] = BlackmoreTRo11PC...
                (N,T,Delta,mean_X_sans_input,stoc_x0_flag,xtarget,...
                ulim,hbig,gbig,Ad,Bd,mean_w,cov_X_sans_input,state_offset,Q,R);
        end


    %% Plotting trajectories of each method: 

        plot_markersize = 10;
        plot_fontSize = 9;
        fig3 = figure(1);
        clf
        hold on
        h1 = plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        h2 = plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        h3 = scatter(x0(1),x0(2),5*plot_markersize,'bo','filled');
        h4 = scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'md');
        h5 = scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'bx');
        h6 = scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'r*');
        h7 = scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    10*plot_markersize, 'ks');


        
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        
        xlabel('\textbf{Position ($z_1$)}')
        ylabel('\textbf{Position ($z_2$)}')
%         title('\textbf{Trajectory}')
        set(gca,'FontSize',plot_fontSize)
        legend([h1 h2 h3 h4 h5 h6 h7],{'Safe Set',...
            'Target Set','Initial Position',...
            'Piecewise affine approach - QP',...
            'Ono2008 IRA Method ',...
            'Piecewise affine approach - MIQP',...
            sprintf('Blackmore11 PC Method, %i Particles',...
            N)},...
            'Location','SouthOutside','FontSize',plot_fontSize);
        
        box on;
        axis tight
        axis([-1.15 0.2 -1.15 0.05])
        set(gca,'XTick',sort([-1.15,get(gca,'XTick')]));
        set(gca,'YTick',sort([-1.15,get(gca,'YTick')]));
        set(gca,'FontSize',plot_fontSize)
        set(fig3,'Units','centimeters');
        set(fig3,'Position',[0 0 10 10]);
        fig3 = tightfig(fig3);
        
        handaxes3 = axes('position', [0.18 0.72 0.225 0.225]);
        hold on
        grid on
        plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        scatter(x0(1),x0(2),50,'bo','filled');
        scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    5, 'bx');
        scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    5, 'ks');
        scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    5, 'md');
        scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    5, 'r*');
                
        
        hgexport(fig3,'Figure3b',hgexport('factorystyle'),'Format', 'png')
        hgexport(fig3,'Figure3b',hgexport('factorystyle'),'Format', 'eps')
        saveas(gcf,'Figures/Fgiure3b.fig','fig');
    

    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5;  
        Q = diag(Q);
        Q = repmat(Q,1,n_mcarlo_sims);
        xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
        collection_of_input_vectors = [onopwl_opt_input_vector,...
                                       ono_opt_input_vector,...
                                       pwa_opt_input_vector,...
                                       blackmore_opt_input_vector];
        collection_of_costs = [onopwl_opt_val,...
                               ono_opt_val,...
                               pwa_opt_val,...
                               blackmore_opt_val];
        % Have a collection of string
        collection_of_method_strings = {'PWA OnoCDC2008  ',...
                                        'OnoCDC2008      ',...
                                        'PWA MILP-based  ',...
                                        'BlackmoreTRO11_A',...
                                        'BlackmoreTRO11_B',...
                                        'BlackmoreTRO11_C'};
        collection_of_time_to_solve = [onopwl_time_to_solve,...
                                       ono_time_to_solve,...
                                       pwa_time_to_solve,...
                                       blackmore_time_to_solve];
    % SReachTools for Monte-Carlo simulation
        max_rel_abserror = 0.1;
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',1-Delta,max_rel_abserror);
        for input_vec_indx = 1:6
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
            particlewise_result = all(hbig*X_mcarlo_sans_init_state <= gbig);
            prob_estim(input_vec_indx) = sum(particlewise_result)/n_mcarlo_sims;
            cost_estim(input_vec_indx) = mean(1/n_mcarlo_sims*sum(sum((X_mcarlo_sans_init_state(1:state_offset:end,:)-xtarget_mcarlo(1:state_offset:end,:)).^2.*Q))+U_vector'*R*U_vector);
            relative_abserror_in_cost(input_vec_indx) = abs(cost_estim(input_vec_indx) - true_cost)/true_cost;
            if prob_estim(input_vec_indx) >= 1-Delta && relative_abserror_in_cost(input_vec_indx) <= max_rel_abserror
                fprintf('PASSD: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                        method_str,... 
                        n_mcarlo_sims,...
                        prob_estim(input_vec_indx),...
                        relative_abserror_in_cost(input_vec_indx));
            else
                fprintf('ERROR: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                        method_str,... 
                        n_mcarlo_sims,...
                        prob_estim(input_vec_indx),...
                        relative_abserror_in_cost(input_vec_indx));
            end
        end
    % Table information
    disp('QP | IRA | MIQP | Blackmore');
    disp('Solve time');
    disp([collection_of_time_to_solve(1:3),mean(collection_of_time_to_solve(4:6))]);
    disp('Relative absolute error in cost (10^-3)');
    disp([10^3*relative_abserror_in_cost(1:3), 10^3*mean(relative_abserror_in_cost(4:6))]);        
    disp('Safety probability');
    disp([prob_estim(1:3),mean(prob_estim(4:6))]);
    save('Figure3b-0.6.mat')