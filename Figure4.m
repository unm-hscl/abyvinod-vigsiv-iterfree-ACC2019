%% abyvinod-vigsiv-iterfree-ACC2019: Figure 2
    % This code runs comparisons of OnoCDC08, PWLRealizationOno, and
    % BlackmoreTRo11 for a single tracjetory of a double integrator with a 
    % window closing in time and compares how the methods fair with
    % increasing time horizon.
    % computationally. 
    % REQUIRED DEPENDENCIES: - SReachTools 
    %                        - MATLAB symbolic toolbox
    %                        - MATLAB global optimization toolbox

    %% Housekeeping: 
    
        clear;
        close all;
        clc;
        cvx_clear;
        
    %% Load parameters: 

        % Time Horizons: 

            T_array = 10:10:60; 
            
                % Disturbance parameters: 

            cov_mat_diag = 0.0001*diag([1 0;]); 
            mean_w = [0;0];


        % Maximum/minimum bound on input: 

            ulim = 1; 

        % Sampling time of the discrete system:

            delT = 0.25;

        % Number of particles for BlackmorePCApproach: 

            N = 50;

        % Probability of being outside the safe set: 

            Delta = 0.2;

        % Initial conditions:   

            x0 = [0.4;0];
            
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
        
    for i=1:length(T_array)
        % Specify the time horizon: 
        
            T = T_array(i);
        
        % Initial conditions that are time varying:  
            xtarget = linspace(-0.4,0.2,T)';
            xtarget = reshape([xtarget'; zeros(size(xtarget'))],[],1);
            Q = kron(eye(T),diag([10 1])); 
            R = kron(eye(T),0.001);
            
        % Bounds on the safe set: 
        
            h = [-1 0; 1 0;];
            gb = linspace(0.5,0.2, T);
            
        % Generate bounds for Ono and PWA: 
            hbig = kron(eye(T),h);
            gbig = kron(gb,[1,1])';
            state_offset = 1;
  

     %% Prepare system matrices: 

        % Generate a large cov_mat for the optimizaiton problem:
        
            cov_mat = kron(eye(T+1),cov_mat_diag); 

        % Generate nominal x using SReachTools:
        
            sys=getChainOfIntegLtiSystem(2, delT,...
                Polyhedron('lb',-ulim,'lb',ulim),...
                RandomVector('Gaussian',mean_w,cov_mat_diag));    
            [Ad, Bd, Gd] = getConcatMats(sys, T);
    if isa(x0,'RandomVector')
        % TODO: Waiting for a reference --- but essentially propagation of mean
        % and covariance
        mean_X_sans_input = Ad * x0.parameters.mean + Gd *...
            mean_concat_disturb;
        cov_X_sans_input = Ad * x0.parameters.covariance * Ad' + ...
            Gd * cov_concat_disturb * Gd';
    else
        % Computation of mean and covariance of X (sans input) by (17),LCSS 2017
        mean_X_sans_input = Ad * initial_state + Ad * mean_concat_disturb;
        cov_X_sans_input = Gd * cov_concat_disturb * Gd';
    end  
            




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

        [blackmore_time_to_solve,blackmore_total_time,blackmore_opt_input_vector,...
            blackmore_opt_mean_X,blackmore_opt_val] = BlackmoreTRo11PC...
            (N,T,Delta,mean_X_sans_input,stoc_x0_flag,xtarget,...
            ulim,hbig,gbig,Ad,Bd,mean_w,cov_X_sans_input,state_offset,Q,R);

    

    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5; 
        Q = diag(Q);
        Q = repmat(Q,1,n_mcarlo_sims);
        xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
        collection_of_input_vectors = [blackmore_opt_input_vector, ono_opt_input_vector, onopwl_opt_input_vector, pwa_opt_input_vector];
        collection_of_costs = [blackmore_opt_val, ono_opt_val, onopwl_opt_val,pwa_opt_val];
        % Have a collection of string
        collection_of_method_strings = {'BlackmoreTRO11',...
                                        'OnoCDC2008    ',...
                                        'PWA OnoCDC2008',...
                                        'PWA MILP-based'};
    % SReachTools for Monte-Carlo simulation
        max_rel_abserror = 0.1;
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',1-Delta,max_rel_abserror);
        for input_vec_indx = 1:4
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
            particlewise_result = all(hbig*X_mcarlo_sans_init_state(3:end,:) <= gbig);
            prob_estim = sum(particlewise_result)/n_mcarlo_sims;
            cost_estim(input_vec_indx) = mean(1/n_mcarlo_sims*sum(sum((X_mcarlo_sans_init_state(3:end,:)-xtarget_mcarlo).^2.*Q))+U_vector'*R*U_vector); 
            relative_abserror_in_cost = abs(cost_estim(input_vec_indx) - true_cost)/true_cost;
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
        
        % log time to solve for each method: 
        
            OnoTTS(i) = ono_time_to_solve;
            BlackmoreTTS(i) = blackmore_time_to_solve;
            PWAWD(i) = onopwl_time_to_solve;
            PWAWOD(i) = pwa_time_to_solve;
    end
    
%% Plotting the time to solve the given trajectory and marking when
%% each method fails and at what time: 
        plot_markersize = 9;
        plot_fontSize = 9;
        fig2 = figure(2);
        clf
        hold on
        h1 = plot(T_array,PWAWD,'md',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h2 = plot(T_array,OnoTTS,'b^','MarkerSize',...
            plot_markersize,'LineWidth',1);
        h3 = plot(T_array,PWAWOD,'r*',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h4 = plot(T_array,BlackmoreTTS,'ks',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        
        xlabel('Time Horizon')
        ylabel('Time to Solve (seconds)')
        axis([5 60 0 3*10^3])
        yticks([10^0 10^1 10^2 10^3])
%         title('Time Horizon vs. Solve time')
        legend([h1 h2 h3 h4],{'Piecewise affine approach - QP',...
            'Iterative risk allocation (IRA)',...
            'Piecewise affine approach - MIQP',...
            sprintf('Particle control (PC), %i Particles',N)},...
            'Location','EastOutside','FontSize',plot_fontSize);
        box on;
        set(gca,'FontSize',plot_fontSize)
        set(gca, 'YScale', 'log')
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        set(fig2,'Units','centimeters');
        set(fig2,'Position',[0 0 20 4.8]);
        grid on
        fig1 = tightfig(fig2);
        hgexport(fig2,'Figure2',hgexport('factorystyle'),'Format', 'png')
        hgexport(fig1,'Figure2',hgexport('factorystyle'),'Format', 'eps')
        saveas(gcf,'Figures/Fgiure2.fig','fig');
