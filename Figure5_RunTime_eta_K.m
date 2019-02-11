%% abyvinod-vigsiv-iterfree-ACC2019: run-time vs. eta
    % This code runs comparisons of OnoCDC08, PWLRealizationOno, and
    % BlackmoreTRo11 for a single tracjetory of a double integrator with a 
    % window closing in time and compares how the methods fair.
    % computationally. 
    % REQUIRED DEPENDENCIES: - SReachTools 
    %                        - MATLAB symbolic toolbox
    %                        - MATLAB global optimization toolbox

    %% Housekeeping: 
    
        clear;
%         close all;
        clc;
        cvx_clear;
        cd('SReachTools-private');
        srtinit
        cd('../');
        
    %% Load parameters: 

        % Time Horizon: 

            T = 10; 

        % Probability of being outside the safe set: 

            Delta = 0.2;

        % Initial conditions:   

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
            
        % Generate bounds for Blackmore: 

            
        % Disturbance parameters: 

            cov_mat_diag = 0.0001*eye(2); 
            mean_w = [0;0];
            cov_mat = kron(eye(T+1),cov_mat_diag); 
            mean_x = [0.4;0];
            x0 = RandomVector('Gaussian',mean_x,cov_mat_diag);


        % Maximum/minimum bound on input: 

            ulim = 1; 

        % Sampling time of the discrete system:

            delT = 0.25;

        % Number of particles for BlackmorePCApproach: 

            N = 100;
            
        % Generate nominal x using SReachTools:
        
            sys=getChainOfIntegLtiSystem(2, delT,...
                Polyhedron('lb',-ulim,'lb',ulim),...
                RandomVector('Gaussian',mean_w,cov_mat_diag));    
            [Ad, Bd, Gd] = getConcatMats(sys, T);
%             [HH, mean_X_sans_input, cov_X_sans_input] =...
%                 getHmatMeanCovForXSansInput(sys,x0, T);
    mean_concat_disturb = kron(ones(T,1), ...
                               sys.dist.parameters.mean);
    cov_concat_disturb  = kron(eye(T), ...
                               sys.dist.parameters.covariance);
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
            
            
    %% Fitting phiinv, logphi, and phi: 
        % Specify max piecewise interpolation errors (eta) and K values: 
            
        maxlierror_phiinv = [1E-2 1E-3 1E-4 1E-5];
        lower_bound_phiinv = [1E-5 1E-6 1E-7 1E-8];
        K = [5, 10, 15, 20, 25]; 
        maxlierror_logphi = [1E-3 1E-4 1E-5];
        maxlierror_log1minusx = maxlierror_logphi; 
        
    % Generate phiinv fits for various max errors:     
    for i = 1:length(maxlierror_phiinv)
        if Delta <= 0.5
            % Generate the PW realization of the distribution for PWLRA: 
            phiinv = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
            fun_monotone_phiinv = 'mono-inc';
            upper_bound_phiinv = Delta; 
            function_handle = @(z) -phiinv(z);
            [~,~,PWA_negphiinv_underapprox_m{i},...
                PWA_negphiinv_underapprox_c{i}] =...
                getPWAOverAndUnderApprox(lower_bound_phiinv(i),...
                upper_bound_phiinv,maxlierror_phiinv(i),...
                function_handle,fun_monotone_phiinv);
            PWA_phiinv_overapprox_m{i} = - PWA_negphiinv_underapprox_m{i};
            PWA_phiinv_overapprox_c{i} = - PWA_negphiinv_underapprox_c{i};
        else
            PWA_phiinv_overapprox_m{i} = zeros(T);
            PWA_phiinv_overapprox_c{i} = zeros(T);
            lower_bound_phiinv = 0;
            upper_bound_phiinv = 0;
        end
    end
    
    % Generate logphi and log(1-x) for various max error values: 
    
    for i = 1:length(maxlierror_logphi)
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        K_e = 5;
        fun_monotone_logphi = 'mono-inc';
        lower_bound_logphi = -K_e;
        upper_bound_logphi(i) = norminv(exp(-maxlierror_logphi(i))); 
        function_handle = logphi;
        [~, ~, PWA_logphi_underapprox_m{i}, PWA_logphi_underapprox_c{i}] =...
            getPWAOverAndUnderApprox(lower_bound_logphi,...
            upper_bound_logphi(i),maxlierror_logphi(i),function_handle,...
            fun_monotone_logphi);

        % Compute underapproximation for log(1-x)
        log1minusx = @(z) log(1-z);
        fun_monotone_log1minusx = 'mono-dec';
        lower_bound_log1minusx = log(1-Delta);
        upper_bound_log1minusx = log(normcdf(K_e)); 
        function_handle = log1minusx;
        [PWA_log1minusx_overapprox_m{i}, PWA_log1minusx_overapprox_c{i},~,~] =...
            getPWAOverAndUnderApprox(lower_bound_log1minusx,...
            upper_bound_log1minusx,maxlierror_log1minusx(i),...
            function_handle,fun_monotone_log1minusx);
    end
    

    % Generate logphi and log(1-x) for various K values: 
    
    for i = 1:length(K)
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        fun_monotone_logphiK = 'mono-inc';
        lower_bound_logphiK(i) = -K(i);
        upper_bound_logphiK = norminv(exp(-maxlierror_logphi(1))); 
        function_handle = logphi;
        [~, ~, PWA_logphi_underapprox_mK{i}, PWA_logphi_underapprox_cK{i}] =...
            getPWAOverAndUnderApprox(lower_bound_logphiK(i),...
            upper_bound_logphiK,maxlierror_logphi(1),function_handle,...
            fun_monotone_logphiK);

        % Compute underapproximation for log(1-x)
        log1minusxK = @(z) log(1-z);
        fun_monotone_log1minusxK = 'mono-dec';
        lower_bound_log1minusxK = log(1-Delta);
        upper_bound_log1minusxK(i) = log(normcdf(K(i))); 
        function_handle = log1minusxK;
        [PWA_log1minusx_overapprox_mK{i}, PWA_log1minusx_overapprox_cK{i},~,~] =...
            getPWAOverAndUnderApprox(lower_bound_log1minusxK,...
            upper_bound_log1minusxK(i),maxlierror_log1minusx(1),...
            function_handle,fun_monotone_log1minusxK);
    end


        %% Run the following scripts (which should be in the same directory) 
        %% with parameters above:
disp('Running Ono-PWL with varying max errors\n')
    for i = 1:length(maxlierror_phiinv)
            [onopwl_time_to_solve{i},onopwl_total_time{i},onopwl_opt_input_vector{i},...
                onopwl_opt_mean_X{i},onopwl_opt_val{i}] = PiecewiseAffineWithDeltaAssum...
                (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,cov_X_sans_input,...
                PWA_phiinv_overapprox_m{i},PWA_phiinv_overapprox_c{i},lower_bound_phiinv(i),upper_bound_phiinv,state_offset,Q,R);
    end
disp('Running Ono-PWL with varying max errors\n')
    for i = 1:length(maxlierror_logphi)
        [pwa_time_to_solve{i},pwa_total_time,pwa_opt_input_vector{i},...
            pwa_opt_mean_X{i},pwa_opt_val{i}] = PiecewiseAffineNoDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
            cov_X_sans_input,PWA_logphi_underapprox_m{i}, PWA_logphi_underapprox_c{i},...
            maxlierror_logphi(i),lower_bound_logphi,PWA_log1minusx_overapprox_m{i},...
            PWA_log1minusx_overapprox_c{i},maxlierror_log1minusx(i),lower_bound_log1minusx,state_offset,Q,R);
    end
    disp('Running Ono-PWL with varying max errors\n')
    for i = 1:length(K)
        [pwa_time_to_solveK{i},pwa_total_timeK{i},pwa_opt_input_vectorK{i},...
            pwa_opt_mean_XK{i},pwa_opt_valK{i}] = PiecewiseAffineNoDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
            cov_X_sans_input,PWA_logphi_underapprox_mK{i}, PWA_logphi_underapprox_cK{i},...
            maxlierror_logphi(1),lower_bound_logphiK(i),PWA_log1minusx_overapprox_mK{i},...
            PWA_log1minusx_overapprox_cK{i},maxlierror_log1minusx(1),lower_bound_log1minusx,state_offset,Q,R);
    end
    
    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5;  
        Q = diag(Q);
        Q = repmat(Q,1,n_mcarlo_sims);
        xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
        collection_of_input_vectors = [onopwl_opt_input_vector{1},...
                                       pwa_opt_input_vector{1}];
        collection_of_costs = [onopwl_opt_val{1},...
                               pwa_opt_val{1}];
        % Have a collection of string
        collection_of_method_strings = {'PWA OnoCDC2008  ',...
                                        'PWA MILP-based  '};
        collection_of_time_to_solve = [onopwl_time_to_solve{1},...
                                       pwa_time_to_solve{1}];
    % SReachTools for Monte-Carlo simulation
        max_rel_abserror = 0.1;
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',1-Delta,max_rel_abserror);
        for input_vec_indx = 1:2
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
            prob_estim(input_vec_indx) = sum(particlewise_result)/n_mcarlo_sims;
            cost_estim(input_vec_indx) = mean(1/n_mcarlo_sims*sum(sum((X_mcarlo_sans_init_state(3:end,:)-xtarget_mcarlo).^2.*Q))+U_vector'*R*U_vector);
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
        
%% Plotting the time to solve the given trajectory and marking when
%% each method fails and at what time: 
%         onopwl_time_to_solve = cell2mat(onopwl_time_to_solve);
%         pwa_time_to_solve = cell2mat(pwa_time_to_solve); 
        

        plot_markersize = 12;
        plot_fontSize = 12;
        fig2 = figure(2);
        clf
        hold on
        h1 = plot(maxlierror_phiinv,onopwl_time_to_solve,'md',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h2 = plot(maxlierror_logphi,pwa_time_to_solve,'b^','MarkerSize',...
            plot_markersize,'LineWidth',1);
%         h3 = plot(T_array,PWAWOD,'r*',...
%             'LineWidth',1,'MarkerSize',plot_markersize);
        
        xlabel('Max Error, $\eta$')
        ylabel('Time to Solve (s)')
        axis([1E-6 1E-2 0 10^0])
%         yticks([10^-3 10^-2 10^-1 10^0])
%         title('Time Horizon vs. Solve time')
        legend([h1 h2],{'Piecewise affine approach - QP',...
            'Piecewise linear approach - MIQP'},...
            'Location','EastOutside','FontSize',plot_fontSize);
        box on;
        set(gca,'FontSize',plot_fontSize)
%         set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
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

        