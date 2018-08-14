%% abyvnod-vigsiv-iterfree-ACC2019: Figure 1
    % This code runs comparisons of OnoCDC08, PWLRealizationOno, and
    % BlackmoreTRo11 for a single tracjetory of a double integrator with a 
    % window closing in time and compares how the methods fair.
    % computationally. 
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
        
    %% Load parameters: 

        % Time Horizon: 

            T = 50; 

        % Probability of being outside the safe set: 

            Delta = 0.2;

        % Initial conditions:   

            x0 = [0.4;0];
            xtarget = linspace(-0.4,0.2,T)'; 
            
        % Bounds on the safe set: 
        
            h = [-1 0; 1 0;];
            gb = linspace(0.5,0.2, T);
            
        % Disturbance parameters: 

            cov_mat_diag = 0.0001*diag([1 0;]); 
            mean_w = [0;0];


        % Maximum/minimum bound on input: 

            ulim = 1; 

        % Sampling time of the discrete system:

            delT = 0.25;

        % Number of particles for BlackmorePCApproach: 

            N = 20;
        

     %% Prepare system matrices: 

        % Generate a large cov_mat for the optimizaiton problem:
        
            cov_mat = kron(eye(T+1),cov_mat_diag); 

        % Generate nominal x using SReachTools:
        
            sys=getChainOfIntegLtiSystem(2, delT,...
                Polyhedron('lb',-ulim,'lb',ulim),...
                RandomVector('Gaussian',mean_w,cov_mat_diag));    
            [Ad, Bd, Gd] = getConcatMats(sys, T);
            [~, mean_X_sans_input, cov_X_sans_input] =...
                getHmatMeanCovForXSansInput(sys, x0, T);  
            
        if Delta <= 0.5
            % Generate the PW realization of the distribution for PWLRA: 
            maxlierror_phiinv=1e-2;
            phiinv = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
            fun_monotone_phiinv = 'mono-inc';
            lower_bound_phiinv = 1E-5;
            upper_bound_phiinv = Delta; 
            function_handle = @(z) -phiinv(z);
            [PWA_overapprox_phiinv_m, PWA_overapprox_phiinv_c, ~, ~] = getPWAOverAndUnderApprox(lower_bound_phiinv,upper_bound_phiinv,maxlierror_phiinv,function_handle,fun_monotone_phiinv);
        end
        
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        maxlierror_logphi=1e-3;
        K = 5;
        fun_monotone_logphi = 'mono-inc';
        lower_bound_logphi = -K;
        upper_bound_logphi = norminv(exp(-maxlierror_logphi)); 
        function_handle = logphi;
        [~, ~, PWA_logphi_underapprox_m, PWA_logphi_underapprox_c] = getPWAOverAndUnderApprox(lower_bound_logphi,upper_bound_logphi,maxlierror_logphi,function_handle,fun_monotone_logphi);

        % Compute underapproximation for log(1-x)
        log1minusx = @(z) log(1-z);
        maxlierror_log1minusx=1e-3;
        fun_monotone_log1minusx = 'mono-dec';
        lower_bound_log1minusx = log(1-Delta);
        upper_bound_log1minusx = log(normcdf(K)); 
        function_handle = log1minusx;
        [PWA_log1minusx_overapprox_m, PWA_log1minusx_overapprox_c,~,~] = getPWAOverAndUnderApprox(lower_bound_log1minusx,upper_bound_log1minusx,maxlierror_log1minusx,function_handle,fun_monotone_log1minusx);

        % Cost ratios b/n input and state --- scalarization term:
        input_state_ratio = 0.0001;

        %% Run the following scripts (which should be in the same directory) 
        %% with parameters above: 

        try
            Ono08_IRA
        catch
            disp('Ono''s method failed');
        end
        disp(' ');
        disp(' ');
        PiecewiseAffineWithDeltaAssum
        disp(' ');
        disp(' ');
        PiecewiseAffineNoDeltaAssum
        disp(' ');
        disp(' ');
        BlackmoreTRo11PCOno08Mod
        
    %% Plotting trajectories of each method: 

        plot_markersize = 9;
        plot_fontSize = 14;
        fig1 = figure(1);
        clf
        hold on
        polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
            [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
        P = Polyhedron('V',polyvertex);
        h1 = plot(P,'alpha',0.1);
        h11 = scatter(1,x0(1),plot_markersize*10,'filled');
        h2 = plot(2:(T+1),xtarget,'go','MarkerSize',...
            plot_markersize,'LineWidth',2);
        h3 = plot(2:(T+1),ono_opt_mean_X(1:2:end),'bx',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h4 = plot(2:(T+1),blackmore_opt_mean_X(1:2:end,1),...
            'ks','MarkerSize',plot_markersize);
        h5 = plot(2:(T+1),onopwl_opt_mean_X(1:2:end),'md',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h6 = plot(2:(T+1),pwa_opt_mean_X(1:2:end),'r*',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        
        xlabel('\textbf{Time (seconds)}')
        ylabel('\textbf{Position (x)}')
        title('\textbf{Trajectory}')
        legend([h1 h11 h2 h3 h4 h5 h6],{'Safe Set',...
            'Initial state','Target Trajectory',...
            sprintf('Ono2008 IRA Method (Cost: %1.3f)',...
            ono_opt_val),sprintf('Blackmore11 PC Method, %i Particles (Cost: %1.3f)',...
            N,blackmore_opt_val),...
            sprintf('Piecewise affine approach - cvx (Cost: %1.3f)',...
            onopwl_opt_val),...
            sprintf('Piecewise linear approach - milp (Cost: %1.3f)',...
            pwa_opt_val)},'Location','SouthOutside','FontSize',plot_fontSize);
        box on;
%         set(fig1,'PaperUnits','centimeters');
%         set(fig1,'PaperPosition',[0 0 8.8 8.8]);
        set(gca,'FontSize',plot_fontSize)
%         fig1 = tightfig(fig1);
        hgexport(fig1,'Figure1',hgexport('factorystyle'),'Format', 'png')
    

    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5;
        % FIXME: Shouldn't be redefining this again!
        hbig = kron(eye(T),h);
        gbig = kron(gb,[1,1])';    
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
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',Delta,max_rel_abserror);
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
        particlewise_result = all(hbig*X_mcarlo_sans_init_state <= gbig);
        prob_estim = sum(particlewise_result)/n_mcarlo_sims;
        cost_estim = mean(sum((X_mcarlo_sans_init_state(1:2:end,:)-xtarget_mcarlo).^2));    
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
