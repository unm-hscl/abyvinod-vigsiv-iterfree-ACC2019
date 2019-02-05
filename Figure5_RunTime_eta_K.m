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
            
    %% Fitting phiinv, logphi, and phi: 
        % Specify max piecewise interpolation errors (eta) and K values: 
            
        maxlierror_phiinv = [1E-2 1E-3 1E-4 1E-5];
        K = [5, 10, 15, 20, 25]; 
        maxlierror_logphi = [1E-2 1E-3 1E-4 1E-5];
        maxlierror_log1minusx = maxlierror_logphi; 
        
    % Generate phiinv fits for various max errors:     
    for i = 1:length(maxlierror_phiinv)
        if Delta <= 0.5
            % Generate the PW realization of the distribution for PWLRA: 
            phiinv = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
            fun_monotone_phiinv = 'mono-inc';
            lower_bound_phiinv = 1E-5;
            upper_bound_phiinv = Delta; 
            function_handle = @(z) -phiinv(z);
            [~,~,PWA_negphiinv_underapprox_m{i},...
                PWA_negphiinv_underapprox_c{i}] =...
                getPWAOverAndUnderApprox(lower_bound_phiinv,...
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
    
    for i = 1:length(maxlierror_phiinv)
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        K = 5;
        fun_monotone_logphi = 'mono-inc';
        lower_bound_logphi = -K;
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
        upper_bound_log1minusx = log(normcdf(K)); 
        function_handle = log1minusx;
        [PWA_log1minusx_overapprox_m{i}, PWA_log1minusx_overapprox_c{i},~,~] =...
            getPWAOverAndUnderApprox(lower_bound_log1minusx,...
            upper_bound_log1minusx,maxlierror_log1minusx(i),...
            function_handle,fun_monotone_log1minusx);
    end
    

    % Generate logphi and log(1-x) for various K values: 
    
    for i = 1:length(maxlierror_phiinv)
        % Compute underapproximation for log(Phi(x))
        logphi = @(z) log(normcdf(z));
        K = 5;
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

    for i = 1:length(maxlierror_phiinv)
            [onopwl_time_to_solve(:,i),onopwl_total_time(:,i),onopwl_opt_input_vector(:,i),...
                onopwl_opt_mean_X(:,i),onopwl_opt_val(:,i)] = PiecewiseAffineWithDeltaAssum...
                (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,cov_X_sans_input,...
                PWA_phiinv_overapprox_m(:,i),PWA_phiinv_overapprox_c(:,i),lower_bound_phiinv(:,i),upper_bound_phiinv(:,i),state_offset,Q,R);
    end

    for i = 1:length(maxlierror_phiinv)
        [pwa_time_to_solve(:,i),pwa_total_time,pwa_opt_input_vector(:,i),...
            pwa_opt_mean_X(:,i),pwa_opt_val(:,i)] = PiecewiseAffineNoDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
            cov_X_sans_input,PWA_logphi_underapprox_m(:,i), PWA_logphi_underapprox_c(:,i),...
            maxlierror_logphi(i),lower_bound_logphi,PWA_log1minusx_overapprox_m(:,i),...
            PWA_log1minusx_overapprox_c(:,i),maxlierror_log1minusx(i),lower_bound_log1minusx,state_offset,Q,R);
        
        [pwa_time_to_solve,pwa_total_time,pwa_opt_input_vector,...
            pwa_opt_mean_X,pwa_opt_val] = PiecewiseAffineNoDeltaAssum...
            (Delta,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
            cov_X_sans_input,PWA_logphi_underapprox_m(:,i), PWA_logphi_underapprox_c(:,i),...
            maxlierror_logphi(i),lower_bound_logphi,PWA_log1minusx_overapprox_m(:,i),...
            PWA_log1minusx_overapprox_c(:,i),maxlierror_log1minusx(i),lower_bound_log1minusx,state_offset,Q,R);
    end
        