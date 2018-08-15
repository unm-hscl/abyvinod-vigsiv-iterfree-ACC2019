function [pwa_time_to_solve,pwa_total_time,pwa_opt_input_vector,...
    pwa_opt_mean_X,pwa_opt_val] = PiecewiseAffineNoDeltaAssum...
    (Delta,x0,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,...
    cov_X_sans_input,PWA_logphi_underapprox_m, PWA_logphi_underapprox_c,...
    maxlierror_logphi,lower_bound_logphi,PWA_log1minusx_overapprox_m,...
    PWA_log1minusx_overapprox_c,maxlierror_log1minusx,lower_bound_log1minusx)
    %% PiecewiseAffineNoDeltaAssump
    % Coder: Abraham Vinod and Vignesh Sivaramakrishnan

    disp('---------Piecewise-Affine-No-Delta-----------')
    disp(' ')

    %% House keeping
    myeps_MILP = 1e-10;


    %% Compute M --- the number of polytopic halfspaces to worry about
    pwa_n_lin_const = size(hbig,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = sqrt(diag(hbig*cov_X_sans_input*hbig'));

    %% Compute bounds for the MILP
    Mlb = PWA_log1minusx_overapprox_m * Delta + PWA_log1minusx_overapprox_c + maxlierror_logphi;
    Mub = PWA_log1minusx_overapprox_c - log( 1 - Delta);

    %% Solve the feasibility problem
    PWA_n_log1minusx = length(PWA_log1minusx_overapprox_m);
    tstart = tic;
    cvx_begin quiet
        variable pwa_U_vector(size(Bd,2),1);
        variable pwa_mean_X(length(mean_X_sans_input), 1);
        variable pwa_deltai(pwa_n_lin_const, 1) nonnegative;
        variable pwa_slackvar(pwa_n_lin_const, 1); 
        variable pwa_bin_log1minusx(PWA_n_log1minusx,pwa_n_lin_const) binary;
        minimize (trace(cov_X_sans_input) + (xtarget-pwa_mean_X)'*(xtarget-pwa_mean_X))
        subject to
            % (17b)
            sum(pwa_deltai) <= Delta;
            pwa_mean_X == Ad*x0+  Bd * pwa_U_vector;
            abs(pwa_U_vector) <= ulim;

            % (17c)
            pwa_deltai <= Delta;
            pwa_slackvar >= log(1-Delta);
            pwa_slackvar <= - maxlierror_logphi;

            % (17d), (17e)
            for pwa_deltai_indx=1:pwa_n_lin_const
                PWA_logphi_underapprox_m * ((gbig(pwa_deltai_indx) - hbig(pwa_deltai_indx,:) * pwa_mean_X)/sigma_vector(pwa_deltai_indx)) +...
                    PWA_logphi_underapprox_c >= pwa_slackvar(pwa_deltai_indx); 
            end
            (gbig - hbig * pwa_mean_X)>= lower_bound_logphi.*sigma_vector;

            % (17f), (17g)
            for pwa_deltai_indx=1:pwa_n_lin_const
                PWA_log1minusx_overapprox_m * pwa_deltai(pwa_deltai_indx) +...
                    PWA_log1minusx_overapprox_c - repmat(pwa_slackvar(pwa_deltai_indx),1,PWA_n_log1minusx) <= Mub.*(ones(1,PWA_n_log1minusx) - pwa_bin_log1minusx(:,pwa_deltai_indx)'); 
                PWA_log1minusx_overapprox_m * pwa_deltai(pwa_deltai_indx) +...
                    PWA_log1minusx_overapprox_c - repmat(pwa_slackvar(pwa_deltai_indx),1,PWA_n_log1minusx) >= repmat(myeps_MILP,1,PWA_n_log1minusx) + (Mlb - myeps_MILP).* pwa_bin_log1minusx(:,pwa_deltai_indx)'; 
            end

            % (17i) --- sum works columnwise
            sum(pwa_bin_log1minusx,1) >= 1;
     t1 = toc(tstart);
     cvx_end;
     t2 = toc(tstart);


    %% Overwrite the solutions
    if strcmpi(cvx_status, 'Solved')
        pwa_opt_input_vector = pwa_U_vector;
        pwa_opt_mean_X = pwa_mean_X;
        pwa_opt_val = cvx_optval;
        pwa_time_to_solve = t2 - t1;
        pwa_total_time = cvx_cputime;
    else
        pwa_opt_val = nan;
        pwa_opt_mean_X = nan(length(mean_X_sans_input), 1);
        pwa_opt_input_vector = nan(size(Bd,2),1);
        pwa_time_to_solve = 0;
        pwa_total_time = 0;
        warning('Piecewise MILP failed');
    end
    

    disp('------------------------------------')
    fprintf('Total CVX Solve Time: %1.4f seconds\n\n',pwa_time_to_solve)
    disp('------------------------------------')
    fprintf('Total Run Time: %1.4f seconds\n',pwa_total_time)
end
