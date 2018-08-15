function [onopwl_time_to_solve,onopwl_total_time,onopwl_opt_input_vector,...
    onopwl_opt_mean_X,onopwl_opt_val] = PiecewiseAffineWithDeltaAssum...
    (Delta,x0,xtarget,ulim,hbig,gbig,Ad,Bd,mean_X_sans_input,cov_X_sans_input,...
    PWA_phiinv_overapprox_m,PWA_phiinv_overapprox_c,lower_bound_phiinv,upper_bound_phiinv)
    %% PiecewiseAffine Code code
    % Coder: Abraham Vinod and Vignesh Sivaramakrishnan

    disp('---------Piecewise-Affine-With-Delta-----------')
    disp(' ')

        %% House keeping
        onopwl_opt_val = nan;
        onopwl_opt_mean_X = nan(length(mean_X_sans_input), 1);
        onopwl_opt_input_vector = nan(size(Bd,2),1);


    if Delta>0.5
        warning('Skipping piecewise linear approximation (Ono''s formulation since Delta is not <0.5');
    else   



        %% Compute M --- the number of polytopic halfspaces to worry about
        onopwl_n_lin_const = size(hbig,1);

        %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
        sigma_vector = sqrt(diag(hbig*cov_X_sans_input*hbig'));


        %% Solve the feasibility problem
        tstart = tic;
        cvx_begin quiet
            variable onopwl_U_vector(size(Bd,2),1);
            variable onopwl_mean_X(length(mean_X_sans_input), 1);
            variable onopwl_deltai(onopwl_n_lin_const, 1);
            variable onopwl_norminvover(onopwl_n_lin_const, 1);

            minimize (trace(cov_X_sans_input(1:2:end,1:2:end)) + (xtarget-onopwl_mean_X(1:2:end))'*(xtarget-onopwl_mean_X(1:2:end)))
            subject to
                onopwl_mean_X == Ad*x0+  Bd * onopwl_U_vector;
                abs(onopwl_U_vector) <= ulim;
                for onopwl_deltai_indx=1:onopwl_n_lin_const
                    % Implements piecewise maximization
                    onopwl_norminvover(onopwl_deltai_indx) >=...
                       PWA_phiinv_overapprox_m.* onopwl_deltai(onopwl_deltai_indx) + PWA_phiinv_overapprox_c;
                end
                hbig*onopwl_mean_X + sigma_vector.* onopwl_norminvover <= gbig;
                onopwl_deltai >= lower_bound_phiinv; 
                onopwl_deltai <= upper_bound_phiinv; 
                sum(onopwl_deltai) <= Delta;
         t1 = toc(tstart);
         cvx_end;
         t2 = toc(tstart);
         onopwl_time_to_solve = t2 - t1;
         onopwl_total_time = cvx_cputime;

        %% Overwrite the solutions
        if strcmpi(cvx_status, 'Solved')
            onopwl_opt_input_vector = onopwl_U_vector;
            onopwl_opt_mean_X = onopwl_mean_X;
            onopwl_opt_val = cvx_optval;
        else
            error('PWA-based Ono failed');
        end

    disp('------------------------------------')
    fprintf('Total CVX Solve Time: %1.4f seconds\n\n',onopwl_time_to_solve)
    disp('------------------------------------')
    fprintf('Total Run Time: %1.4f seconds\n',onopwl_total_time)
    end
end