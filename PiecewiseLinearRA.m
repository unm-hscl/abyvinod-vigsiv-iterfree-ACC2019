%% Ono_IRA 2008 code
% Coder: Abraham Vinod and Vignesh Sivaramakrishnan

disp('---------Piecewise-linear-----------')
disp(' ')

if Delta>0.5
    warning('Skipping piecewise linear approximation (Ono''s formulation since Delta is not <0.5');
else   


% Generate bounds: 
    hbig = kron(eye(T),h);
    gbig = kron(g,[1,1])';


    %% House keeping
    onopwl_opt_val = nan;
    onopwl_opt_mean_X = nan(length(mean_X_sans_input), 1);
    onopwl_opt_input_vector = nan(size(Bd,1),1);

    %% Compute M --- the number of polytopic halfspaces to worry about
    onopwl_n_lin_const = size(hbig,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = sqrt(diag(hbig*cov_X_sans_input(3:end,3:end)*hbig'));

    % TODO: Translate desired_accuracy to piecewise_count
    [onopwl_invcdf_approx_m, onopwl_invcdf_approx_c, onopwl_lb_deltai, max_error_onopwl]=...
        computeNormCdfInvOverApprox();

    % onopwl approach introduces an artifical conservativeness of max_gap *
    % n_lin_const
    onopwl_artificial_error = max(max_error_onopwl, onopwl_lb_deltai) * onopwl_n_lin_const;
    if  onopwl_artificial_error > desired_accuracy 
        warning(sprintf('Required accuracy: %1.3e | Error due to onopwl: %1.3e', desired_accuracy, onopwl_artificial_error));
    end

    %% Solve the feasibility problem
    tstart = tic;
    cvx_begin quiet
        variable onopwl_U_vector(size(Bd,2),1);
        variable onopwl_mean_X(length(mean_X_sans_input), 1);
        variable onopwl_deltai(onopwl_n_lin_const, 1);
        variable onopwl_norminvover(onopwl_n_lin_const, 1);
        minimize (input_state_ratio*sum(abs(onopwl_U_vector))/(ulim*T) + sum(abs(onopwl_mean_X(3:2:end)-xtarget))/(2*g(1)*T));
        subject to
            onopwl_mean_X == mean_X_sans_input + Bd * onopwl_U_vector;
            abs(onopwl_U_vector) <= ulim;
            for onopwl_deltai_indx=1:onopwl_n_lin_const
                onopwl_norminvover(onopwl_deltai_indx) >= onopwl_invcdf_approx_m.*...
                    onopwl_deltai(onopwl_deltai_indx) + onopwl_invcdf_approx_c; 
            end
            hbig*onopwl_mean_X(3:end) + sigma_vector.* onopwl_norminvover <= gbig;
            onopwl_deltai >= onopwl_lb_deltai;
            onopwl_deltai <= 0.5;
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
    end

disp('------------------------------------')
fprintf('Total CVX Solve Time: %1.4f seconds\n\n',onopwl_time_to_solve)
disp('------------------------------------')
fprintf('Total Run Time: %1.4f seconds\n',onopwl_total_time)
end
