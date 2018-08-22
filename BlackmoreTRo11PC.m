function [blackmore_time_to_solve,blackmore_total_time,blackmore_opt_input_vector,...
    blackmore_opt_mean_X,blackmore_opt_val] = BlackmoreTRo11PC...
    (N,T,Delta,x0,xtarget,ulim,hbig,gbig,Ad,Bd,mean_w,cov_X_sans_input,state_offset)
    %% Blackmore TRo 2011 Code to stay in a feasible set. 
    % Coder: Vignesh Sivaramakrishnan
    
    % The following housekeeping if statement does not run Blackmore for
    % longer than 40 seconds for it takes way to long to compute. 
    


    % System matrices: 

        disp('---------BlackmorePC11-----------')
        disp(' ')
        fprintf('No. of particles: %d\n',N);

        large_constant = 5000;

    % Vectorize the target trajectory for the optimization problem. 
        xtargetbig = repmat(xtarget,1,N);

    % Randomly generate the disturbance vector from the standard normal.

        mean_GdTimesW = repmat(mean_w,T,N);

        GdTimesW = mvnrnd(mean_GdTimesW', cov_X_sans_input)';
        
    if T >= 40
        
        blackmore_opt_mean_X = nan(length(mean_GdTimesW),1);
        blackmore_opt_val = nan;
        blackmore_opt_input_vector = nan(size(Bd,2),1); 
        blackmore_time_to_solve = nan;
        blackmore_total_time = nan;
        warning('NOTE: BlackmoreTRo11 takes a very long time for long time horizons!!')
        return;
    else


    %% Run optimization problem for an optimal control policy
    % We run an optimization problem to determine the control policy over the
    % time horizon T.
        tstart = tic;
        cvx_clear
            cvx_precision BEST
        cvx_begin quiet
            variable U_vector(size(Bd,2),1);
            variable xBl(size(mean_GdTimesW,1),N);
            variable mean_X(size(mean_GdTimesW,1),1);
            variable d(N) binary;

            minimize (sum(sum((xBl(1:state_offset:end,:)-xtargetbig(1:state_offset:end,:)).^2))/N);

            subject to
              mean_X == Ad*x0+ Bd*U_vector;

              xBl(1:end,1:N) == GdTimesW+repmat(mean_X,1,N);

              abs(U_vector) <= ulim;

              for i = 1:N
                  hbig*xBl(:,i) - gbig <= large_constant*(d(i));
              end
              1/N*sum(d)<=Delta;

        t1 = toc(tstart);
        cvx_end;
        t2 = toc(tstart);
        blackmore_time_to_solve = t2 - t1;
        blackmore_total_time = cvx_cputime;

        if strcmpi(cvx_status,'Solved')
            blackmore_opt_mean_X = mean_X;
            blackmore_opt_val = cvx_optval;
            blackmore_opt_input_vector = U_vector;            
        else
            blackmore_opt_mean_X = nan(length(mean_GdTimesW),1);
            blackmore_opt_val = nan;
            blackmore_opt_input_vector = nan(size(Bd,2),1);         
        end
    end

        fprintf('Total CVX Run Time for %1i particles: %1.4f seconds\n',...
            N,cvx_cputime)
        disp('------------------------------------')
        fprintf('Total CVX Solve Time for %1i particles: %1.4f seconds\n'...
            ,N,blackmore_time_to_solve)

        d = full(d);
end
