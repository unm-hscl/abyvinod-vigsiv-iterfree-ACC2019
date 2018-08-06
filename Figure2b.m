%% abyvnod-vigsiv-iterfree-ACC2019
% This code runs comparisons of OnoCDC08, PWLRealizationOno, and for a
% single tracjetory and compares how they fair as the time horizon is
% increased.
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

    %% Load specifications (Ones which do not depend on the Time Horzion): 

        % Time Horizon array: 

            T_array = 20:10:80;
            T_toplot = 40;

        % Probability of being outside the safe set: 

            Delta = 0.2;

        % Initial condition:   

            x0 = [0.4;0];
            xtar = [0 0];
            
        % Disturbance parameters: 

            cov_mat_diag = 0.0001*diag([1 0;]); 
            mean_w = [0;0];
            
        % Bounds on the safe set (Non-time dependent vars): 

            h = [-1 0; 1 0;];
            gdes = [0.5 0.2];
        % Maximum/minimum bound on input: 

            ulim = 1; 

        % Sampling time of the discrete system:

            delT = 0.25;
      
        % Generate the PW realization of the distribution for PWLRA: 
            maxlierror = 1E-3;
            lbdelta = 1E-4;
            [onopwl_invcdf_approx_m, onopwl_invcdf_approx_c,...
            onopwl_lb_deltai] = RolleLerpClosedForm(Delta,lbdelta,...
            maxlierror);

        % Cost ratios b/n input and state --- scalarization term:
            input_state_ratio = 0.0001;
            
for i = 1:length(T_array)
    cvx_clear
    T = T_array(i);
    
    %% Load specifications (Ones which do depend on the Time Horizon): 
            
        % Desired target:
        
            xtarget = linspace(xtar(1),xtar(2),T)'; 

        % Bounds on the safe set (Non-time dependent vars):

            g = linspace(gdes(1),gdes(2), T);


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


    %% Run the following scripts (which should be in the same directory) 
    %% with parameters above: 

        try
            Ono08_IRA
            OnoIRAFlag(i) = 1;
        catch
            disp('Ono''s method failed');
            OnoIRAFlag(i) = 0; 
        end
        disp(' ');
        disp(' ');
        try
            PiecewiseLinearRA
            PWLAFlag(i) = 1;
        catch
            disp('PWLA method failed');
            PWLAFlag(i) = 0; 
        end

    %% Monte-Carlo simulation using SReachTools
        n_mcarlo_sims = 1e5;
        % FIXME: Shouldn't be redefining this again!
        hbig = kron(eye(T),h);
        gbig = kron(g,[1,1])';    
        xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
        collection_of_input_vectors = [ono_opt_input_vector,...
            onopwl_opt_input_vector];
        collection_of_costs = [ono_opt_val, onopwl_opt_val];
        % Have a collection of strings:
        collection_of_method_strings = {'OnoCDC2008    ',...
                                        'PWL OnoCDC2008'};
    % SReachTools for Monte-Carlo simulation:
        max_rel_abserror = 0.1;
        fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',Delta,max_rel_abserror);
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
        
        
    % Log time to solve for each method and indicate if it passed or
    % failed:
    if T_toplot == T
        Ono_x= ono_opt_mean_X;
        PWLA_x= onopwl_opt_mean_X;
    end
    Ono08TTS(i) = ono_time_to_solve;
    PWLRATTS(i) = onopwl_time_to_solve;
    
    Ono08OPT(i) = ono_opt_val;
    PWLRAOPT(i) = onopwl_opt_val;
    
end

    %% Plotting trajectories of each method for a specified Time Horizon: 
        hbig = kron(eye(T_toplot),h);
        g = linspace(gdes(1),gdes(2), T_toplot);
        gbig = kron(g,[1,1])';  
        xtarget = linspace(xtar(1),xtar(2),T_toplot)';
        plot_markersize = 10;
        plot_fontSize = 10;
        figure(1)
        clf
        hold on
        polyvertex =[1:T_toplot+1,1:T_toplot+1;[-gbig(1),-gbig(1:2:end)'],...
            [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
        P = Polyhedron('V',polyvertex);
        h1 = plot(P,'alpha',0.1);
        h11 = scatter(1,x0(1),plot_markersize*10,'filled');
        h2 = plot(2:(T_toplot+1),xtarget,'go','MarkerSize',...
            plot_markersize,'LineWidth',2);
        h3 = plot(2:(T_toplot+1),Ono_x(1:2:end),'bx',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h4 = plot(2:(T_toplot+1),Ono_x(1:2:end),'md',...
            'LineWidth',1,'MarkerSize',plot_markersize);    
        xlabel('Time (seconds)')
        ylabel('X Position')
        title(sprintf('Trajectory at T = %i seconds',T_toplot))
        legend([h1 h11 h2 h3 h4],{'Safe Set',...
            'Initial state','Target Trajectory',...
            sprintf('Ono2008 IRA Method (Cost: %1.3f)',...
            Ono08OPT(find(T_array==T_toplot))),...
            sprintf('Piecewise linear approach (Cost: %1.3f)',...
            PWLRAOPT(find(T_array==T_toplot)))});
        box on;
        set(gca,'FontSize',plot_fontSize)

        figure(1);
        print -dpng Figure2b1.png

        
    %% Plotting the time to solve the given trajectory and marking when
    %% each method fails and at what time: 
        plot_markersize = 10;
        plot_fontSize = 10;
        figure(2)
        clf
        hold on
        h5 = plot(T_array,Ono08TTS,'go','MarkerSize',...
            plot_markersize,'LineWidth',2);
        h6 = plot(T_array,PWLRATTS,'b^',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        legFlag = 0;
        Altleg = 0;
        if all(OnoIRAFlag) < 1 || all(PWLAFlag) < 1
            for i = 1:length(T_array)
               if  OnoIRAFlag(i) == 0
                   h7 = plot(T_array(i),Ono08TTS(i),'rx','MarkerSize',...
                        plot_markersize,'LineWidth',2);
                    Altleg = 1;
               end

               if  PWLAFlag(i) == 0
                   h7 = plot(T_array(i),PWLRATTS(i),'rx','MarkerSize',...
                        plot_markersize,'LineWidth',2);
                    Altleg = 1;
               end
                if Altleg == 1
                   legend([h5 h6 h7],'Ono2008 IRA Method',...
                    'Piecewise linear approach','Time horizon at which the method failed','Location','best');
                end
               legFlag = 1;
            end
        end
        
        xlabel('Time Horizon (seconds)')
        ylabel('Time to Solve (seconds)')
        title('Time Horizon vs. Solve time')
        if legFlag ~= 1
            legend([h5 h6],'Ono2008 IRA Method',...
                'Piecewise linear approach','Location','best');
        end
        box on;
        set(gca,'FontSize',plot_fontSize)
        print -dpng Figure2b2.png
        
