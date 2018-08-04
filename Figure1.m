%% abyvnod-vigsiv-iterfree-ACC2019
    % This code runs comparisons of the 


    clear;
    close all;
    clc;
    cvx_clear;
    
    cd('SReachTools-private');
    srtinit
    cd('../');
    % Load parameters: 

    % Time Horizon: 

        T = 30; 

    % Probability of being outside the safe set: 

        Delta = 0.2;
%       Delta = 0.6;

    % Initial conditions:   

%         x0 = [0.2;0];
%         xtarget = linspace(0.2,0.2,T)';   
        x0 = [0.4;0];
        xtarget = linspace(-0.4,0.2,T)'; 
    %% Bounds on the safe set
    h = [-1 0; 1 0;];
%     g = linspace(0.5,0.5, T);
    g = linspace(0.5,0.2, T);
    % Disturbance parameters: 

        cov_mat_diag = 0.0001*diag([1 0;]); 
        mean_w = [0;0];


    % Maximum/minimum bound on input: 

        ulim = 1; 

    % Sampling time of the discrete system:

        delT = 0.25;

    % Number of particles for BlackmorePCApproach: 

        N = 50;
    % Desired accuracy
    desired_accuracy = 0.001; 

    % Prep
    % System matrices: 

% Generate a large cov_mat for the optimizaiton problem.
    cov_mat = kron(eye(T+1),cov_mat_diag); 

% Generate nominal x (Note this is a code snippet taken from SReachTools):
    sys=getChainOfIntegLtiSystem(2, delT, Polyhedron('lb',-ulim,'lb',ulim),RandomVector('Gaussian',mean_w,cov_mat_diag));    
    [Ad, Bd, Gd] = getConcatMats(sys, T);
    [~, mean_X_sans_input, cov_X_sans_input] = getHmatMeanCovForXSansInput(sys, x0, T);        
        
    %% Cost ratios b/n input and state --- scalarization term
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
        PiecewiseLinearRA
        disp(' ');
        disp(' ');
        BlackmoreTRo11PCOno08Mod


%% Plotting
%     plot_markersize = 15;
%     plot_fontSize = 20;
    plot_markersize = 10;
    plot_fontSize = 10;
    figure(1)
    clf
    hold on
    polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
        [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
    P = Polyhedron('V',polyvertex);
    h1 = plot(P,'alpha',0.1);
    h11 = scatter(1,x0(1),plot_markersize*10,'filled');
    h2 = plot(2:(T+1),xtarget,'go','MarkerSize',plot_markersize,'LineWidth',2);
    h3 = plot(2:(T+1),ono_opt_mean_X(1:2:end),'bx','LineWidth',1,'MarkerSize',plot_markersize);
    h4 = plot(2:(T+1),blackmore_opt_mean_X(1:2:end,1),'ks','MarkerSize',plot_markersize);
    h5 = plot(2:(T+1),onopwl_opt_mean_X(1:2:end),'md','LineWidth',1,'MarkerSize',plot_markersize);    
%   h4 = plot(1:T,x(1:2:end,1),'.');
%     plot(1:T,x(1:2:end,2:N),'.')  
    xlabel('time')
    ylabel('x')
    title('Trajectory')
    legend([h1 h11 h2 h3 h4 h5],{'Safe Set','Initial state','Target Trajectory',...
        sprintf('Ono2008 IRA Method (Cost: %1.3f)',ono_opt_val),sprintf('Blackmore11 PC Method (Cost: %1.3f)',blackmore_opt_val),sprintf('Piecewise linear approach (Cost: %1.3f)',onopwl_opt_val)});
    box on;
    set(gca,'FontSize',plot_fontSize)
    
%     % Figure 2 shows the cost of Ono's method with each iteration.
%     figure(2);  
%     clf
%     hold on
%     plot(ono_opt_value_array)
%     xlabel('Iteration Count')
%     ylabel('Cost, J')
%     title('Final J cost for OnoIRA2008');
    
    figure(1);
    

%% Monte-Carlo simulation using SReachTools
n_mcarlo_sims = 1e5;
% FIXME: Shouldn't be redefining this again!
hbig = kron(eye(T),h);
gbig = kron(g,[1,1])';    
xtarget_mcarlo = repmat(xtarget, 1, n_mcarlo_sims);
collection_of_input_vectors = [blackmore_opt_input_vector, ono_opt_input_vector, onopwl_opt_input_vector];
collection_of_costs = [blackmore_opt_val, ono_opt_val, onopwl_opt_val];
% Have a collection of string
collection_of_method_strings = {'BlackmoreTRO11',...
                                'OnoCDC2008    ',...
                                'PWL OnoCDC2008'};
% SReachTools for Monte-Carlo simulation
max_rel_abserror = 0.1;
fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',Delta,max_rel_abserror);
for input_vec_indx = 1:3
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
%     cost_estim = (input_state_ratio*sum(abs(U_vector))/(ulim*T) +...
%             sum(sum(abs(X_mcarlo(1:2:end,:)-xtarget_mcarlo)))/(2*g(1)*T)/n_mcarlo_sims);
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
