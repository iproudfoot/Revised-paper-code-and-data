function total_squared_error = total_square_error_12_params_vary3(fit_params,exp_data_descriptor,init_params,logic_vec)


%exp_data_descriptor format
% time series data
% concentration measured
% initial taxol concentration
% volume of medium
% initial cell number
% time of washout (if applicable)
% binary value for inclusion of nocodazole in experiment
% indices for fitting variables
% label of experiment
%Exp1a = {A1(:,1),A1(:,2),initial_taxol,V_medium,ICN,washout_time,nocod};

% Give known parameters
V_onecell=2.425*10^(-6);
B0=6.3*1000;
T0 = 3*B0;

% Define model parameters.
% 
init_params2 = init_params - logic_vec.*init_params;
init_params(init_params2 == 0) = fit_params;
k1=init_params(1);
k2=init_params(2);
k3=init_params(3)*1000;
% k4 = init_params(4);
k4=init_params(4);
a2=init_params(5);
a3 = init_params(6);
b1 = init_params(7);
b2 = init_params(8);
b3 = init_params(9)*1000;
gamma2 = init_params(10)/100;
gamma3 = init_params(11);
alpha = init_params(12)*10^(-3);
k5 = init_params(13);
k6 = init_params(14);

gamma = 1000;
Cn = 0;

a1 = a3*(T0 - B0);
k5 = 72;
k6 = 0.11;
Kbar = 781;
% d =10;
% a2 = 10;
% a3 = 15000;

%alpha = 0.0034;
gamma1 = log(2)/20;

% b2 = 2;
% b3 = 10000;
b4 = 2;
b5 = 0.00001;
%b1 = 11.52;


num_data_points = 0.0;
total_squared_error=0.0;
num_experimental_time_series_measurements=length(exp_data_descriptor);
for i=1:num_experimental_time_series_measurements
 
    
    %% Load data and metadata
    time=cell2mat(exp_data_descriptor{i}(1));
    exp_conc_measurement=cell2mat(exp_data_descriptor{i}(2));

    initial_condition=cell2mat(exp_data_descriptor{i}(3));
    k_cellnumber=-.0053+0.04*exp(-0.02*initial_condition);
    V_medium=cell2mat(exp_data_descriptor{i}(4));
    ICN=cell2mat(exp_data_descriptor{i}(5));
    taxol_washout_time=cell2mat(exp_data_descriptor{i}(6));
    nocadozole_binary=cell2mat(exp_data_descriptor{i}(7));
    label = cell2mat(exp_data_descriptor{i}(9));
    
%     %% Inserting Gaussian noise
%     Mean = mean(exp_conc_measurement);
%     n3 = length(exp_conc_measurement);
%     R = Mean/10*randn(n3,1);
%     exp_conc_measurement = exp_conc_measurement + R;
    
    model_variable_fit_ind=cell2mat(exp_data_descriptor{i}(8)); % indices for solution vector corresponding to experimental data
    
    
    
    % solve model
    %% Generating the solution for the particular iteration
    %% if loop for nocodozole
       if nocadozole_binary == 0
    %% if loop for a washout experiment
    
    if taxol_washout_time < 1000
        y_initl = [ICN 0 0 0 initial_condition 0 0 0];
        
        sol2 =  ode15s(@model1T1,[0 taxol_washout_time],y_initl,[],k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3,k_cellnumber,gamma,Cn,V_medium);
        
        g = deval(sol2,taxol_washout_time); % getting variables at time of washout
        
        y_initl = [g(1) g(2) g(3) g(4) 0 0 0 0]; % updated initial conditions
    else
        y_initl = [ICN 0 0 0 initial_condition 0 0 0]; %ICs with no washout
    end

    sol = ode15s(@model1T1,[0 time(end)],y_initl,[],k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3,k_cellnumber,gamma,Cn,V_medium);
    
    sol1 = deval(sol,time);
    end
    if nocadozole_binary == 1
        B0 = 0;
        k1 =0;
          Cn = 1000;
          gamma = 100;
           %% if loop for a washout experiment
    
    if taxol_washout_time < 1000
        y_initl = [ICN 0 0 0 initial_condition 0 T0 B0];
        
        sol2 =  ode15s(@model1T1,[0 taxol_washout_time],y_initl,[],k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3,k_cellnumber,gamma,Cn,V_medium);
        
        g = deval(sol2,taxol_washout_time); % getting variables at time of washout
        
        y_initl = [g(1) g(2) g(3) g(4) 0 0 g(7) g(8)]; % updated initial conditions
    else
        y_initl = [ICN 0 0 0 initial_condition 0 T0 B0]; %ICs with no washout
    end

    sol = ode15s(@model1T1,[0 time(end)],y_initl,[],k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3,k_cellnumber,gamma,Cn,V_medium);
    
    sol1 = deval(sol,time);
    k1 = init_params(1);
    B0 = 6300;
    end

    % get total concentration from model
    
    if length(model_variable_fit_ind) == 1
        model_measureable_quantity = sol1(model_variable_fit_ind,:);
    else
        model_measureable_quantity=sum(sol1(model_variable_fit_ind,:));
    end
%m2 = max([exp_conc_measurement' model_measureable_quantity]);
m = mean(model_measureable_quantity);

% m2 = max(exp_conc_measurement');
% m1 = min(exp_conc_measurement');

% error = sum(((model_measureable_quantity-exp_conc_measurement')./(m2-m1)).^2);

     % compute error
    w = ones(1,length(exp_conc_measurement'));
    out1 = (model_measureable_quantity-exp_conc_measurement');
    out2 = (model_measureable_quantity-m);% percentage error of model to each data point
    % out(exp_conc_measurement'==0)=0; %avoid dividing by 0 - gives NaN
    %out;
    error = sqrt(sum(out1.^2)/sum(out2.^2)); %mean percentage error of the model to the data for given experiment
%      
    % compute total error
    num_data_points = num_data_points + length(exp_conc_measurement');
    total_squared_error = total_squared_error+error;
end