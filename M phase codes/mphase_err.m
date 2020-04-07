function err = mphase_err(x);
%%
fsize = 16;
lwidth = 2;
%% Give known parameters (in nM and muL)

%% Model parameters
B0 = 6300;
T0 = 3*B0;
k_cellnumber = 0;
k1 = 4;
k2 = 16;
k3 = 700;
k4 = 14;
k5 = 72;
k6 = 0.11;
alpha = 2.7*10^(-3);
a2 = 0;
a3 = 12;
a1 = a3*(T0-B0);
V_medium = 1000;
V_onecell = 2.25*10^(-6);

b2 =0.0023;
b3 = 0.01;
b1 = 21.6;
gamma = 1000000000;
Cn = 0;
alpha2 = 0;
gamma1 = log(2)/24;
gamma2 = gamma1;
gamma3 = 6000;
Kbar = 781;

k_cellnumber = 0;

k7 = 0.035;
u1 = 2.2;
% u2 = 0.41;
% v1 = 0.0071;
v2 = 0.0001;

v1 = x(1);
u2 = x(2);
% v2 = x(3);

param_sets(1,:) = [k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3,k7,u1,u2,v1,v2];
P1 = num2cell(param_sets(1,:));

err = 0;

save_file = 'C:\Users\isaac\Documents\neat_figs3\';
%% Formatting
set(groot,'defaultLineLineWidth',1.5)
%% Loading m phase data
mydir = 'C:\Users\isaac\Documents\MATLAB\lab_data';
file_name = [mydir '\GV0433_TaxDoses_summary.xlsx'];
B = xlsread(file_name,1,'A2:BW51');
num_experiments = 10;
t_stop = 7;
pc_exited = zeros(1,7);
concs = [80 40 20 10 5 2.5 1.25 0.62 0.31 100];
n1 = 24*60;

for i = 1:num_experiments
    entry_time = B(:,8*(i-1)+1)/15;
    exit_time = B(:,8*(i-1)+2)/15;
    arrest_time = B(:,8*(i-1)+3)/15;
    
    perma_arrest_entry = entry_time(arrest_time>= t_stop);
    arrest_time(arrest_time >= t_stop) = NaN;
x = entry_time(entry_time <= 17);
y = arrest_time(entry_time <= 17);

%% CDF of arrested cells
str = sprintf('%g',i);

full_arrest_pred_entry = perma_arrest_entry;
num_arrested(i) = sum(perma_arrest_entry./perma_arrest_entry);
sf = 100*sum(isnan(arrest_time))/50;
[xdata, ydata] = cdfplot_scaled_data(full_arrest_pred_entry,sf);

%% Simulate
ICN = 0.5*10^6;
T = 19; t = 2:0.01:T;
t1 = 0:0.01:17;
init_conc = concs(i);
y_0 = [ICN 0 0 0 init_conc 0 T0 B0 0.95*ICN 0.05*ICN/3 0.05*ICN/3 0.05*ICN/3 0];
sol = ode15s(@model1T3mphase,[0 2], y_0, [], P1{:}, k_cellnumber,gamma,Cn,V_medium);
g = deval(sol,2);
y_1 = g';

sol1 = ode15s(@model1T3mphase,[0 T], y_1, [], P1{:}, k_cellnumber,gamma,Cn,V_medium);
g = deval(sol1,xdata);

G1_S_G2_phase = g(9,:);
M1 = g(10,:);
M2 = g(11,:);
M3 = g(12,:);
R = g(13,:);

entries = 0.035*G1_S_G2_phase;
avg_entries = mean(entries);
tot_entries = 17*avg_entries;
trap_ents = trapz(entries);
removed_cells = g(13,:);
removed_cells = removed_cells - removed_cells(1);
pc_removed = 100*removed_cells/tot_entries; 

tot_sq_err = sum((pc_removed - ydata').^2)/length(ydata');
err = err + tot_sq_err;
end
end