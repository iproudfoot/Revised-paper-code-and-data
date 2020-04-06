clear all
close all

%% Fitting to Adrian's Data
% mydir = '/Users/idmproudfoot/Box Sync/Adrian''sDataForSharing/lab_data/'; % Specify directory
% save_file = '/Users/idmproudfoot/Documents/fig9/'; % Place to save figures

mydir = 'C:\Users\isaac\Documents\MATLAB\lab_data\'; % Specify directory
save_file = 'C:\Users\isaac\Documents\neat_figs2\'; % Place to save figures

%%
V_onecell=2.425*10^(-6);
B0=6.93*1000;
T0=100*B0/33;
Bmaxm = 3940;
Kbar = 781;
Tmax = 4*B0;
k_cellnumber = 0;
%% Loading the Data
%% Dec 2016

file_name = [mydir 'TaxolDataDec2016.xlsx'];

% with/without nocadozole data
initial_taxol = 100;
V_medium = 1*10^3;
ICN = 1*10^6;
washout_time = 1000;
A = xlsread(file_name,1); 

A1 = [A(:,1) A(:,2)]; % without nocod
nocod = 0;
fit_inds = [2 3 4];

Exp1a = {A1(:,1)/60,20*A1(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

A2 = [A(:,1)/60 A(:,3)]; % with nocod
nocod = 1;
fit_inds = [2 3 4];

Exp1b = {A2(:,1), A2(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

% washout expt
B = xlsread(file_name,2);
initial_taxol = 100;
V_medium = 1*10^3;
ICN = 1*10^6;
nocod = 0; 
washout_time = 2;
fit_inds = [2 3 4];

Exp2 = {B(:,1)/60,B(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds}; 

% effect of incubation time (odd values for low concs)
C = xlsread(file_name,3);
t_data = [0 2 16];
initial_concs = [0.3 0.6 1.2 2.5 5 10 20 40 80 100];
initial_taxol = 100;
V_medium = 1*10^3;
ICN = 1*10^6;
washout_time = 1000;
nocod = 0;
fit_inds = [2 3 4];

Exp3a = {t_data',[0 C(1,2:3)]',0.3,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3b = {t_data',[0 C(2,2:3)]',0.6,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3c = {t_data',[0 C(3,2:3)]',1.2,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3d = {t_data',[0 C(4,2:3)]',2.5,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3e = {t_data',[0 C(5,2:3)]',5,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3f = {t_data',[0 C(5,2:3)]',10,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3g = {t_data',[0 C(6,2:3)]',20,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3h = {t_data',[0 C(7,2:3)]',40,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3i = {t_data',[0 C(8,2:3)]',80,V_medium,ICN,washout_time,nocod,fit_inds};
Exp3j = {t_data',[0 C(9,2:3)]',100,V_medium,ICN,washout_time,nocod,fit_inds};

%% Feb 2017

file_name = [mydir 'TaxolDataFeb2017'];

% Exp 1
D = xlsread(file_name,1);
ICN = 0.7*10^6;
V_medium = 1*10^3;
washout_time = 2;
initial_taxol = 100;
fit_inds = [2 3 4];

D1 = [D(:,1) D(:,2)];
D2 = [D(:,1) D(:,3)];
D3 = [D(:,1) D(:,4)];

nocod = 0;
Exp4a = {D1(:,1)/60,D1(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

initial_taxol = 100;
nocod = 1;
Exp4b = {D2(:,1)/60,D2(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

initial_taxol = 100;
nocod = 1;
Exp4c = {D3(:,1)/60,D3(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

% Exp 2
E = xlsread(file_name,2);
t_data = [0 4 8]';
ICN = 0.7*10^6;
V_medium = 1*10^3;
washout_time = 2;
nocod = 0;
fit_inds = [2 3 4];

Exp5a = {t_data,E(1,2:4)'*1000,800,V_medium,ICN,washout_time,nocod,fit_inds};
Exp5b = {t_data,E(2,2:4)'*1000,400,V_medium,ICN,washout_time,nocod,fit_inds};
Exp5c = {t_data,E(3,2:4)'*1000,200,V_medium,ICN,washout_time,nocod,fit_inds};
Exp5d = {t_data,E(4,2:4)'*1000,100,V_medium,ICN,washout_time,nocod,fit_inds};
Exp5e = {t_data,E(5,2:4)'*1000,50,V_medium,ICN,washout_time,nocod,fit_inds};
Exp5f = {t_data,E(6,2:4)'*1000,25,V_medium,ICN,washout_time,nocod,fit_inds};

ICN = 0.8*10^6;
Exp5g = {t_data,E(7,2:4)'*1000,12.5,V_medium,ICN,washout_time,nocod,fit_inds};

ICN = 0.7*10^6;
Exp5h = {t_data,E(8,2:4)'*1000,6.25,V_medium,ICN,washout_time,nocod,fit_inds};

% Exp 3
F = xlsread(file_name,3);
t_data = [0 2]';
ICN = 0.7*10^6;
V_medium = 1*10^3;
washout_time = 1000;
nocod = 0;
fit_inds = [2 3 4];

Exp6a = {t_data,[0 F(1,2)*1000]',40,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6b = {t_data,[0 F(2,2)*1000]',30,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6c = {t_data,[0 F(3,2)*1000]',25,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6d = {t_data,[0 F(4,2)*1000]',20,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6e = {t_data,[0 F(5,2)*1000]',10,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6f = {t_data,[0 F(6,2)*1000]',7.5,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6g = {t_data,[0 F(7,2)*1000]',6.25,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6h = {t_data,[0 F(8,2)*1000]',5,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6i = {t_data,[0 F(9,2)*1000]',3.75,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6j = {t_data,[0 F(10,2)*1000]',2.5,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6k = {t_data,[0 F(11,2)*1000]',1.8,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6l = {t_data,[0 F(12,2)*1000]',1.25,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6m = {t_data,[0 F(13,2)*1000]',0.9,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6n = {t_data,[0 F(14,2)*1000]',0.6,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6o = {t_data,[0 F(15,2)*1000]',0.45,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6p = {t_data,[0 F(16,2)*1000]',0.3,V_medium,ICN,washout_time,nocod,fit_inds};
Exp6q = {t_data,[0 F(17,2)*1000]',0.225,V_medium,ICN,washout_time,nocod,fit_inds};

%% May 2017

file_name = [mydir 'TaxolDataMay2017'];

% Exp 1
G = xlsread(file_name,1);
initial_taxol = 1000;
washout_time = 2;
ICN = 0.5*10^6;
V_medium = 1*10^3;
nocod = 0;
fit_inds = [2 3 4];

Exp7 = {G(:,1)/60,G(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

% Exp 2.1
H = xlsread(file_name,2);
initial_taxol = 100;
washout_time = 2;
ICN = 0.5*10^6;
V_medium = 1*10^3;
nocod = 0;
fit_inds = [2 3 4];

Exp8 = {H(:,1)/60,H(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

% Exp 2.2
I = xlsread(file_name,3);
I(1,3) = I(9,2);
I(9,3) = 0.006*10^3;
initial_taxol = 100;
ICN = 0.5*10^6;
V_medium = 1*10^3;
nocod = 0;
fit_inds = [2 3 4];
washout_time = 2; 

Exp9a = {I(:,1)/60,I(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

washout_time = 0;
Exp9b = {I(:,1)/60,I(:,3)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

% Exp 3
J = xlsread(file_name,4);
ICN = 0.5*10^6;
V_medium = 1*10^3;
t_data = [0 8];
washout_time = 2;
nocod = 0;
fit_inds = [2 3 4];

Exp10a = {t_data,J(1,2:3)'*1000/1.4,3200,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10b = {t_data,J(2,2:3)'*1000/1.4,1600,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10c = {t_data,J(3,2:3)'*1000/1.4,800,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10d = {t_data,J(4,2:3)'*1000/1.4,400,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10e = {t_data,J(5,2:3)'*1000/1,4,200,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10f = {t_data,J(6,2:3)'*1000/1.4,100,V_medium,ICN,washout_time,nocod,fit_inds};
Exp10g = {t_data,J(7,2:3)'*1000/1.4,50,V_medium,ICN,washout_time,nocod,fit_inds};

%% June 2017
file_name = [mydir 'TaxolDataJun17'];

% Exp 1
K = xlsread(file_name,1);

V_medium = 1*10^3;
ICN = 0.9*10^6;
nocod = 0;
washout_time = 1000;
fit_inds = [2 3 4];

Exp11a = {[0 2],[0 K(2,2)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11b = {[0 2],[0 K(3,2)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11c = {[0 2],[0 K(4,2)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11d = {[0 2],[0 K(5,2)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11e = {[0 2],[0 K(6,2)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11f = {[0 2],[0 K(7,2)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds};

V_medium = 2*10^3;
Exp11g = {[0 2],[0 K(2,3)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11h = {[0 2],[0 K(3,3)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11i = {[0 2],[0 K(4,3)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11j = {[0 2],[0 K(5,3)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11k = {[0 2],[0 K(6,3)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11l = {[0 2],[0 K(7,3)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds};

V_medium = 10^10^3;
Exp11m = {[0 2],[0 K(2,4)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11n = {[0 2],[0 K(3,4)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11o = {[0 2],[0 K(4,4)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11p = {[0 2],[0 K(5,4)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11q = {[0 2],[0 K(6,4)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp11r = {[0 2],[0 K(7,4)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds};

%% Sept 2017

file_name = [mydir 'TaxolDataSept2017'];

% Exp 1
L = xlsread(file_name,1);
V_medium = 10*10^3;
ICN = 1.1*10^6;
nocod = 0;
washout_time = 1000;
initial_taxol = 100;
fit_inds = [2 3 4];

Exp12 = {L(:,1),L(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};

%% Feb 2018

% Exp 1 (no cells)
file_name = [mydir 'TaxolDataFeb2018.xlsx'];
M = xlsread(file_name,1);
V_medium = 1*10^3;
initial_taxol = 150;
ICN = 0;
nocod = 0;
washout_time = 1000;
fit_inds = [5]; 

Exp13 = {M(:,1),M(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds};


%% jordan data 
int_conc = [3 10 100 1000];
concs = [1.8 4.8 40.5 111];
N = [6.3 8.4 19.1 28.8]*1000;
ICN = 0.27*10^6;
V_medium = 1*10^3;
nocod = 0;
washout_time = 1000;

fit_inds = [8];

Exp14a = {[0 20]',[B0 N(1)]',int_conc(1),V_medium,ICN,washout_time,nocod,fit_inds};
Exp14b = {[0 20]',[B0 N(2)]',int_conc(2),V_medium,ICN,washout_time,nocod,fit_inds};
Exp14c = {[0 20]',[B0 N(3)]',int_conc(3),V_medium,ICN,washout_time,nocod,fit_inds};
Exp14d = {[0 20]',[B0 N(4)]',int_conc(4),V_medium,ICN,washout_time,nocod,fit_inds};

%%
%% Model parameters
n = 20;

k1 = 4*ones(n,1); 
k2 = 16*ones(n,1);
k3 = 350*ones(n,1);
k4 = 10*ones(n,1);
k5 = 0.082*ones(n,1);
k6 = 0.106*ones(n,1);
alpha = 2.7*10^(-3)*ones(n,1);
a3 = 0.8*12*ones(n,1) + 0.4*12*rand(n,1);
a2 = zeros(n,1);
a1 = a3*(T0-B0);
V_medium = 1000;
V_onecell = 2.25*10^(-6)*ones(n,1);

b2 =0.0023*ones(n,1);
b3 = 0.01*ones(n,1);
b1 = 0.8*21.6*ones(n,1) + 0.4*21.6*rand(n,1);

% k1 = 4;
% k2 = 16;
% k3 = 700;
% k4 = 14;
% k5 = 72;
% k6 = 0.11;
% alpha = 2.7*10^(-3);
% a2 = 0;
% a3 = 12;
% a1 = a3*(T0-B0);
% V_medium = 1000;
% V_onecell = 2.25*10^(-6);
% 
% b2 =0.0023;
% b3 = 0.01;
% b1 = 21.6;
% gamma = 1000000000;
% Cn = 0;
% alpha2 = 0;
gamma1 = log(2)/24*ones(n,1);
gamma2 = gamma1;
gamma3 = 6000*ones(n,1);
Kbar = 781*ones(n,1);

k_cellnumber = 0;

param_sets = zeros(n,18);
param_sets = [k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3];
P1 = num2cell(param_sets);

%% Plot uptake profile
ICN = 0.7*10^6;
y_0 = [ICN 0 0 0 100 0 T0 B0];
t = 0:0.01:24;
h = zeros(8,length(t),n);
Cn = 0;
gamma = 0;
for i = 1:n
sol = ode15s(@model1T3, [0 24], y_0, [], P1{i,:}, k_cellnumber,gamma,Cn,V_medium);
h(:,:,i) = deval(sol,t);
end


figure(2)
plot(A1(:,1)/60,25*A1(:,2)/1000,'go')
hold on
% plot(A2(1:8,1)+1,30*A2(1:8,2)/1000,'bx')
% hold on
h1 = (h(2,:,:) + h(3,:,:) + h(4,:,:))/1000;
h1 = squeeze(h1);
plot_distribution(t, h1');
hold on
plot(L(:,1),1.2*L(:,2),'ro')
%hold on 
%p2 =plot(t+1, (h1(2,:) + h1(3,:) + h1(4,:))/1000,'k--','linewidth',1.5);

xlabel('time of incubation (h)','fontsize',16)
ylabel('intracellular concentration (\mu M)','fontsize',16)

%title('time series data for incubation (100 nM)','fontsize',16)
file_name = [save_file 'pdeuptake_vary.png'];
export_fig(figure(2), file_name);

%%
model_vol = zeros(3,21,n);
V = [1 2 10]*1000;
concs = 0:10:200;
for j = 1:3
    for i = 1:21 
        V_medium = V(j);
        y_0 = [0.7*10^6 0 0 0 concs(i) 0 T0 B0];
        for k = 1:n
        sol = ode15s(@model1T3, [0 2], y_0, [], P1{k,:},k_cellnumber,gamma,Cn,V_medium);
        g = deval(sol,2);
        model_vol(j,i,k) = (g(2)+g(3)+g(4))/1000;
        end
    end
end

m1 = squeeze(model_vol(1,:,:));
m2 = squeeze(model_vol(2,:,:));
m3 = squeeze(model_vol(3,:,:));

figure(3)
plot_distribution(concs,m1','color',[1 0 0])
hold on
plot_distribution(concs,1.15*m2','color',[0 0 1])
hold on
plot_distribution(concs,1.3*m3','color',[0 0 0])
hold on
plot(K(:,1),K(:,2),'rx')
hold on
plot(K(:,1),K(:,3),'bx')
hold on
plot(K(:,1),K(:,4),'kx')
xlabel('initial concentration (nM)','fontsize',16)
ylabel('intracellular concentration (\mu M)','fontsize',16)
%legend({'1ml volume','2ml volume','10ml volume'},'location','southeast')
%title('volume variation data','fontsize',16)
file_name = [save_file 'pdevol_vary.png'];
export_fig(figure(3), file_name);

%%  Concentration after 2 hrs

Feb_concs = fliplr(E(:,1)'/1000);
May_concs = fliplr(J(:,1)');
Feb_washout_concs = fliplr(E(:,2)');
May_washout_concs = fliplr(J(:,2)');
Feb_eq_concs = fliplr(E(:,4)');
May_eq_concs = fliplr(J(:,3)');
washout_mod = zeros(n,31);
eq_mod = zeros(n,31);
concs = logspace(-4, 1, 31);
for i = 1:31
    for k = 1:n
        ICN = 0.7*10^6;
        y_0 = [ICN 0 0 0 1000*concs(i) 0 T0 B0];
        sol = ode15s(@model1T3, [0 2], y_0, [], P1{k,:}, k_cellnumber,gamma,Cn,V_medium);
        h = deval(sol,2);
        washout_mod(k,i) = (h(2)+h(3)+h(4))/1000;
        y_1 = [h(1) h(2) h(3) h(4) 0 0 h(7) h(8)];
        sol1 = ode15s(@model1T3, [0 8], y_1, [], P1{k,:}, k_cellnumber,gamma,Cn,V_medium);
        g = deval(sol1, 8);
        eq_mod(k,i) = (g(2)+g(3)+g(4))/1000;
    end
end

figure(5)
f1 = loglog(C(5:end,1),20*C(5:end,2)/1000,'gx');
hold on
f2 = loglog(Feb_concs*1000, Feb_washout_concs, 'rx');

hold on
f3 = plot(F(:,1),F(:,2),'rx');

hold on
f4 = plot(May_concs*1000,May_washout_concs/1.4, 'bx');

hold on
f5 = plot(K(:,1),K(:,2),'kx');

hold on
plot_distribution(1000*concs, washout_mod/1.2);
hold on
plot_distribution(1000*concs, 2.15*eq_mod);
hold on
f6 = plot(1000*Feb_concs, Feb_eq_concs, 'ro');

hold on
f7 = plot(1000*May_concs,May_eq_concs/1.4, 'bo');

%legend([d1,d2],{'Concentration after 2 h incubation','Concentration 8 h after washout'},'location','northwest')
xlabel('initial concentration (nM)','fontsize',16)
ylabel('intracellular concentration (\mu M)','fontsize',16)
%title('intracellular concentrations after 2 h incubation','fontsize',16)
file_name = [save_file 'pdeintra_vary.png'];
export_fig(figure(5), file_name);

%% 100 nM washout
t1 = 0:0.01:2;
t = 0.01:0.01:24;
tot_concpre =zeros(n,length(t1));
tot_conc = zeros(n,length(t));
%alpha = 1.9*10^(-3);
%k5 = 16;
%k1 =3.5;
%gamma = 10000;
% %Kbar = 80;
% param_sets(1,:) = [k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3];
% P1 = num2cell(param_sets(1,:));
V_medium = 1000;
% Cn =0;
ICN = 0.5*10^6;
y_0 = [ICN 0 0 0 100 0 T0 B0];
for k = 1:n
sol = ode15s(@model1T3, [0 2], y_0, [], P1{k,:}, k_cellnumber,gamma,Cn,V_medium);
h = deval(sol,2);

h1 = deval(sol,t1);
tot_concpre(k,:) = (h1(2,:) + h1(3,:) + h1(4,:))/1000;
y_1 = [h(1) h(2) h(3) h(4) 0 0 h(7) h(8)];
sol1 = ode15s(@model1T3, [0 24], y_1, [], P1{k,:}, k_cellnumber,gamma,Cn,V_medium);

v = deval(sol1, t);
tot_conc(k,:) = (v(2,:) + v(3,:) + v(4,:))/1000;

end
figure(6)
plot_distribution([t1 t+2],[tot_concpre, tot_conc])
hold on
% plot([5*t2 t+2],1.25*[tot_conc1pre tot_conc1],'b','linewidth',1.5)
% hold on
% plot([t1 t+2],[tot_concpre tot_conc2-2],'k','linewidth',1.5)
% hold on
plot(2+D1(:,1)/60,D1(:,2),'rx','Markersize',10)
% hold on
% plot(2+D2(:,1)/60,D2(:,2),'bx','Markersize',10)
% hold on
% plot(2+D3(:,1)/60,D3(:,2),'kx','Markersize',10)
% hold on
plot(2+D1(:,1)/60,D1(:,2),'ro','Markersize',10)
hold on
plot(2+H(:,1)/60,H(:,2),'r^','Markersize',10)
hold on
plot(2+I(1:9,1)/60,I(1:9,2)/1.4,'rd','Markersize',10)
xlabel('time (h)','fontsize',16)
ylabel('intracellular concentration (\mu M)','fontsize',16)
file_name = [save_file 'PDEwashtimeseries_vary.png'];
export_fig(figure(6), file_name);