clear all
close all

%% Fitting to Adrian's Data
mydir = 'C:\Users\isaac\Documents\MATLAB\lab_data\'; % Specify directory
save_file = 'C:\Users\isaac\Documents\neat_figs3\'; % Place to save figures

%%
global exp_data_descriptor initial_guess logic_vec fixed_params

%%
V_onecell=2.425*10^(-6);
B0=6.93*1000;
T0=100*B0/33;
Bmaxm = 3940;
Kbar = 781;
Tmax = 4*B0;
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
label = 'Dec 16 Exp 1 (no noc) - Exp1a';
Exp1a = {A1(:,1)/60,30*A1(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

A2 = [A(:,1)/60 A(:,3)]; % with nocod
nocod = 1;
fit_inds = [2 3 4];
label = 'Dec 16 Exp 1 (with noc) - Exp1b';

Exp1b = {A2(:,1), 30*A2(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

% washout expt
B = xlsread(file_name,2);
initial_taxol = 100;
V_medium = 1*10^3;
ICN = 1*10^6;
nocod = 0; 
washout_time = 2;
fit_inds = [2 3 4];
label = 'Dec 16 Exp 2 washout - Exp2';

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
label = 'Dec 16 Exp 3 Uptake on conc - Exp3a-3j';

Exp3a = {t_data',[0 C(1,2:3)]',0.3,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3b = {t_data',[0 C(2,2:3)]',0.6,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3c = {t_data',[0 C(3,2:3)]',1.2,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3d = {t_data',[0 C(4,2:3)]',2.5,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3e = {t_data',[0 C(5,2:3)]',5,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3f = {t_data',[0 C(5,2:3)]',10,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3g = {t_data',[0 C(6,2:3)]',20,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3h = {t_data',[0 C(7,2:3)]',40,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3i = {t_data',[0 C(8,2:3)]',80,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp3j = {t_data',[0 C(9,2:3)]',100,V_medium,ICN,washout_time,nocod,fit_inds,label};

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
label = 'Feb 17 Exp 1.1 (washout) - Exp 4a';
Exp4a = {D1(:,1)/60,D1(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

initial_taxol = 100;
nocod = 1;
label = 'Feb 17 Exp 1.2 (incub + wash with noc) - Exp 4b';
Exp4b = {D2(:,1)/60,D2(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

initial_taxol = 100;
nocod = 1;
label = 'Feb 17 Exp 1.3 (washout into noc) - Exp 4c';
Exp4c = {D3(:,1)/60,D3(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

% Exp 2
E = xlsread(file_name,2);
t_data = [0 4 8]';
ICN = 0.7*10^6;
V_medium = 1*10^3;
washout_time = 2;
nocod = 0;
fit_inds = [2 4];
label = 'Feb 17 Exp 2 (washouts) - Exp 5a-5g';

Exp5a = {t_data,E(1,2:4)'*1000,800,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp5b = {t_data,E(2,2:4)'*1000,400,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp5c = {t_data,E(3,2:4)'*1000,200,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp5d = {t_data,E(4,2:4)'*1000,100,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp5e = {t_data,E(5,2:4)'*1000,50,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp5f = {t_data,E(6,2:4)'*1000,25,V_medium,ICN,washout_time,nocod,fit_inds,label};

ICN = 0.8*10^6;
Exp5g = {t_data,E(7,2:4)'*1000,12.5,V_medium,ICN,washout_time,nocod,fit_inds,label};

ICN = 0.7*10^6;
Exp5h = {t_data,E(8,2:4)'*1000,6.25,V_medium,ICN,washout_time,nocod,fit_inds,label};

% Exp 3
F = xlsread(file_name,3);
t_data = [0 2]';
ICN = 0.7*10^6;
V_medium = 1*10^3;
washout_time = 1000;
nocod = 0;
fit_inds = [2 3 4];
label = 'Feb 2017 Exp 2 (incubation for diff. concs) - Exp 6a-6q';

Exp6a = {t_data,[0 F(1,2)*1000]',40,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6b = {t_data,[0 F(2,2)*1000]',30,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6c = {t_data,[0 F(3,2)*1000]',25,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6d = {t_data,[0 F(4,2)*1000]',20,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6e = {t_data,[0 F(5,2)*1000]',10,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6f = {t_data,[0 F(6,2)*1000]',7.5,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6g = {t_data,[0 F(7,2)*1000]',6.25,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6h = {t_data,[0 F(8,2)*1000]',5,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6i = {t_data,[0 F(9,2)*1000]',3.75,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6j = {t_data,[0 F(10,2)*1000]',2.5,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6k = {t_data,[0 F(11,2)*1000]',1.8,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6l = {t_data,[0 F(12,2)*1000]',1.25,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6m = {t_data,[0 F(13,2)*1000]',0.9,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6n = {t_data,[0 F(14,2)*1000]',0.6,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6o = {t_data,[0 F(15,2)*1000]',0.45,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6p = {t_data,[0 F(16,2)*1000]',0.3,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp6q = {t_data,[0 F(17,2)*1000]',0.225,V_medium,ICN,washout_time,nocod,fit_inds,label};

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
label = 'May 2017 Experiment 1 (washout at 1000nM) - Exp 7';
Exp7 = {G(:,1)/60,G(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

% Exp 2.1
H = xlsread(file_name,2);
initial_taxol = 100;
washout_time = 2;
ICN = 0.5*10^6;
V_medium = 1*10^3;
nocod = 0;
fit_inds = [2 3 4];
label = 'May 2017 Experiment 2.1 (washout at 100nM) - Exp 8';

Exp8 = {H(:,1)/60,H(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

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
label = 'May 2017 Experiment 2.3 (1st washout at 100nM) - Exp 9a';
Exp9a = {I(1:9,1)/60,I(1:9,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

label = 'May 2017 Experiment 2.2 (2nd washout at 100nM) - Exp 9b';
washout_time = 0;
Exp9b = {I(:,1)/60,I(:,3)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};
% this cant be used in its current form

% Exp 3
J = xlsread(file_name,4);
ICN = 0.5*10^6;
V_medium = 1*10^3;
t_data = [0 8];
washout_time = 2;
nocod = 0;
fit_inds = [2 3 4];
label = 'May 2017 Experiment 3 (washouts) - Exp 10a-10g';

Exp10a = {t_data,J(1,2:3)'*1000/1.4,3200,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10b = {t_data,J(2,2:3)'*1000/1.4,1600,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10c = {t_data,J(3,2:3)'*1000/1.4,800,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10d = {t_data,J(4,2:3)'*1000/1.4,400,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10e = {t_data,J(5,2:3)'*1000/1.4,200,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10f = {t_data,J(6,2:3)'*1000/1.4,100,V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp10g = {t_data,J(7,2:3)'*1000/1.4,50,V_medium,ICN,washout_time,nocod,fit_inds,label};

%% June 2017
file_name = [mydir 'TaxolDataJun17'];

% Exp 1
K = xlsread(file_name,1);

V_medium = 1*10^3;
ICN = 0.9*10^6;
nocod = 0;
washout_time = 1000;
fit_inds = [2 3 4];
label = 'June 17 Exp 1 (incubations at 1ml) - Exp11a-g';
Exp11a = {[0 2],[0 K(1,2)]'*1000,K(1,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11b = {[0 2],[0 K(2,2)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11c = {[0 2],[0 K(3,2)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11d = {[0 2],[0 K(4,2)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11e = {[0 2],[0 K(5,2)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11f = {[0 2],[0 K(6,2)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11g = {[0 2],[0 K(7,2)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds,label};

label = 'June 17 Exp 1 (incubations at 2ml) - Exp11h-n';
V_medium = 2*10^3;
Exp11h = {[0 2],[0 K(1,3)]'*1000,K(1,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11i = {[0 2],[0 K(2,3)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11j = {[0 2],[0 K(3,3)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11k = {[0 2],[0 K(4,3)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11l = {[0 2],[0 K(5,3)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11m = {[0 2],[0 K(6,3)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11n = {[0 2],[0 K(7,3)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds,label};

label = 'June 17 Exp 1 (incubations at 2ml) - Exp11o-u';
V_medium = 10^10^3;
Exp11o = {[0 2],[0 K(1,4)]'*1000,K(1,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11p = {[0 2],[0 K(2,4)]'*1000,K(2,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11q = {[0 2],[0 K(3,4)]'*1000,K(3,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11r = {[0 2],[0 K(4,4)]'*1000,K(4,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11s = {[0 2],[0 K(5,4)]'*1000,K(5,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11t = {[0 2],[0 K(6,4)]'*1000,K(6,1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp11u = {[0 2],[0 K(7,4)]'*1000,K(7,1),V_medium,ICN,washout_time,nocod,fit_inds,label};

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
label = 'Sept 17 Exp 1 incubation at 100 nM - Exp 12'; 

Exp12 = {L(:,1),L(:,2)*1000,initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};

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
label = 'Feb 2018/Sept 17 (extracellular conc, no cells) - Exp 13'; 

Exp13 = {M(:,1),M(:,2),initial_taxol,V_medium,ICN,washout_time,nocod,fit_inds,label};


%% jordan data 
int_concs = [3 10 100 1000];
concs = [1.8 4.8 40.5 111];
N = [6.3 8.4 19.1 28.8]*1000;
ICN = 0.26*10^6;
V_medium = 1*10^3;
nocod = 0;
washout_time = 1000;
label = 'Jordan MT data - Exp 14a-d';

fit_inds = [8];

Exp14a = {[0 20]',[B0 N(1)]',int_concs(1),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp14b = {[0 20]',[B0 N(2)]',int_concs(2),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp14c = {[0 20]',[B0 N(3)]',int_concs(3),V_medium,ICN,washout_time,nocod,fit_inds,label};
Exp14d = {[0 20]',[B0 N(4)]',int_concs(4),V_medium,ICN,washout_time,nocod,fit_inds,label};

% %% for loop for minimisation test
% 
Kbar = 781;


k1 = 2.9676;
k2 = 15;
k3 = 0.7341;
k4 = 14;
k5 = 72;
k6 = 0.11;
alpha = 0.8652;
a2 = 11.8614;
a3 = 21.5105;
a1 = a3*(T0-B0);
V_medium = 1000;
V_onecell = 2.25*10^(-6);

b2 =14.209;
b3 = 0.27;
b1 = 11.52;
gamma = 1000000000;
Cn = 0;
alpha2 = 0;
gamma1 = log(2)/24;
gamma2 = 0.0059;
gamma3 = 1.17*10^3;
Kbar = 781;

initial_guess = [k1 k2 k3 k4 a2 a3 b1 b2 b3 gamma2 gamma3 alpha k5 k6];
l = length(initial_guess);
scaling_factors = ones(l,l);
scaling_factors(3) = 1000;
scaling_factors(9) = 1000;
scaling_factors(10) = 0.01;
scaling_factors(12) = 10^(-3);

exp_data_descriptor = {Exp4a,Exp7,Exp14b};
%exp_data_descriptor = {Exp4b};
logic_vec = ones(1,length(initial_guess));
% exp_data_descriptor = {Exp4a,Exp4b,Exp5g,Exp5h,Exp6a,Exp6b,Exp6c,Exp6d,Exp7,Exp8};
logic_vec = [1 0 1 0 1 0 0 0 0 0 0 1 1 1];
%%
%% Plotting error function landscapes of parameters pairwise
fit_params = [3 0.72 11.73 0.85];
% str4 = str3(~~logic_vec);
str4 = ["$$ \lambda_{2} $$", "$$ \lambda_{3} $$", "$$ \lambda_{9} $$", "$$ \lambda_{1} $$" "$$\eta_{7}$$" "$$ \eta_{8} $$"];
scaling_factors = ones(l,12);
scaling_factors(1) = 0.5*10^2;
scaling_factors(2) = 1/800;
scaling_factors(3) = 100;
scaling_factors(8) = 1;
scaling_factors(9) = 0.2;
scaling_factors(10) = 0.01;
scaling_factors(11) = 1/200;
scaling_factors(12) = 10^(-3)/V_onecell;
scaling_factors = scaling_factors(~~logic_vec);

l = length(fit_params);
n1 = 9;
V = zeros(l,n1);

for k = 1:l
V(k,:) = linspace(fit_params(k) - 0.5*fit_params(k),fit_params(k) + 0.5*fit_params(k),n1); % matrix of various param values
end

%

%     str3 = str(m); 
%     str5 = sprintf('%g',m);


%%
LB = [1 0.4 5 0.5];
UB = [10, 1.25, 20, 1.5];
n = 2001;
n1 = 1;
tic
[M1,E1,L1] = gibbs_sampler_uniform_medparams(LB,UB,n,n1,@err_fun);
toc
%%
str4 = ["\lambda_{2}", "\lambda_{3}", "\lambda_{9}", "\lambda_{1}"];
num_params = 4;
meds = zeros(1,4);
for k = 1:4
    meds(k) = scaling_factors(k).*median(M1(:,k));
end
for i = 1:4
    str1 = sprintf('%g',i);
    for j = i:4
        str2 = sprintf('%g',j);
        if j == i
            v = linspace(0,350,10);
            q = meds(i)*ones(1,10);
            figure(10*i+j)
            box on
            hist(scaling_factors(i).*M1(:,i))
            xlabel(str4(i))
            ylabel('counts')
            hold on
            plot(q,v,'k--','linewidth',2)
            file_name = [save_file 'unifhist' str1 '.eps'];
            export_fig(figure(10*i+j),file_name)
        else
            figure(10*i+j)
            scatter(scaling_factors(i).*M1(:,i),scaling_factors(j).*M1(:,j),'rx')
            hold on
            xticks = get(gca,'Xtick');
            yticks = get(gca,'Ytick');
            v1 = linspace(min(xticks),max(xticks),10);
            q1 = meds(i)*ones(1,10);
            v2 = linspace(min(yticks),max(yticks),10);
            q2 = meds(j)*ones(1,10);
            plot(v1,q2,'k--','linewidth',2)
            hold on
            plot(q1,v2,'k--','linewidth',2)
            box on
            xlabel(str4(i))
            ylabel(str4(j))
            file_name = [save_file 'unifscatter' str1 str2 '.eps'];
            export_fig(figure(10*i+j),file_name)
        end
    end
end

%%
iters = 1:1:length(E1);

figure(50)
box on
plot(iters,E1)
xlabel('iteration')
ylabel('Error function (E)')
file_name = [save_file 'uniferr_trace.eps'];
export_fig(figure(50),file_name)

%%
for i = 1:4
    str3 = sprintf('%g',i);
    figure(50+i)
    box on
    plot(iters, scaling_factors(i).*M1(:,i))
    xlabel('iteration')
    ylabel(str4(i))
    file_name = [save_file 'uniftrace' str3 '.eps'];
    export_fig(figure(50+i),file_name)
end

%%
mu1 = [202 78.3 10.6 410];
sigma1 = [60 8 1.1 80];
% for i = 1:4
%     str1 = sprintf('%g',i);
%     V =  linspace(min(scaling_factors(i)*M1(:,i)),max(scaling_factors(i)*M1(:,i)),1001);
%     figure(60+i)
%     histogram(scaling_factors(i)*M1(:,i),'Normalization','probability')
%     hold on
%     plot(V,scaling_factors(i)*gauss_distribution(V,mu1(i),sigma1(i)),'r','linewidth',2)
%     xlabel(str4(i))
%     ylabel('normalized counts')
%     file_name = [save_file 'unifnormhist' str1 '.eps'];
%     export_fig(figure(60+i),file_name)
% end
%%
str1 = sprintf('%g',1);
    V =  linspace(min(scaling_factors(1)*M1(:,1)),max(scaling_factors(1)*M1(:,1)),1001);
    figure(61)
    box on
    histogram(scaling_factors(1)*M1(:,1),'Normalization','probability')
    hold on
    plot(V,30*gauss_distribution(V,mu1(1),sigma1(1)),'r','linewidth',2)
    hold on
    yticks = get(gca,'Ytick');
    v = linspace(0,max(yticks),10);
    q = meds(1)*ones(1,10);
    plot(q,v,'k--','linewidth',2)
    xlabel(str4(1))
    ylabel('probability density')
    file_name = [save_file 'unifnormhist' str1 '.eps'];
    export_fig(figure(61),file_name)
 %%   
    
    str1 = sprintf('%g',2);
    V =  linspace(min(scaling_factors(2)*M1(:,2)),max(scaling_factors(2)*M1(:,2)),1001);
    figure(62)
    box on
    histogram(scaling_factors(2)*M1(:,2),'Normalization','probability')
    hold on
    plot(V,2.5*gauss_distribution(V,mu1(2),sigma1(2)),'r','linewidth',2)
    hold on
    yticks = get(gca,'Ytick');
    v = linspace(0,max(yticks),10);
    q = meds(2)*ones(1,10);
    plot(q,v,'k--','linewidth',2)
    xlabel(str4(2))
    ylabel('probability density')
    file_name = [save_file 'unifnormhist' str1 '.eps'];
    export_fig(figure(62),file_name)
    
    %%
    str1 = sprintf('%g',3);
    V =  linspace(min(scaling_factors(3)*M1(:,3)),max(scaling_factors(3)*M1(:,3)),1001);
    figure(63)
    box on
    histogram(scaling_factors(3)*M1(:,3),'Normalization','probability')
    hold on
    plot(V,0.7*gauss_distribution(V,mu1(3),sigma1(3)),'r','linewidth',2)
    hold on
    yticks = get(gca,'Ytick');
    v = linspace(0,0.3,10);
    q = meds(3)*ones(1,10);
    plot(q,v,'k--','linewidth',2)
    xlabel(str4(3))
    ylabel('probability density')
    file_name = [save_file 'unifnormhist' str1 '.eps'];
    export_fig(figure(63),file_name)
    
    %%
    str1 = sprintf('%g',4);
    V =  linspace(min(scaling_factors(4)*M1(:,4)),max(scaling_factors(4)*M1(:,4)),1001);
    figure(64)
    box on
    histogram(scaling_factors(4)*M1(:,4),'Normalization','probability')
    hold on
    plot(V,50*gauss_distribution(V,mu1(4),sigma1(4)),'r','linewidth',2)
    hold on
    yticks = get(gca,'Ytick');
    v = linspace(0,max(yticks),10);
    q = meds(4)*ones(1,10);
    plot(q,v,'k--','linewidth',2)
    xlabel(str4(4))
    ylabel('probability density')
    file_name = [save_file 'unifnormhist' str1 '.eps'];
    export_fig(figure(64),file_name)
%%
figure(65)
box on
plot(iters,L1)
xlabel('iteration')
ylabel('Marginal Likelihood (L)')

%% Credible intervals
for k = 1:4
M2 = zeros(n-n1,4);
M2 = sort(M1(:,k));
LB1 = M2(51)*scaling_factors(k)
UB1 = M2(1950)*scaling_factors(k)
end

%%
function f = err_fun(x)
global exp_data_descriptor initial_guess logic_vec fixed_params
    fit_params = fixed_params;
    fit_params = x;
    f = 30*total_square_error_12_params_vary3(fit_params,exp_data_descriptor,initial_guess,logic_vec)/100;
end