clear all
close all

save_file = 'C:\Users\isaac\Documents\neat_figs2\'; % Place to save figures
%%
DX = [0.02 0.01 0.005 0.002];
for k = 1:4
str1 = sprintf('%g',k);
v = 0;
L = 1;
V_medium = L*1*10^3;
dx = DX(k);
N = L/dx + 1;
C_initial = 100;
D = 0.48; % Diffusion coefficent in medium
Bmaxc = 6930;
ICN = 0.9*10^6;
V_onecell = 4*10^(-6);
%% string
str = sprintf('1%g',i);

%% Loop over with and without noc
    % Initial conditions
    B0 = Bmaxc;
    T0 = (100/33)*Bmaxc; % Initial total tubulin
Cn = 0;
    n = 0.1/dx;
    for j = 1:n
        C_ef1(j) = C_initial;
        C_eb1(j) = 0;
        P1(j) = ICN/n;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = B0;
        C_n1(j) = Cn;
        T_tot1(j) = T0;
    end
    for j = (n+1):N
        C_ef1(j) = C_initial;
        C_eb1(j) = 0;
        P1(j) = 0;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
    end

    y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
    Amt_e_int = trapz(0:dx:L,C_ef1);
    %% Independent time variable for ode integration
    t0 = 0;
    T = 30;
    dt = 0.01;
    
    x = 0:dx:L;
    t1 = 0:0.1:T;
    t2 = 0:1:7*24;


    %% Solving explicit FDs
    sol = ode15s(@pde_23,[t0 T],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t1);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    C_ef2 = zeros(4,N);
    g2 = deval(sol,1);
    C_ef2(k,:) = g2(1:N); % free extra conc over first 2 hrs
    C_ebnd2 = g2(N+1:2*N); % extra bnd conc over first 2 hrs
    P2 = g2(2*N+1:3*N);
    C_if2 = g2(3*N+1:4*N);
    C_is2 = g2(4*N+1:5*N);
    C_ins2 = g2(5*N+1:6*N);
    B2 = g2(6*N+1:7*N);
    C_n2 = g2(7*N+1:8*N);
    T_tot2 = g2(8*N+1:9*N);
    C_tot2 = zeros(4,N);
    C_tot2(k,:) = C_is2+C_if2+C_ins2;
    C_fin2 = C_tot(:,end);
    
    c = size(C_efree);
    C_efree_avg = zeros(1,c(1));
    C_ebnd_avg = zeros(1,c(1));
%%    

    figure(10)
    plot(t2,C_efree(21,1:169),'b')
    hold on
    plot(t2,C_ebnd(21,1:169),'r')
    xlim([0 25])
    xlabel('incubation time (h)')
    ylabel('extracellular concnetration (nM)')
    file_name = [save_file 'med_avg.eps'];
    export_fig(figure(10),file_name)
  %%  
    [X, T] = meshgrid(0.1*x,t1);
    
    figure(1)
    surf(X,T,C_efree')
    shading interp
    ylim([0 4])
    xlabel('z (mm)')
    ylabel('t (h)')
    zlabel('C_{ef} (nM)')
    file_name = [save_file 'Cefdynam' str1 '.eps'];
    %export_fig(figure(1),file_name)
%%    
    figure(2)
    surf(X,T,5.4*C_is'/1000)
    shading interp
 xlabel('z (mm)')
    ylabel('t (h)')
    zlabel('C_{is} (\mu M)')
    xlim([0, 0.01])
    file_name = [save_file 'Cisdynam.eps'];
    %export_fig(figure(2),file_name)
  %%  
    figure(3)
    surf(X,T,C_if')
    shading interp
    xlabel('z (mm)')
    ylabel('t (h)')
    zlabel('C_{if} (nM)')
    xlim([0, 0.01])
    file_name = [save_file 'Cifdynam.eps'];
    %export_fig(figure(3),file_name)
    %%
    figure(4)
    surf(X,T,C_ins')
    shading interp
    xlabel('z (mm)')
    ylabel('t (h)')
    zlabel('C_{is} (\mu M)')
    xlim([0, 0.01])
    file_name = [save_file 'Cinsdynam.eps'];
    %export_fig(figure(4),file_name)
    %%
    figure(5)
    surf(X,T,C_ebnd')
    shading interp
    xlabel('z (mm)')
    ylabel('t (h)')
    zlabel('C_{ens} (nM)')
    file_name = [save_file 'Censdynam' str(1) '.eps'];
    %export_fig(figure(5),file_name)
    %%
    figure(20)
    plot(x,C_ef2(k,:))
    hold on
    xlabel('z (mm)')
    ylabel('C_{ef}(t) (nM)')
    %%
    figure(21)
    plot(x,C_tot2(k,:))
    hold on
    xlabel('z (mm)')
    ylabel('C(t) (nM)')
end