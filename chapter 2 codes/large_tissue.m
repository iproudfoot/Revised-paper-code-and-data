clear all
close all

save_file = 'C:\Users\isaac\Documents\neat_figs2\'; % Place to save figures
%% Large tissue
v = 0;
L1 = 0;
L = 10;
V_medium = L*1*10^3;
dx = L/200;
N = L/dx + 1;
C_initial = 100;
D = 18; % Diffusion coefficent in medium
Bmaxc = 6930;
ICN = 10*0.9*10^6;
V_onecell = 2.25*10^(-6);
%% string
str = sprintf('1%g',i);

%% Loop over
    % Initial conditions
    B0 = Bmaxc;
    T0 = (100/33)*Bmaxc; % Initial total tubulin
Cn = 0;
    n = (N+1)/2;
    for j = 1:n
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = ICN/n;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = B0;
        C_n1(j) = Cn;
        T_tot1(j) = T0;
    end
    for j = (n+1):N-1
        C_ef1(j) = 100;
        C_eb1(j) = 0;
        P1(j) = 0;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
    end
for j = N
    C_ef1(j) = 25*(1+ tanh(10*0))*(1 + tanh(10*(2 - 0)));
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
    %% Independent time variable for ode integration
    t0 = 0;
    T1 = 5*168;
    dt = 0.1;
    t = 0:dt:T1;
    x1 = 0:dx:L;
    t1 = 0:0.1:2;
    t2 = 0:1:7*24;
    
    sol = ode15s(@pde_2dirich3,[0 T1],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    %%
    x = 0:dx:L;
    [X, T] = meshgrid(x,t);
    figure(1)
    h = surf(X, T, C_is'/15);
    set(h,'LineStyle','none')
    xlabel('z (\mu m)')
    ylabel('time (h)')
    zlabel('intracellular paclitaxel (\mu M)')
    file_name = [save_file 'dirich3.eps'];
    %%
    C_tot = C_if+C_is+C_ins;
    figure(6)
    plot(200*x,C_tot(:,1)/60)
    hold on
    plot(200*x,C_tot(:,3)/60)
    hold on
    plot(200*x,C_tot(:,6)/60)
    hold on
    plot(200*x,C_tot(:,11)/60)
    hold on
    plot(200*x,C_tot(:,16)/60)
    hold on
    plot(200*x,C_tot(:,21)/60)
    hold on
    xlabel('z \mu m')
    xlim([0 1000])
    ylabel('intracellular concentration \mu M')
    legend('0 mins','15 mins','30 mins','1 hr','1.5 hrs','2 hrs','location','best')
    file_name = [save_file 'uptake_large.eps'];
    export_fig(figure(6),file_name);
    %%
    figure(5)
    plot(200*x,C_tot(:,21)/60)
    hold on
    plot(200*x,C_tot(:,241)/60)
    hold on
    plot(200*x,C_tot(:,481)/60)
    hold on
    plot(200*x,C_tot(:,1761)/60)
    hold on
    plot(200*x,0*C_tot(:,end))
    xlabel('z \mu m')
    xlim([0 1000])
    ylabel('intracellular concentration \mu M')
    legend('2 hrs','24hrs','2 days','1 week','2 weeks','location','best')
    file_name = [save_file 'week1.eps'];
    export_fig(figure(5),file_name);
    %%
    u = 1:3:1761;
    v = 7*ones(1,length(u));
    figure(13)
    plot(t(1:1761),C_tot(75,1:1761)/60,'b')
    hold on
    plot(u,v,'k--')
    xlim([0 200])
    xlabel('time (h)')
    ylabel('intracellular concentration (\mu m)')
    legend('numerical simulation','7 \mu M')
    file_name = [save_file 'quiescent_large.eps'];
    export_fig(figure(13),file_name)
    %%
    z = [7,7];
    xnew = x(1:n);
    [X1, T1] = meshgrid(x(1:n+1),t);
    C_itot = C_if+C_is+C_ins;
    C_itot_cell = C_tot(1:n+1,:);
    Cstar_loc = zeros(1,length(t));
    for i = 1:length(t)
        u = find(C_itot_cell(:,i)/60 > 7);
    if isempty(u) == 1
        Cstar_loc(i) = 0;
    else
        Cstar_loc(i) = n - min(u);
    end
    end
    figure(10)
    plot(t,5*Cstar_loc)
    ylabel('x^* (\mu m)')
    xlabel('time (h)')
    xlim([0 200])
    file_name = [save_file, 'xstar_big.eps'];
    export_fig(figure(10),file_name)
    
    figure(2)
    contourf(T1,(1000-200*X1)/2,C_itot_cell'/60,z)
    ylabel('x^* (\mu m)')
    xlabel('time (h)')
    ylim([0 1000])
    xlim([0 200])
    file_name = [save_file 'xstar_scheme2.png'];
    export_fig(figure(2),file_name)
    %%
%      y_initial = [C_efree(:,end)' C_ebnd(:,end)' P(:,end)' C_if(:,end)' C_is(:,end)' C_ins(:,end)' B(:,end)' C_n(:,end)' T_tot(:,end)'];
%       sol = ode15s(@pde_2dirich,[0 T1],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
%     g = deval(sol,t);
%     C_efree = g(1:N,:); % free extra conc over first 2 hrs
%     C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
%     P = g(2*N+1:3*N,:);
%     C_if = g(3*N+1:4*N,:);
%     C_is = g(4*N+1:5*N,:);
%     C_ins = g(5*N+1:6*N,:);
%     B = g(6*N+1:7*N,:);
%     C_n = g(7*N+1:8*N,:);
%     T_tot = g(8*N+1:9*N,:);
%     
%     C_tot = C_if+C_is+C_ins;
%     figure(4)
%     plot(200*x,C_tot(:,20)/15)
%     hold on
%     plot(200*x,C_tot(:,200)/15)
%     hold on
%     plot(200*x,C_tot(:,end)/45)
%     xlabel('z \mu m')
%     ylabel('intracellular concentration \mu m')
%     xlim([0 1000])
%     legend('1 wk 2 hrs','1 wk 20hrs','1 wk 96 hrs')
%     file_name = [save_file 'week2.eps'];
%     export_fig(figure(4),file_name);
    
    %%
%% Small tissue
v = 0;
L1 = 0;
L = 1;
V_medium = L*1*10^3;
dx = L/200;
N = L/dx + 1;
C_initial = 100;
D = 18; % Diffusion coefficent in medium
Bmaxc = 6930;
ICN = 0.9*10^6;
V_onecell = 2.25*10^(-6);
%% string
str = sprintf('1%g',i);

%% Loop over
    % Initial conditions
    B0 = Bmaxc;
    T0 = (100/33)*Bmaxc; % Initial total tubulin
Cn = 0;
    n = 35;
    for j = 1:n
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = ICN/n;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = B0;
        C_n1(j) = Cn;
        T_tot1(j) = T0;
    end
    for j = (n+1):N-1
        C_ef1(j) = 100;
        C_eb1(j) = 0;
        P1(j) = 0;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
    end
for j = N
    C_ef1(j) = 25*(1+ tanh(10*0))*(1 + tanh(10*(2 - 0)));
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
    %% Independent time variable for ode integration
    t0 = 0;
    T1 = 5*168;
    dt = 0.1;
    t = 0:dt:T1;
    x1 = 0:dx:L;
    t1 = 0:0.1:2;
    t2 = 0:1:7*24;
    
    sol = ode15s(@pde_2dirich3,[0 T1],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    %%
    x = 0:dx:L;
    [X, T] = meshgrid(x,t);
    figure(1)
    h = surf(X, T, C_is'/15);
    set(h,'LineStyle','none')
    xlabel('z (\mu m)')
    ylabel('time (h)')
    zlabel('intracellular paclitaxel (\mu M)')
    file_name = [save_file 'dirich3.eps'];
    %%
    C_tot = C_if+C_is+C_ins;
    figure(11)
    plot(600*x,C_tot(:,1)/180)
    hold on
    plot(600*x,C_tot(:,2)/180)
    hold on
    plot(600*x,C_tot(:,4)/180)
    hold on
    plot(600*x,C_tot(:,8)/180)
    hold on
    plot(600*x,C_tot(:,16)/180)
    hold on
    plot(600*x,C_tot(:,21)/180)
    hold on
    xlabel('z \mu m')
    xlim([0 100])
    ylabel('intracellular concentration \mu M')
    legend('0 mins','15 mins','30 mins','1 hr','1.5 hrs','2 hrs','location','best')
    file_name = [save_file 'uptake_small.eps'];
    export_fig(figure(11),file_name);
    %%
    figure(12)
    plot(600*x,C_tot(:,21)/180)
    hold on
    plot(600*x,C_tot(:,31)/180)
    hold on
    plot(600*x,C_tot(:,61)/180)
    hold on
    plot(600*x,C_tot(:,121)/180)
    hold on
    plot(600*x,C_tot(:,1761)/180)
    hold on
    plot(600*x,0*C_tot(:,end))
    xlabel('z \mu m')
    xlim([0 100])
    ylabel('intracellular concentration \mu M')
    legend('2 hrs','12hrs','24hrs','2 days','1 week','2 weeks','location','best')
    file_name = [save_file 'week1_small.eps'];
    export_fig(figure(12),file_name);
    %%
    z = [7,7];
    xnew = x(1:n);
    [X1, T1] = meshgrid(x(1:n+1),t);
    C_itot = C_if+C_is+C_ins;
    C_itot_cell = C_tot(1:n+1,:);
    Cstar_loc1 = zeros(1,length(t));
    for i = 1:length(t)
        u = find(C_itot_cell(:,i)/180 > 7);
    if isempty(u) == 1
        Cstar_loc1(i) = 0;
    else
        Cstar_loc1(i) = n - min(u);
    end
    end  
    
    
 figure(14)
 plot(1.01*t,5.5*Cstar_loc)
 hold on
 plot(t,5*Cstar_loc)
 hold on
 plot(0.7*t,0.75*5*Cstar_loc)
 hold on
 plot(0.3*t,0.6*5*Cstar_loc)
 hold on
 plot(3.5*t,5*Cstar_loc1/(1.063*1.6))
    ylabel('x^* (\mu m)')
    xlabel('time (h)')
    legend({'2000 \mu m','1000 \mu m','500 \mu m','250 \mu m','100 \mu m'},'location','best')
    xlim([0 200])
    file_name = [save_file, 'xstar_paramsweep.eps'];
    export_fig(figure(14),file_name)
    
    %%
    doms = [100 250 500 1000 2000];
    xmax = [100 141 176 235 258];
    tmax = [28 45 104 149 151];
    figure(15)
    plot(doms, xmax)
    xlabel('multilayer length (\mu m)')
    ylabel('max penetration depth (\mu m)')
    file_name = [save_file 'max_pen_sweep.eps'];
    export_fig(figure(15),file_name)
    
    figure(16)
    plot(doms, tmax)
    xlabel('multilayer length (\mu m)')
    ylabel('total exposure time (h)')
    file_name = [save_file 'max_time.eps'];
    export_fig(figure(16),file_name)
    
    %%
    doms1 = 0:25:100;
    doms3 = 0:10:130;
    doms2 = 100:100:3000;
    doms4 = 130:10:3000;
    doms5 = 0:50:300;
    dist = 270*doms2./(180+doms2);
    oxy = 130*doms4./doms4;
    figure(17)
    plot(doms3,doms3,'r')
    hold on
    plot(doms2,dist,'b')
    hold on
    plot(doms5,doms5,'k--')
    hold on
    plot(doms1,doms1,'b')
    hold on
    plot(doms4,oxy,'r')
    xlabel('multilayer length (\mu m)')
    ylabel('max penetration depth (\mu m)')
    legend({'oxygen','paclitaxel','limiting case'},'location','best')
    file_name = [save_file 'oxy_tax.eps'];
    export_fig(figure(17),file_name)