clear all
close all

save_file = 'C:\Users\isaac\Documents\neat_figs3\'; % Place to save figures
%% Large tissue
v = 0;
L1 = 0;
L = 10;
C_max = 0.0010;
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
    C_ef1(j) = C_max*(1+ tanh(10*0))*(1 + tanh(10*(2 - 0)));
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
    
    sol = ode15s(@pde_2dirich3per,[0 T1],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v,C_max);
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
%     figure(6)
%     plot(200*x,C_tot(:,1)/60)
%     hold on
%     plot(200*x,C_tot(:,3)/60)
%     hold on
%     plot(200*x,C_tot(:,6)/60)
%     hold on
%     plot(200*x,C_tot(:,11)/60)
%     hold on
%     plot(200*x,C_tot(:,16)/60)
%     hold on
%     plot(200*x,C_tot(:,21)/60)
%     hold on
%     xlabel('z \mu m')
%     xlim([0 1000])
%     ylabel('intracellular concentration \mu M')
%     legend('0 mins','15 mins','30 mins','1 hr','1.5 hrs','2 hrs','location','best')
%     file_name = [save_file 'uptake_large.eps'];
%     export_fig(figure(6),file_name);
    %%
    figure(5)
    plot(200*x,C_tot(:,1761)/120)
    hold on
    plot(200*x,C_tot(:,1761)/480+C_tot(:,20)/60)
    hold on
    plot(200*x,C_tot(:,1761)/480+C_tot(:,241)/60)
    hold on
    plot(200*x,C_tot(:,1761)/480+C_tot(:,481)/60)
    hold on
    plot(200*x,C_tot(:,1761)/480+C_tot(:,1761)/120)
    xlabel('z \mu m')
    xlim([0 1000])
    ylabel('intracellular concentration \mu M')
    legend('1 week','1 wk 2 hrs',' 1wk 24hrs',' 1 wk 2 days','2 weeks','location','best')
    file_name = [save_file 'weeklyperprof.eps'];
    export_fig(figure(5),file_name);
    
    %%
       z = [7,7];
    xnew = x(1:n);
    [X1, T1] = meshgrid(x(1:n+1),t);
    C_itot = C_if+C_is+C_ins;
    C_itot_cell = C_tot(1:n+1,:);
    Cstar_loc = zeros(1,length(t1));
    t1 = 0:dt:168;
    for i = 1:length(t1)
        u = find(C_itot_cell(:,i)/60 > 7);
    if isempty(u) == 1
        Cstar_loc(i) = 0;
    else
        Cstar_loc(i) = n - min(u);
    end
    end
    t2 = 0:dt:182;
    figure(10)
    plot(t1,5*Cstar_loc,'b')
    hold on
    plot(168+t1,5*Cstar_loc,'b')
    ylabel('x^* (\mu m)')
    xlabel('time (h)')
    xlim([0 400])
    file_name = [save_file, 'xstar_big_per2.eps'];
    export_fig(figure(10),file_name)
    
    figure(2)
    contourf(T1,(1000-200*X1)/2,C_itot_cell'/60,z)
    ylabel('x^* (\mu m)')
    xlabel('time (h)')
    ylim([0 1000])
    xlim([0 400])
    file_name = [save_file 'xstar_scheme5.png'];
    export_fig(figure(2),file_name)