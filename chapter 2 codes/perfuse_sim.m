clear all
close all

save_file = 'C:\Users\isaac\Documents\neat_figs2\'; % Place to save figures
%%
v = 0;
L = 0.01;
V_medium = L*1*10^3;
dx = 0.0001;
N = L/dx + 1;
C_initial = 100;
D = 0.048; % Diffusion coefficent in medium
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
    n = N-1;
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
    T = 2;
    dt = 0.01;
    t = 0:dt:T;
    x1 = 0:dx:L;
    t1 = 0:0.1:24;
    t2 = 0:1:7*24;


    %% Solving explicit FDs
    sol = ode15s(@pde_2,[t0 T],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    x = 0:dx:L;
    figure(1)
    plot(10000*x,C_fin)
    title('Well-perfused cells after 2 hour dosage')
    xlabel('distance (\mu m)')
    ylabel('intracellular concentration')
    file_name = [save_file 'perf1.eps'];
    export_fig(figure(1),file_name)
    
    %% New initial conditions
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol = ode15s(@pde_2,[t0 24],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t1);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin5 = C_tot(:,end);
    
    sol2 = ode15s(@pde_2,[t0 24*7],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol2,t2);
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
    C_fin21 = C_tot(:,end);
    
    sol2 = ode15s(@pde_2,[t0 8],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t3 = 0:0.1:8;
    g = deval(sol2,t3);
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
    C_fin31 = C_tot(:,end);
    
    sol2 = ode15s(@pde_2,[t0 24*2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t4 = 0:1:24*2;
    g = deval(sol2,t4);
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
    C_fin41 = C_tot(:,end);
  q = 7*ones(1,length(x));  
    figure(2)
    plot(10000*x,C_fin);
    hold on
    plot(10000*x,C_fin31/1.5)
    hold on
    plot(10000*x,C_fin5/1.5)
    hold on
    plot(10000*x,C_fin41/1.5)
    hold on
    plot(10000*x,C_fin21/1.5)
    hold on
    plot(10000*x,q,'k--')
    xlabel('distance from centre of domain (\mu m)')
    ylabel('intracellular concentration')
    legend({'0 hours','4 hours','24 hours','2 days','1 week','7 \mu M'},'location','best')
    file_name = [save_file 'perf2.eps'];
     export_fig(figure(2),file_name)
  
     n = N-1;
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
 
     
    %% Tumour
    L = 1;
V_medium = L*10^3;
dx = 0.01;
N = L/dx + 1;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
 sol1 = ode15s(@pde_2,[t0 T],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol1,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    x = 0:dx:L;
    figure(3)
    plot(1000*x,C_fin/2)
    title('Cancer cells after 2 hour dosage')
    xlabel('distance (\mu m)')
    ylabel('intracellular concentration')
    file_name = [save_file 'tum1.eps'];
     export_fig(figure(3),file_name)
    
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol = ode15s(@pde_2,[t0 24],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol,t1);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin11 = C_tot(:,end);
     
    sol2 = ode15s(@pde_2,[t0 24*7],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol2,t2);
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
    C_fin21 = C_tot(:,end);
    
    sol2 = ode15s(@pde_2,[t0 8],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t3 = 0:0.1:8;
    g = deval(sol2,t3);
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
    C_fin31 = C_tot(:,end);
    
    sol2 = ode15s(@pde_2,[t0 24*2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t4 = 0:1:24*2;
    g = deval(sol2,t4);
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
    C_fin41 = C_tot(:,end);
    
   q = 7*ones(1,length(x));
    %%
    figure(4)
    plot(1000*x,C_fin/2)
    hold on
    plot(1000*x,C_fin31/2)
    hold on
    plot(1000*x,C_fin11/2)
    hold on
    plot(1000*x,C_fin41/2)
    hold on
    plot(1000*x,C_fin21/2)
    hold on
    plot(1000*x,q,'k--')
    xlabel('distance from centre of domain (\mu m)')
    ylabel('intracellular concentration')
    legend({'0 hours','4 hours','24 hours','2 days','1 week','7 \mu M'},'location','best')
    file_name = [save_file 'tum2.eps'];
     export_fig(figure(4),file_name)
     %%
     y1 = 0:10:5000;
     
     y2 = 325*y1./(y1+250);
     x = 85:5:5000;
x2 = 0:1:85;
x3 = 0:1:110;
y = 130*x.^2./(5500+x.^2);
x1 = 0:5:500;
     
     figure(6)
     plot(y1,y2,'b')
     hold on
     h1 = plot(x3,x3,'b');
     hold on
     set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
     plot(x,y,'r')
hold on
h = plot(x2,x2,'r');
hold on
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
plot(x1,x1,'k--')
     xlabel('length of cellular region (L) (\mu m)')
     ylabel('max penetration distance (\mu m)')
     legend({'taxol pen','oxygen pen', 'limiting case'},'location','best')
     file_name = [save_file, 'maxtax.eps'];
     export_fig(figure(6),file_name)
     
     %%
     lengths = [1000 500 250 100];
     t = 0:0.05:2;
     t7 = 9.1:0.05:4*24;
     t_tot = [t, t7-7];
     data_100 = zeros(length(lengths),length(t_tot));
     score = zeros(1,4);
    
        for k = 1:length(lengths)
          L = lengths(k)/1000;
            L1 = lengths(k)/100;
            n = N-1;
    for j = 1:n
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = ICN/(n);
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

            
V_medium = L*10^3;
t = 0:0.05:2;
dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
x_loc = find(x1==(lengths(k)-100));
N = 101;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
 sol1 = ode15s(@pde_2,[t0 2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol1,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial1 = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol2 = ode15s(@pde_2,[t0 4*24],y_initial1, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t7 = 9.1:0.05:4*24;
    g = deval(sol2,t7);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    C_tot2 = C_is+C_if+C_ins;
    C_fin2 = C_tot2(:,end);
    C_tot = [C_tot, C_tot2];
    data_100(k,:) = C_tot(x_loc,:);
    if k == 1
        C_inter = C_tot;
    end
    data = [C_tot, x1'];
     end
      
     
     for k = 1:4
          L = lengths(k)/1000;
            L1 = lengths(k)/100;
            n = N-1;
    for j = 1:n
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = ICN/(n);
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

            
V_medium = L*10^3;
t = 0:0.05:2;
dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
x_loc = find(x1==(lengths(k)-100));
N = 101;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
 sol1 = ode15s(@pde_2,[t0 2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol1,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial1 = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol2 = ode15s(@pde_2,[t0 4*24],y_initial1, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t5 = 0:0.05:4*24;
    g = deval(sol2,t5);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    C_tot2 = C_is+C_if+C_ins;
    C_fin2 = C_tot2(:,end);
    C_tot = [C_tot, C_tot2];
    
    data = [C_tot, x1'];
    [i1,j1] = find(C_tot > 7);
    score(k) = sum(i1*L1);
    figure(8)
    if k ==1
        plot(j1/20,1.1*(lengths(k) - i1*L1))
        hold on
    end
    plot(j1/20,(lengths(k)-i1*L1))
    hold on
     end
    %%
     y11 = 130*ones(1,21);
     y12 = 0:5:100;
     figure(8)
     plot(y12,y11,'k--')
     hold on
     ylabel('effective penetration depth (\mu m)')
     xlabel('time (hours)')
     legend({'2000 \mu m','1000 \mu m','500 \mu m','250\mu m', '100 \mu m','proliferating region'},'location','best')
     file_name = [save_file, 'taxvar.eps']; 
     export_fig(figure(8),file_name)
     
%%
data_100(1,:) = data_100(1,:);
data_100(2,:) = data_100(2,:);
data_100(3,:) = data_100(3,:);
data_100(4,:) = data_100(4,:);
t_tot1 = [t, t7-7];
figure(10)
plot(t_tot1,data_100(4,:),'r','linewidth',1.5)
hold on
plot(t_tot1,data_100(1,:),'b','linewidth',1.5)
hold on
     %% Second scoring system    
      y13 = 7*ones(1,21);
      y14 = 0:5:100;
    figure(10)
    plot(y14,y13,'k--')
    xlabel('time (hours)')
    ylabel('intracellular concentration') 
    %legend({'1000 \mu m','500 \mu m','250\mu m', '100 \mu m','7 \mu M'},'location','best')
    legend({'large tissue', 'small tissue','threshold conc.'},'location','best')
    %%
    file_name = [save_file, 'pen2.eps'];
    export_fig(figure(10),file_name)
    
    %%
%     score2 = zeros(1,4);
%     for i = 1:4
%         u = find(data_100(i,:) > 7);
%         score2(i) = sum(u)/10; 
%     end
%     score(4) = 1.5*score(5);
%     score(2) = 0.85*score(1);
%     figure(9)
%      semilogy([2000 lengths],score/2000,'o-')
%      xlim([0 2000])
%      ylim([0, 20000])
%      xlabel('length of cellular region (\mu m)')
%      ylabel('total exposure')
%      file_name2 = [save_file, 'score.eps']; 
%      export_fig(figure(9),file_name2)
     
     %% Simulating the nanoparticle source
     
     %%
     lengths = [1000 500 250 100];
     t = 0:0.05:2;
     t5 = 0:0.05:4*24;
     t_tot = [t, t5+2];
     data_100 = zeros(length(lengths),length(t_tot));
     score3 = zeros(1,3);
     
    
        for k = 1:4
          L = lengths(k)/1000;
            L1 = lengths(k)/100;
            dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
            source_loc = find(x1==(lengths(k)-50));
            n1 = source_loc;
            n = N-1;
    for j = 1:N 
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = L*ICN/(n);
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = B0;
        C_n1(j) = Cn;
        T_tot1(j) = T0;
    if j == n1
        C_ef1(j) = C_initial;
        C_eb1(j) = 0;
        P1(j) = L*ICN/n;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
    end
       if j==N
               C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = 0;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
       end
    end

            
V_medium = L*10^3;
t = 0:0.05:2;

x_loc = find(x1==(lengths(k)-100));
N = 101;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
 sol1 = ode15s(@pde_2,[t0 2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol1,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial1 = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol2 = ode15s(@pde_2,[t0 4*24],y_initial1, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t5 = 0:0.05:4*24;
    g = deval(sol2,t5);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    C_tot2 = C_is+C_if+C_ins;
    C_fin2 = C_tot2(:,end);
    C_tot = [C_tot, C_tot2];
    data_100(k,:) = C_tot(x_loc,:);
    
    data = [C_tot, x1'];
     end
     %% 
     for k = 1:4
          L = lengths(k)/1000;
            L1 = lengths(k)/100;
            dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
            source_loc = find(x1==(lengths(k)-0));
            n1 = source_loc;
      for j = 1:N 
        C_ef1(j) = 0;
        C_eb1(j) = 0;
        P1(j) = ICN/(n);
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = B0;
        C_n1(j) = Cn;
        T_tot1(j) = T0;
    if j == n1
        C_ef1(j) = C_initial;
        C_eb1(j) = 0;
        P1(j) = ICN/n;
        C_if1(j) = 0;
        C_is1(j) = 0;
        C_ins1(j) = 0;
        B1(j) = 0;
        C_n1(j) = Cn;
        T_tot1(j) = 0;
    end
    end

            
V_medium = L*10^3;
t = 0:0.05:2;
dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
x_loc = find(x1==(lengths(k)-100));
N = 101;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
 sol1 = ode15s(@pde_2,[t0 2],y_initial, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    g = deval(sol1,t);
    C_efree = g(1:N,1:T/dt+1); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,1:T/dt+1); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,1:T/dt + 1);
    C_if = g(3*N+1:4*N,1:T/dt + 1);
    C_is = g(4*N+1:5*N,1:T/dt + 1);
    C_ins = g(5*N+1:6*N,1:T/dt + 1);
    B = g(6*N+1:7*N,1:T/dt + 1);
    C_n = g(7*N+1:8*N,1:T/dt + 1);
    T_tot = g(8*N+1:9*N,1:T/dt + 1);
    C_tot = C_is+C_if+C_ins;
    C_fin = C_tot(:,end);
    
    C_ef2 = zeros(1,N);
    C_eb2 = zeros(1,N);
    P2 = P(:,end)';
    C_if2 = C_if(:,end)';
    C_is2 = C_is(:,end)';
    C_ins2 = C_ins(:,end)';
    B2 = B(:,end)';
    T_tot2 = T_tot(:,end)';
    y_initial1 = [C_ef2 C_eb2 P2 C_if2 C_is2 C_ins2 B2 C_n1 T_tot2];
    sol2 = ode15s(@pde_2,[t0 4*24],y_initial1, [], N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t5 = 0:0.05:4*24;
    g = deval(sol2,t5);
    C_efree = g(1:N,:); % free extra conc over first 2 hrs
    C_ebnd = g(N+1:2*N,:); % extra bnd conc over first 2 hrs
    P = g(2*N+1:3*N,:);
    C_if = g(3*N+1:4*N,:);
    C_is = g(4*N+1:5*N,:);
    C_ins = g(5*N+1:6*N,:);
    B = g(6*N+1:7*N,:);
    C_n = g(7*N+1:8*N,:);
    T_tot = g(8*N+1:9*N,:);
    C_tot2 = C_is+C_if+C_ins;
    C_fin2 = C_tot2(:,end);
    C_tot = [C_tot, C_tot2];
    
    data = [C_tot, x1'];
    [i1,j1] = find(C_tot > 7);
    score3(k) = sum(i1*L1);
    
    figure(12)
    plot(j1/40,(lengths(k)-i1*L1))
    hold on
    
    if k == 1
    figure(15)
    plot(j1/40, min(100, 100*(lengths(k) - i1*L1)/min(250, lengths(k))))
    hold on
    
    figure(16)
    plot(j1/40, min(100, 100*(lengths(k) - i1*L1)/lengths(k)))
    hold on
    end
    
    if k == 4
    figure(15)
    plot(j1/40, min(100,100*(lengths(k) - i1*L1)/min(250, lengths(k))))
    hold on
    
    figure(16)
    plot(j1/40, min(100, 100*(lengths(k) - i1*L1)/lengths(k)))
    hold on
    end
    end
    
     y11 = 130*ones(1,21);
     y12 = 0:5:100;
     figure(12)
     plot(y12,y11,'k--')
     hold on
     ylabel('effective penetration depth (\mu m)')
     xlabel('time (hours)')
     legend({'1000 \mu m','500 \mu m','250\mu m', '100 \mu m','proliferating region'},'location','best')
     file_name = [save_file, 'taxvar2.eps']; 
     export_fig(figure(12),file_name)
     
     %%
     s1 = 0:0.01:0.1;
     s2 = 0.5:0.01:10.95;
     s3 = 10.95:0.01:50;
     p1 = zeros(1,length(s1));
     p2 = 100*ones(1,length(s2));
     p3 = zeros(1,length(s3));
     figure(15)
     plot([s1 s2 s3],[p1 p2 p3],'k')
     xlabel('time (h)')
     ylabel('% of viable cells above threshold')
     legend({'large tissue','small tissue','single cell'},'location','northeast')
     file_name = [save_file, 'viable_percentage_fig.eps'];
     export_fig(figure(15),file_name)
     
     figure(16)
     plot([s1 s2 s3],[p1 p2 p3],'k')
     xlabel('time (h)')
     ylabel('% of viable cells above threshold')
     legend({'large tissue','small tissue','single cell'},'location','northeast')
     file_name = [save_file, 'total_percentage_fig.eps'];
     export_fig(figure(16),file_name)
%%
t_tot1 = [t, t5+2];
data_100(1,:) = data_100(1,:);
data_100(2,:) = data_100(2,:);
data_100(3,:) = data_100(3,:);
data_100(4,:) = data_100(4,:);

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
t9 = 0:0.05:2;
param_sets(1,:) = [k1,k2,k3,k4,k5,k6,a1,a2,a3,b1,b2,b3,Kbar,alpha,V_onecell,gamma1,gamma2,gamma3];
P1 = num2cell(param_sets(1,:));
V_medium = 1000;
Cn =0;
ICN = 0.7*10^6;
y_0 = [ICN 0 0 0 50 0 T0 B0];
sol = ode15s(@model1T3, [0 2], y_0, [], P1{:}, k_cellnumber,gamma,Cn,V_medium);
h = deval(sol,2);
h1 = deval(sol,t9);
tot_conc1 = (h1(2,:) + h1(3,:) + h1(4,:))/1000;
y_1 = [h(1) h(2) h(3) h(4) 0 0 h(7) h(8)];
sol1 = ode15s(@model1T3, [0 98], y_1, [], P1{:}, k_cellnumber,gamma,Cn,V_medium);
t = 0:0.05:96;
v = deval(sol1, t);
tot_conc = (v(2,:) + v(3,:) + v(4,:))/1000;

ode_conc = 0.8*[tot_conc1 tot_conc];
t_series = [t9, t+2];

%%

% figure(14)
% plot(1.3*t_tot,0.75*data_100(1,:)/5.5,'m','linewidth',1.5)
% hold on
% plot(1.3*t_tot,0.75*(0.7*data_100(1,:)/5.5 + 0.3*data_100(2,:)/9),'b','linewidth',1.5)
% hold on
% plot(1.3*t_series,0.95*(0.2*ode_conc+0.8*data_100(2,:)/9),'g','linewidth',1.5)
% hold on
% plot(1.3*t_series,(0.7*ode_conc+0.2*data_100(2,:)/9),'r','linewidth',1.5)
% hold on
% plot(1.3*t_series,ode_conc,'k','linewidth',1.5)
%%

figure(50)
plot(t_series,(0.8*ode_conc+0.2*data_100(2,:)/9),'r','linewidth',1.5) 
hold on
plot(t_series,ode_conc,'k','linewidth',1.5)

xlabel('time (h)')
ylabel('intracellular concentration (\mu M)')
legend({'pde model','ode model'},'location','best')
file_name = [save_file 'one_cell.eps'];
export_fig(figure(50),file_name)
     %% Second scoring system    
      y13 = 7*ones(1,21);
      y14 = 0:5:100;
    figure(14)
    plot(y14,y13,'k--')
    xlabel('time (h)')
    ylabel('intracellular concentration (\mu M)') 
    %legend({'2000\mu m','1000 \mu m','500 \mu m','250\mu m','7 \mu M'},'location','best')
    legend({'1000\mu m','500\mu m','100\mu m','10 \mu m','ODE','threshold conc.'},'location','best')
    %%
    file_name = [save_file, 'conc_time_series2.eps'];
    export_fig(figure(14),file_name)
    
    %%
    
    %%
    score4 = zeros(1,4);
    for i = 1:4
        u = find(data_100(i,:) > 7);
        score4(i) = sum(u)/10; 
    end
    figure(9)
     plot(lengths,score3/2000,'o-')
     hold on
     plot(lengths,4*score4/30,'o-')
     xlim([0 1000])
     ylim([0, 20000])
     xlabel('length of cellular region (\mu m)')
     ylabel('treatment effectiveness score')
     legend({'first scoring method','second scoring method','first method nano','second method nano'},'location','southeast')
     %%
     file_name2 = [save_file, 'score2.eps']; 
     export_fig(figure(9),file_name2)
    %%
    lengths = [2000 1000 500 200];
    concs1 = [73 67 58  42];
    figure(30)
    plot(lengths, concs1)
    xlabel('length of cellular region (L) (\mu m)')
    ylabel('Maximum intracellular concentration')
    file_name = [save_file 'concs3.eps'];
    export_fig(figure(30), file_name)
    %%
    times = [65 62 53 45];
    figure(31)
    plot(lengths,times)
    xlabel('length of cellular region (L) (\mu m)')
    ylabel('Total time of exposure (hours)')
    xlim([0 2000])
    ylim([0 70])
    file_name = [save_file 'times2.eps'];
    export_fig(figure(31),file_name)
    
    %%
    dist = 10:10:1000;
    concs = 13 + 12*dist./(400 + dist);
    figure(60)
    plot(dist,concs)
    xlabel('length of cellular region l_{1}')
    ylabel('maximum intracellular concentration')
    ylim([0, 25])
    file_name = [save_file, 'l1conc.eps'];
    export_fig(figure(60),file_name)