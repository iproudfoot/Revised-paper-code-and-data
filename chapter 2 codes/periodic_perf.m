clear all
close all
save_file = 'C:\Users\isaac\Documents\neat_figs2\'; % Place to save figures
%%
x1 = 0:0.2:400;
y = 25*(1+ tanh(100*mod(x1,24*7))).*(1 + tanh(100*(2 - mod(x1,24*7))));

figure(3)
plot(x1,y)
xlabel('time (hours)')
ylabel('taxol concentration at boundary (nM)')
file_name = [save_file 'periodicBC.eps'];
export_fig(figure(3),file_name)

%%
Bmaxc = 6930;
ICN = 0.9*10^6;
V_onecell = 4*10^(-6);
D = 0.024;
B0 = 6300;
T0 = 3*B0;
Cn = 0;
t0 = 0;
v = 0;
lengths = [1000 500 250 100];
   for k = 1:4
          L = lengths(k)/1000;
            L1 = lengths(k)/100;
            N = 101;
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
    for j = N
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

            
V_medium = L*10^3;
t = 0:0.05:2*24*7;
dx = L/100;
dt = 0.05;
x = 0:dx:L;
x1 = 1000*x;
x_loc = find(x1==(lengths(k)-100));
N = 101;
D = 0.048;
C_initial = 100;
y_initial = [C_ef1 C_eb1 P1 C_if1 C_is1 C_ins1 B1 C_n1 T_tot1];
opts = odeset('MaxStep',1);
    sol2 = ode15s(@pde_2,[t0 48*7],y_initial, opts, N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t5 = 0:0.05:2*24*7;
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
    C_tot = C_tot2;
    
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
    %%
%%
y_initial1 = [C_ef1 C_ebnd(:,end)' P1(:,end) C_if1(:,end)' C_is1(:,end)' C_ins1(:,end)' B1(:,end)' C_n1(:,end)' T_tot1(:,end)'];

opts = odeset('MaxStep',1);
    sol2 = ode15s(@pde_2,[t0 24*7],y_initial, opts, N,V_medium,D,dx,Bmaxc,n,ICN,v);
    t5 = 0:0.05:24*7;
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
    C_tot3 = C_is+C_if+C_ins;
    C_fin2 = C_tot2(:,end);
    C_tot = [C_tot2 C_tot3];
    
    data = [C_tot, x1'];
    [i2,j2] = find(C_tot3 > 7);
    times1 = [j1/20; 168; j2/20+168];
    prof1 = [1.1*(lengths(k) - i1*L1); 0; 1.1*(lengths(k) - i2*L1)];
    prof2 = [(lengths(k) - i1*L1); 0; (lengths(k) - i2*L1)];
     figure(9)
    if k ==1
        plot(times1,prof1)
        hold on
    end
    plot(times1,prof2)
    hold on
   end
    %%
    y11 = 130*ones(1,21);
     y12 = 0:5:100;
     figure(8)
     plot(y12,y11,'k--')
     ylabel('effective penetration depth (\mu m)')
     xlabel('time (hours)')
     legend({'2000 \mu m','1000 \mu m','500 \mu m','250\mu m', '100 \mu m','proliferating region'},'location','best')
     file_name = [save_file, 'taxvar.eps']; 
export_fig(figure(8),file_name)
       %%
     y13 = 130*ones(1,71);
     y14 = 0:5:350;
     figure(9)
     plot(y14,y13,'k--')
     hold on
     ylabel('effective penetration depth (\mu m)')
     xlabel('time (hours)')
     legend({'2000 \mu m','1000 \mu m','500 \mu m','250\mu m', '100 \mu m','proliferating region'},'location','best')
     file_name = [save_file, 'taxvarper.eps']; 
 export_fig(figure(9),file_name)