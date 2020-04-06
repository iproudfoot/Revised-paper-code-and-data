function f = pde_nondim(t,y,N,V_medium,D,dx,Bmaxc,n,ICN,v)
v =0;
Bmaxm = 3.94*10^3;
Bmaxc = 6300;
Kbar=781;
K=4.93;
lambda_1 = 2*10^(-3);
L = 1;

rho = 0.0016;
V_onecell = 2.25*10^(-6);
NSB= 10;
lambda_1 = 80000;
lambda_2 = 16;
lambda_3 = 70;
lambda_4 = 0.076;
lambda_5 = 0.1;
K = 781;
B0 = Bmaxc;
T0 = 100*B0/33; 

K = 0;
B0 = Bmaxc;
T0 = 100*B0/33; 
% Estimated parameters

lambda_7 = 0.02;
lambda_6 = lambda_7*(T0-B0);
lambda_8 = 21.6;
lambda_9 = 0.0023;
b3 = 5000;
b4 = 0;
b5 = 0;
d = 0;
gamma = 1000;
k_fast = 10^6;

%
C_ef = y(1:N);
C_eb = y(N+1:2*N);
P = y(2*N+1:3*N);
C_if = y(3*N+1:4*N);
C_is = y(4*N+1:5*N);
C_ins = y(5*N+1:6*N);
B = y(6*N+1:7*N);
C_n = y(7*N+1:8*N);
T_tot = y(8*N+1:9*N);

% PDE discretisation and RHS for ODE solver for t


for i = 1:n
    if i == 1 % BC at x = 0
        Bt(i) = (lambda_8 + lambda_9*C_is(i))*(T_tot(i)-B(i)) - 2*lambda_8*B(i) - gamma*C_n(i)*B(i);
        C_eft(i) = D/dx^2*(C_ef(i+1)-C_ef(i)) + lambda_1*rho/(1 - rho)*(C_if(i) - C_ef(i))*P(i) - lambda_4*C_ef(i) + lambda_5*C_eb(i);
        C_ebt(i) = lambda_4*C_ef(i) - lambda_5*C_eb(i);
        Pt(i) = 0;
        C_ift(i) =  lambda_1*(C_ef(i) - C_if(i))/V_onecell - k_fast*lambda_1*(B(i)-C_is(i))*C_if(i) + k_fast*lambda_2*C_is(i) - k_fast*lambda_3*C_if(i) + k_fast*k4*C_ins(i) - Pt(i)*C_if(i)/P(i);
        C_ist(i) = k_fast*lambda_1*(B(i)-C_is(i))*C_if(i) - k_fast*lambda_2*C_is(i) - Pt(i)*C_is(i)/P(i);
        C_inst(i) = k_fast*lambda_3*C_if(i) - k_fast*C_ins(i) - Pt(i)*C_ins(i)/P(i);
        C_nt(i) = 0;
        T_tot(i) = lambda_6 - lambda_7*(T_tot(i)-B(i));
    else 
        Bt(i) = (lambda_8+lambda_9*C_is(i))*(T_tot(i)-B(i)) - 2*lambda_8*B(i) - gamma*C_n(i)*B(i);
        C_eft(i) = D/dx^2*(C_ef(i-1) - 2*C_ef(i) + C_ef(i+1)) + lambda_1*rho/(1 - rho)*(C_if(i) - C_ef(i))*P(i) - lambda_4*C_ef(i) + lambda_5*C_eb(i);
        C_ebt(i) = D/dx^2*(C_eb(i-1) - 2*C_eb(i) + C_eb(i+1)) + lambda_4*C_ef(i) - lambda_5*C_eb(i);
        Pt(i) = 0;
        C_ift(i) =  lambda_1*(C_ef(i) - C_if(i))/V_onecell - k_fast*lambda_1*(B(i)-C_is(i))*C_if(i) + k_fast*lambda_2*C_is(i) - k_fast*lambda_3*C_if(i) + k_fast*k4*C_ins(i)- Pt(i)*C_if(i)/P(i);
        C_ist(i) = k_fast*lambda_1*(B(i)-C_is(i))*C_if(i) - k_fast*lambda_2*C_is(i)- Pt(i)*C_is(i)/P(i);
        C_inst(i) = k_fast*lambda_3*C_if(i) - k_fast*k4*C_ins(i)- Pt(i)*C_ins(i)/P(i);
        C_nt(i) = 0;
        T_tot(i) = lambda_6 - lambda_7*(T_tot(i)-B(i));
    end
    for i = (n+1):N
    if i == N
        C_eft(i) = D/dx^2*(C_ef(i-1)-C_ef(i)) - lambda_4*C_ef(i) + lambda_5*C_eb(i);
        C_ebt(i) =  lambda_4*C_ef(i) - lambda_5*C_eb(i);
        Pt(i) = 0;
        C_ift(i) = 0;
        C_ist(i) = 0;
        C_inst(i) = 0;
        Bt(i) = 0;
        C_nt(i) = 0;
        Tt(i) = 0;
    else
    C_eft(i) = D/dx^2*(C_ef(i-1) - 2*C_ef(i) + C_ef(i+1)) - lambda_4*C_ef(i) + lambda_5*C_eb(i);
    C_ebt(i) =  lambda_4*C_ef(i) - lambda_5*C_eb(i);
    Pt(i) = 0;
    C_ift(i) = 0;
    C_ist(i) = 0;
    C_inst(i) = 0;
    Bt(i) = 0;
    C_nt(i) = 0;
    Tt(i) = 0;
    end
end
f = [C_eft C_ebt Pt C_ift C_ist C_inst Bt C_nt Tt];
f = f';
end