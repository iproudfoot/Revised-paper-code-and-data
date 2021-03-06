function f = pde_2(t,y,N,V_medium,D,dx,Bmaxc,n,ICN,v)

% Put in known parameters
Bmaxm = 3.94*10^3;
Bmaxc = 6930;
Kbar=781;
K=4.93;
Clf = 2*10^(-3);
L = 1;
k_fast = 10^6;

if ICN == 0
    Clf = 0;
end
V_onecell = 4*10^(-6);
NSB= 10;
k1 = 4;
k2 = 16;
k3 = 700;
k4 = 12;
k5 = 0.72;
k6 = 0.11;
K = 781;
B0 = Bmaxc;
T0 = 100*B0/33; 
v = 0;


% Estimated parameters
a2 = 0;
a3 = 12;
a4 = 0.02;
a1 = a4*(T0-B0);
b1 = 21.6;
b2 = 0.023;
b3 = 5000;
b4 = 0;
b5 = 0;
d = 0;
gamma = 1000;

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
        Bt(i) = (b1 + b2*C_is(i))*(T_tot(i)-B(i)) - 2*b1*B(i) - gamma*C_n(i)*B(i);
        C_eft(i) = D/dx^2*(C_ef(i+1)-C_ef(i))+ v/dx*(C_ef(i+1)-C_ef(i)) + Clf*(C_if(i) - C_ef(i))*P(i) - k5*C_ef(i) + k6*C_eb(i);
        C_ebt(i) = D/dx^2*(C_eb(i+1)-C_eb(i))+ v/dx*(C_eb(i+1)-C_eb(i)) + k5*C_ef(i) - k6*C_eb(i);
        Pt(i) = 0;
        C_ift(i) =  Clf*(C_ef(i) - C_if(i))/V_onecell - k1*(B(i)-C_is(i))*C_if(i) + k2*C_is(i) - k3*C_if(i) + k4*C_ins(i) - Pt(i)*C_if(i)/P(i);
        C_ist(i) = k_fast*(k1*(B(i)-C_is(i))*C_if(i) - k2*C_is(i));
        C_inst(i) = k_fast*(k3*C_if(i) - k4*C_ins(i));
        C_nt(i) = 0;
        T_tott(i) = a1 - a4*(T_tot(i)-B(i));
    else 
        Bt(i) = (b1+b2*C_is(i))*(T_tot(i)-B(i)) - 2*b1*B(i) - gamma*C_n(i)*B(i);
        C_eft(i) = D/dx^2*(C_ef(i-1) - 2*C_ef(i) + C_ef(i+1)) + Clf*(C_if(i) - C_ef(i))*P(i) - k5*C_ef(i) + k6*C_eb(i);
        C_ebt(i) = D/dx^2*(C_eb(i-1) - 2*C_eb(i) + C_eb(i+1)) + k5*C_ef(i) - k6*C_eb(i);
        Pt(i) = 0;
        C_ift(i) =  Clf*(C_ef(i) - C_if(i))/V_onecell - k_fast*k1*(B(i)-C_is(i))*C_if(i) + k_fast*k2*C_is(i) - k_fast*k3*C_if(i) + k_fast*k4*C_ins(i)- Pt(i)*C_if(i)/P(i);
        C_ist(i) = k_fast*(k1*(B(i)-C_is(i))*C_if(i) - k2*C_is(i)) - Pt(i)*C_is(i)/P(i);
        C_inst(i) = k_fast*(k3*C_if(i) - k4*C_ins(i)) - Pt(i)*C_ins(i)/P(i);
        C_nt(i) = 0;
        T_tott(i) = a1 - a4*(T_tot(i)-B(i));
    end
end
for i = (n+1):N
    if i == N
        C_eft(i) = D/dx^2*(C_ef(i-1)-C_ef(i)) - k5*C_ef(i) + k6*C_eb(i);
        C_ebt(i) =  k5*C_ef(i) - k6*C_eb(i);
        Pt(i) = 0;
        C_ift(i) = 0;
        C_ist(i) = 0;
        C_inst(i) = 0;
        Bt(i) = 0;
        C_nt(i) = 0;
        Tt(i) = 0;
    else
    C_eft(i) = D/dx^2*(C_ef(i-1) - 2*C_ef(i) + C_ef(i+1)) - k5*C_ef(i) + k6*C_eb(i);
    C_ebt(i) =  k5*C_ef(i) - k6*C_eb(i);
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
