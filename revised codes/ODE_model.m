function f = ODE_model(t,y,til_lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda_6,lambda_7,lambda_8,lambda_9,V_onecell,V_medium,tilde_rho)
% y(1) = rho (constant)
% y(2) = C_if
% y(3) = C_is
% y(4) = C_ins
% y(5) = C_ef
% y(6) = C_es
% y(7) = T
% y(8) = B 
k_fast = 10^6; % QSSA
f = zeros(8,1);

f(1) = 0;
f(2) = til_lambda_1*(y(5)-y(2)) - kfast*((y(8)-y(3))*y(2) - lambda_2*y(3)) - k_fast*(lambda_3*y(2) - y(4));
f(3) = k_fast*((y(8)-y(3))*y(2) - lambda_2*y(3));
f(4) = k_fast*(lambda_3*y(2) - y(4));
f(5) = til_lambda_1*tilde_rho/(1 - tilde_rho)*(y(2)-y(5)) - lambda_4*y(5) + lambda_5*y(6);
f(6) = lambda_4*y(5) - lambda_5*y(6);
f(7) = lambda_6 - lambda_7*(y(7)-y(8));
f(8) = (lambda_8 + lambda_9*y(3))*(y(7)-y(8)) - 2*lambda_8*y(8);