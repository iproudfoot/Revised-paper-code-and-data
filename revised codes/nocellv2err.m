function f = nocellv2err(coeffs,t_data,conc_data)
k1 = coeffs(1);
k2 = coeffs(2);
Kbar = 781;
sol = ode15s(@nocellv2,[0 24],[150 0],[],k1,k2,Kbar);
model = deval(sol,t_data);
err = model(1,:) - conc_data;
f = sum(err.^2);
end

