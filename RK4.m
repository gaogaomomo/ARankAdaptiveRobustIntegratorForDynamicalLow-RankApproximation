function Y1 = RK4(h, Z0, f)
% in our RK4, the function is independent of time
k1 = f(Z0); % here Z0 is our matrix that we want to evaluate in each step
k2 = f(Z0 + h*k1/2);
k3 = f(Z0 + h*k2/2);
k4 = f(Z0 + h*k3);
Y1 = Z0 + h/6 * (k1 + 2*k2 + 2*k3 + k4);