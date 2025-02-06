function [err_norm, err_energy] = err(Yn, Y0, D, V_cos)
err_norm = abs(norm(Yn, 'fro') - norm(Y0, 'fro'));
H0 = disc_schroedinger(D, V_cos, Y0);
Hn = disc_schroedinger(D, V_cos, Yn);
% inner product of matrices is the Frobenius inner product
err_energy = abs(trace(Yn' * Hn) - trace(Y0' * H0));
% absolute errors between our solution and the approximation with 'ode45'
% options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
% An = ode45(@(t, Yn)Hn, [0,1], Y0, options);
% abs_err = norm(Yn - An, 'fro');

function H = disc_schroedinger(D, V_cos, Y)
H = 1/2 * (D * Y + Y * D') + V_cos * Y * V_cos';