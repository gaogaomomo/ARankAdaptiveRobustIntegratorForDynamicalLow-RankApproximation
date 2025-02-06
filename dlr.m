function [U1, S1, V1] = dlr(t0, h, U0, S0, V0, r0, theta, n, rmax)
Y0 = U0 * S0 * V0'; % record the original Y0 for error estimation
for i = 1 : n
    r(i) = r0;
    % [U0, S0, V0] is the SVD decomposition for Y0, S0 ist r×r matrix
    % K-step
    K0 = U0 * S0; % m×r matrix U0, K0
    % compute K(t1) with classic Runge–Kutta method
    % evaluate K1 = K0 + h * (F(t0, K0 * V0') * V0) with input function;
    Fk_tilde = @(K0) F_k(K0, V0);
    K1 = RK4(h, K0, Fk_tilde);
    [U_hat, ~] = qr([K1, U0], 0); % U_hat is the orth. basis of m×2r matrix
    M = U_hat' * U0; % 2r×r matrix M

    % L-step
    L0 = V0 * S0'; % n×r matrix V0, L0
    % compute L(t1) with classic Runge–Kutta method
    % evaluate L1 = L0 + h * (F(t0, U0 * L0')' * U0) with input function;
    Fl_tilde = @(L0) F_l(U0, L0);
    L1 = RK4(h, L0, Fl_tilde);
    [V_hat, ~] = qr([L1, V0], 0); % V_hat is the orth. basis of n×2r matrix
    N = V_hat' * V0; % 2r×r matrix N

    % S-step
    S0 = M * S0 * N'; % compute S_hat(t0) which equals S0 here
    % update S0 to S_hat(t1)
    % S_hat = S0 + h * (U_hat' * F(t0, U_hat * S0 * V_hat') * V_hat);
    Fs_tilde = @(S0) F_s(U_hat, S0, V_hat);
    S_hat = RK4(h, S0, Fs_tilde);

    % Truncation
    [P_hat, Sig, Q_hat] = svd(S_hat); % compute the SVD of S_hat(t1)
    sigma = diag(Sig); % diagonal elements of Sig
    r1 = 1; % set an initial value of new rank
    % choose the new rank r1 <= 2r, here is r = r0
    while sqrt(sum(sigma((r1+1):2*r0).^2)) > theta
        r1 = r1 + 1;
    end
    r1 = min(r1, rmax);

    % compute the new factors (after 1 time step) for the approx. of Y
    S1 = diag(sigma(1:r1));
    P1 = P_hat(:, 1:r1);
    Q1 = Q_hat(:, 1:r1);
    U1 = U_hat * P1; % m×r1 matrix
    V1 = V_hat * Q1; % n×r1 matrix

    % absolute error estimation for norms and energies
    Y1 = U1 * S1 * V1';
    [D, V_cos] = const(n);
    [err_norm1, err_energy1] = err(Y1, Y0, D, V_cos);
    err_norm(i) = err_norm1;
    err_energy(i) = err_energy1;

    % accept new variables for next iteration
    t0 = t0 + h; % our process is time independ but we record for time step
    U0 = U1;
    S0 = S1;
    V0 = V1;
    r0 = r1;
end
figure(1);
plot(0:n-1, r)
hold on;
figure(2);
plot(0:n-1, err_norm)
hold on;
figure(3);
plot(0:n-1, err_energy)
hold on;

% Implement the discrete Schrödinger equation as our example function
% evaluate the functions
% H[Y] = 1/2 * (D*Y + Y*D') + V_cos*Y*V_cos'
function Fs_tilde = F_s(U, S, V)
r = size(U, 2);
n = size(U, 1);
[D, V_cos] = const(n);
U_tilde = [U'*(1/2*D*U), U'*(1/2*U), U'*(V_cos*U)];
S_tilde = [S, S, S];
V_tilde = [V'*V, (D*V)'*V, (V_cos*V)'*V];
Fs_tilde = -1i*(U_tilde(:, 1:r) * S_tilde(:, 1:r) * V_tilde(:, 1:r) + ...
           U_tilde(:, r+1:2*r) * S_tilde(:, r+1:2*r) * V_tilde(:, r+1:2*r) + ...
           U_tilde(:, 2*r+1:3*r) * S_tilde(:, 2*r+1:3*r) * V_tilde(:, 2*r+1:3*r));

function Fk_tilde = F_k(K, V)
r = size(K, 2);
n = size(K, 1);
[D, V_cos] = const(n);
K_tilde = [1/2*D*K, 1/2*K, V_cos*K];
V_tilde = [V'*V, (D*V)'*V, (V_cos*V)'*V];
Fk_tilde = -1i*(K_tilde(:, 1:r) * V_tilde(:, 1:r) + ...
           K_tilde(:, r+1:2*r) * V_tilde(:, r+1:2*r) + ...
           K_tilde(:, 2*r+1:3*r) * V_tilde(:, 2*r+1:3*r));

function Fl_tilde = F_l(U, L)
r = size(U, 2);
n = size(U, 1);
[D, V_cos] = const(n);
U_tilde = [(1/2*D*U)'*U, (1/2*U)'*U, (V_cos*U)'*U];
L_tilde = [L, D*L, V_cos*L];
Fl_tilde = -1i*(L_tilde(:, 1:r) * U_tilde(:, 1:r) + ...
           L_tilde(:, r+1:2*r) * U_tilde(:, r+1:2*r) + ...
           L_tilde(:, 2*r+1:3*r) * U_tilde(:, 2*r+1:3*r));

% Input the constants for evaluation functions and F
function [D, V_cos] = const(n)
D = diag(2 * ones(1,n)) + diag(-1 * ones(1,n-1),1) + diag(-1 * ones(1,n-1),-1);
V_cos = diag(1 - cos(2 * pi * (-n/2 : n/2-1) / n));
