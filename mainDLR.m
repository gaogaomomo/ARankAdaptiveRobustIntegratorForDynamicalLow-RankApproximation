% Input the vatiables for function dlr
t0 = 0;
T = 1;
n = 100;
r0 = [4, 8, 12];
theta = 10e-6;
rmax = 30;
h = (T - t0) / n;
for i = 1:3
    [U0, ~] = qr(randn(n,r0(i)),0);
    S0 = diag(10.^(-(1 : r0(i))));
    [V0, ~] = qr(randn(n,r0(i)),0);
    [U, S, V] = dlr(t0, h, U0, S0, V0, r0(i), theta, n,rmax);
end
figure(1); % rank
xlabel('time'); % Label for x-axis
ylabel('rank'); % Label for y-axis
legend({'r0 = 4', 'r0 = 8', 'r0 = 12'}, 'Location', 'best'); % Legend
title('Plot of Three Different r0'); % Title
grid on; % Turn on the grid

figure(2); % error of norms
xlabel('time'); % Label for x-axis
ylabel('error'); % Label for y-axis
legend({'r0 = 4', 'r0 = 8', 'r0 = 12'}, 'Location', 'best'); % Legend
title('Plot of Error of Norms'); % Title
grid on; % Turn on the grid

figure(3); % error of energies
xlabel('time'); % Label for x-axis
ylabel('error'); % Label for y-axis
legend({'r0 = 4', 'r0 = 8', 'r0 = 12'}, 'Location', 'best'); % Legend
title('Plot of Error of Energies'); % Title
grid on; % Turn on the grid