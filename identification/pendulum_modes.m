%% Linearized modes of the double pendulum with the arm clamped
% Compares the natural periods and mode shapes of the model with the
% measurement 'kyvadlo2fixB'. This check revealed that the 2nd pendulum
% encoder measures the ABSOLUTE angle of the 2nd link: the measured ratio
% theta_2,abs / theta_1 of the slow mode matches the model, while the
% relative angle does not.
clear; clc;

p = rdp_params();
y0 = zeros(6,1);
M = rdp_mass_matrix(y0, p);  M = M(2:3, 2:3);
K = p.g * [p.m2*p.d1 + p.m3*p.L1 + p.m3*p.d2, p.m3*p.d2;
           p.m3*p.d2,                         p.m3*p.d2];
[V, W] = eig(M \ K);
[w2, idx] = sort(diag(W)); V = V(:, idx);

fprintf('Model (arm clamped):\n');
for k = 1:2
    fprintf('  mode %d: T = %.3f s, theta2_rel/theta1 = %+.2f, theta2_abs/theta1 = %+.2f\n', ...
            k, 2*pi/sqrt(w2(k)), V(2,k)/V(1,k), 1 + V(2,k)/V(1,k));
end
fprintf('  2nd pendulum alone: T = %.3f s\n', ...
        2*pi*sqrt((p.m3*p.d2^2 + p.I3_zz) / (p.m3*p.g*p.d2)));

% Measurement
e = rdp_load_experiment('kyvadlo2fixB');
i = e.t > 8 & e.t < 20;
th1 = e.theta1(i); th2abs = e.theta2(i) + th1;
[~, loc] = findpeaks(th1, e.t(i), 'MinPeakProminence', 0.05);
fprintf('\nMeasured (kyvadlo2fixB, 8-20 s):\n');
fprintf('  T = %.3f s, theta2_abs/theta1 = %+.2f (least squares)\n', ...
        mean(diff(loc)), th1 \ th2abs);

e = rdp_load_experiment('kyvadlo2fix');
i = e.t > 18.3 & e.t < 22;
[~, loc] = findpeaks(e.theta2(i) + e.theta1(i), e.t(i), 'MinPeakProminence', 0.02);
fprintf('Measured (kyvadlo2fix, 2nd pendulum alone): T = %.3f s\n', mean(diff(loc)));
