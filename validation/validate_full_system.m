%% Validation of the full model against the measurement 'kyvadlo2'
% The arm and both pendulums are free; the arm is driven by a sequence of
% voltage steps. The model is driven by the same measured voltage and
% starts from the measured initial angles.
clear; clc; close all;

p = rdp_params();
e = rdp_load_experiment('kyvadlo2');

fprintf('Simulating %.1f s ...\n', e.t(end));
tic; [t_sim, y_sim] = rdp_simulate(p, e); toc;

err_rms = @(a, b) sqrt(mean((a - b).^2));
early = e.t < 4;
fprintf('RMS error, whole record: phi = %.3f, theta1 = %.3f, theta2 = %.3f rad\n', ...
        err_rms(y_sim(:,1), e.phi), err_rms(y_sim(:,2), e.theta1), err_rms(y_sim(:,3), e.theta2));
fprintf('RMS error, first 4 s:    theta1 = %.3f, theta2 = %.3f rad\n', ...
        err_rms(y_sim(early,2), e.theta1(early)), err_rms(y_sim(early,3), e.theta2(early)));

figure('Name', 'Model validation - kyvadlo2', 'Color', 'w', 'Position', [100 100 900 800]);
ax(1) = subplot(4,1,1);
plot(e.t, e.u, 'Color', [0 0.6 0]); ylabel('u [V]'); title('Input voltage'); grid on;
lbl = {'\phi [rad]', '\theta_1 [rad]', '\theta_2 [rad]'};
ttl = {'Arm', '1st pendulum', '2nd pendulum (relative)'};
meas = [e.phi e.theta1 e.theta2];
for k = 1:3
    ax(k+1) = subplot(4,1,k+1); %#ok<SAGROW>
    plot(e.t, meas(:,k), 'b--', t_sim, y_sim(:,k), 'r', 'LineWidth', 1.2);
    ylabel(lbl{k}); title(ttl{k}); grid on;
    legend('measured', 'model', 'Location', 'best');
end
xlabel('t [s]');
linkaxes(ax, 'x');
