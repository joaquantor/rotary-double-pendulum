%% Motor + arm without pendulums (experiment 'bezzateze')
% 1-DOF check of the DC motor model:
%   J dd(phi) = Kt_out/R * (U_app - Ke_out dphi) - b_phi dphi
clear; clc; close all;

p = rdp_params();
e = rdp_load_experiment('bezzateze');
e.dphi = smoothdata(gradient(e.phi, e.t), 'gaussian', 20);   % measured velocity

J = p.I1_zz;   % arm only (pendulums removed)
u_fun = griddedInterpolant(e.t, e.u, 'previous', 'nearest');
f = @(t, y) [y(2);
    (p.Kt_out*((max(min(u_fun(t), p.U_ctrl_max), -p.U_ctrl_max)/p.U_ctrl_max*p.U_pwm_amp ...
      - p.Ke_out*y(2))/p.R) - p.b_phi*y(2)) / J];
[t_sim, y_sim] = ode45(f, e.t, [e.phi(1); 0], odeset('RelTol', 1e-5, 'AbsTol', 1e-6));

fprintf('RMS error: phi = %.3f rad, dphi = %.3f rad/s\n', ...
        sqrt(mean((y_sim(:,1) - e.phi).^2)), sqrt(mean((y_sim(:,2) - e.dphi).^2)));

figure('Name', 'Arm without load', 'Color', 'w');
ax(1) = subplot(3,1,1); plot(e.t, e.u, 'g'); ylabel('u [V]'); title('Input voltage'); grid on;
ax(2) = subplot(3,1,2); plot(e.t, e.phi, 'b--', t_sim, y_sim(:,1), 'r');
ylabel('\phi [rad]'); title('Arm angle'); legend('measured', 'model'); grid on;
ax(3) = subplot(3,1,3); plot(e.t, e.dphi, 'b--', t_sim, y_sim(:,2), 'r');
ylabel('d\phi/dt [rad/s]'); title('Arm velocity'); xlabel('t [s]'); grid on;
linkaxes(ax, 'x');
