%% Free simulation with a 3-D animation
% Simulates the identified model from a chosen initial state with an
% optional voltage pulse train and animates the result in real time.
clear; clc; close all;

p = rdp_params();

% --- Scenario -------------------------------------------------------------
q0  = [0; pi - 0.1; 0.5];   % [phi; theta1; theta2]  (0 = hanging down)
dq0 = [0; 0; 0];
T_end = 20;
pulse.amp = 0;  pulse.period = 2;  pulse.width = 0.05;  % set amp ~= 0 to drive the arm
animate = true;

u_fun = @(t) pulse.amp * (mod(t, pulse.period) < pulse.width);
[t, y] = ode45(@(t, y) rdp_dynamics(t, y, p, u_fun), [0 T_end], [q0; dq0], ...
               odeset('RelTol', 1e-6, 'AbsTol', 1e-8));

% Energy check (without input the total energy can only decrease)
E = zeros(size(t));
for k = 1:numel(t)
    M = rdp_mass_matrix(y(k,:)', p);
    V = -p.g*(p.m2*p.d1*cos(y(k,2)) + p.m3*(p.L1*cos(y(k,2)) + p.d2*cos(y(k,2)+y(k,3))));
    E(k) = 0.5*y(k,4:6)*M*y(k,4:6)' + V;
end

figure('Name', 'Simulation', 'Color', 'w');
subplot(2,1,1); plot(t, y(:,1:3)); grid on; ylabel('[rad]');
legend('\phi', '\theta_1', '\theta_2'); title('Generalized coordinates');
subplot(2,1,2); plot(t, E); grid on; ylabel('E [J]'); xlabel('t [s]');
title('Total mechanical energy');

if animate, animate_rdp(t, y, p); end

function animate_rdp(t, y, p)
f = figure('Name', 'Rotary double pendulum', 'Color', 'w');
ax = axes(f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal'); view(ax, 45, 30);
lim = p.L0 + p.L1 + p.L2 + 0.05;
xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); zlim(ax, [-lim lim]);
plot3(ax, [0 0], [0 0], [-0.1 0], 'k-', 'LineWidth', 4);
h_arm  = plot3(ax, 0, 0, 0, 'k-o', 'LineWidth', 3, 'MarkerFaceColor', 'k');
h_pen1 = plot3(ax, 0, 0, 0, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
h_pen2 = plot3(ax, 0, 0, 0, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
tic;
while ishandle(f) && toc <= t(end)
    q = interp1(t, y(:,1:3), toc);
    x1 = p.L0*cos(q(1)); y1 = p.L0*sin(q(1));
    r1 = p.L1*sin(q(2));                    r2 = r1 + p.L2*sin(q(2)+q(3));
    z2 = -p.L1*cos(q(2));                   z3 = z2 - p.L2*cos(q(2)+q(3));
    x2 = x1 - r1*sin(q(1)); y2 = y1 + r1*cos(q(1));
    x3 = x1 - r2*sin(q(1)); y3 = y1 + r2*cos(q(1));
    set(h_arm,  'XData', [0 x1],  'YData', [0 y1],  'ZData', [0 0]);
    set(h_pen1, 'XData', [x1 x2], 'YData', [y1 y2], 'ZData', [0 z2]);
    set(h_pen2, 'XData', [x2 x3], 'YData', [y2 y3], 'ZData', [z2 z3]);
    drawnow;
end
end
