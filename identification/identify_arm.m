%% Identification of the arm friction and centre of mass
% Data: kyvadlo2 (whole system driven by a voltage sequence)
% Fitted: b_phi (viscous friction of the arm), d_ram (arm CoM offset)
% Cost:   RMS error of the arm angle phi over the whole record
% Run identify_pendulums.m first - the pendulum parameters affect M11.
clear; clc; close all;

p = rdp_params();
e = rdp_load_experiment('kyvadlo2');
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);

cost = @(x) arm_cost(x, p, e, opts);
fprintf('Current parameters: RMS = %.4f rad\n', cost([p.b_phi p.d_ram]));
x = fminsearch(cost, [0.1 -0.05], optimset('Display', 'iter', 'TolX', 1e-5, 'TolFun', 1e-6));

fprintf('\nIdentified parameters (copy into model/rdp_params.m):\n');
fprintf('  p.b_phi = %.4f;\n  p.d_ram = %.4f;\n', x(1), x(2));
fprintf('RMS error of phi: %.4f rad\n', cost(x));

function J = arm_cost(x, p, e, opts)
    if x(1) <= 0, J = 1e6; return; end
    p.b_phi = x(1); p.d_ram = x(2);
    [~, y] = rdp_simulate(p, e, opts);
    if size(y,1) ~= numel(e.t), J = 1e6; return; end
    J = sqrt(mean((y(:,1) - e.phi).^2));
end
