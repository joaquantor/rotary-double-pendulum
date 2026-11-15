%% Identification of the pendulum inertias and joint friction
% Data: free oscillations with the arm clamped
%   kyvadlo2fix  - 1st pendulum clamped too, only the 2nd one swings
%   kyvadlo2fixB - both pendulums swing
% Fitted: I2_zz, I3_zz, b_theta1, b_theta2, c_theta1, c_theta2
% Cost:   mean RMS error of the pendulum angles over 5 decay segments
clear; clc; close all;

p = rdp_params();
segs = {'kyvadlo2fix', 18.3, 22.5, true;   % {experiment, t_from, t_to, pend1 clamped}
        'kyvadlo2fix', 25.6, 30.0, true;
        'kyvadlo2fixB',  5,  15,  false;
        'kyvadlo2fixB', 31,  41,  false;
        'kyvadlo2fixB', 58,  68,  false};

for s = 1:size(segs,1)
    e = rdp_load_experiment(segs{s,1});
    th = [e.theta1 e.theta2];
    i = find(e.t >= segs{s,2} & e.t <= segs{s,3});
    if segs{s,4}                         % start at a peak -> zero velocity
        [~, j] = max(th(i,2)); i = i(j):i(end);
    end
    dth = gradient(smoothdata(th, 'gaussian', 8)', 0.005)';
    D(s).t = e.t(i) - e.t(i(1));                         %#ok<SAGROW>
    D(s).th = th(i,:);                                   %#ok<SAGROW>
    D(s).lock = segs{s,4};                               %#ok<SAGROW>
    if D(s).lock
        D(s).y0 = [th(i(1),2); dth(i(1),2)];             %#ok<SAGROW>
    else
        D(s).y0 = [th(i(1),:) dth(i(1),:)]';             %#ok<SAGROW>
    end
end

names = {'I2_zz','I3_zz','b_theta1','b_theta2','c_theta1','c_theta2'};
x0 = [7e-5 7e-5 3e-4 3e-4 5e-4 5e-4];          % initial guess = scaling
cost = @(z) segment_cost(z .* x0, names, p, D);
o = optimset('Display', 'iter', 'MaxFunEvals', 3000, 'MaxIter', 3000, ...
             'TolX', 1e-4, 'TolFun', 1e-6);
fprintf('Initial cost: %.4f rad\n', cost(ones(1,6)));
z = fminsearch(cost, ones(1,6), o);
z = fminsearch(cost, z, o);                     % restart helps Nelder-Mead
x = z .* x0;

fprintf('\nIdentified parameters (copy into model/rdp_params.m):\n');
for k = 1:numel(names), fprintf('  p.%-9s = %.4g;\n', names{k}, x(k)); end
fprintf('Mean RMS error: %.4f rad\n', cost(z));

% Plot the fit
pf = apply(p, names, x);
figure('Name', 'Pendulum identification', 'Color', 'w');
for s = 1:numel(D)
    y = rdp_locked_arm(pf, D(s).t, D(s).y0, D(s).lock);
    subplot(numel(D), 1, s);
    plot(D(s).t, D(s).th(:,2), 'r', D(s).t, y(:,2), 'r--'); hold on;
    if ~D(s).lock, plot(D(s).t, D(s).th(:,1), 'b', D(s).t, y(:,1), 'b--'); end
    grid on; ylabel('[rad]');
    title(sprintf('%s, segment %d', segs{s,1}, s), 'Interpreter', 'none');
end
legend('\theta_2 measured', '\theta_2 model', '\theta_1 measured', '\theta_1 model', ...
       'Location', 'northeastoutside');
xlabel('t [s]');

function J = segment_cost(x, names, p, D)
    if any(x < 0), J = 1e3; return; end          % friction/inertia must be >= 0
    p = apply(p, names, x);
    J = 0;
    for s = 1:numel(D)
        y = rdp_locked_arm(p, D(s).t, D(s).y0, D(s).lock);
        if isempty(y), J = 1e3; return; end
        if D(s).lock, err = y(:,2) - D(s).th(:,2); else, err = y - D(s).th; end
        J = J + sqrt(mean(err(:).^2));
    end
    J = J / numel(D);
end

function p = apply(p, names, x)
    for k = 1:numel(names), p.(names{k}) = x(k); end
    p.I2_yy = p.I2_zz; p.I3_yy = p.I3_zz;
end
