function e = rdp_load_experiment(name)
%RDP_LOAD_EXPERIMENT Load one experiment from data/experimenty.mat in SI units.
%   e = RDP_LOAD_EXPERIMENT('kyvadlo2') returns a struct with fields
%     t      time [s] (starting at 0)
%     u      control voltage [V]
%     phi    arm angle [rad]
%     theta1 1st pendulum angle [rad]
%     theta2 2nd pendulum angle RELATIVE to the 1st one [rad]
%
%   The raw signals are encoder counts (x4 quadrature decoding):
%     arm:       500 PPR on the motor shaft, gearbox 128/14
%     pendulums: 1024 PPR
%   NOTE: the 2nd pendulum encoder reads the ABSOLUTE angle of the 2nd link
%   (identified from the mode shapes in 'kyvadlo2fixB'), so the relative
%   angle is obtained as theta_abs - theta1.
%
%   Available experiments: bezzateze, sezavazim, otacky, kyvadlo1vol,
%   kyvadlo1, kyvadlo2fix, kyvadlo2fixB, kyvadlo2vol, kyvadlo2.

S = load('experimenty.mat', name);
d = S.(name);
v = double(d.signals.values);

k_phi   = 2*pi / (4 * 500 * 128/14);
k_theta = 2*pi / (4 * 1024);

e.name   = name;
e.t      = d.time(:) - d.time(1);
e.u      = v(:,1);
e.phi    = v(:,2) * k_phi;
e.theta1 = v(:,3) * k_theta;
e.theta2 = v(:,4) * k_theta - e.theta1;
end
