function y = rdp_locked_arm(p, t, y0, lock_pend1)
%RDP_LOCKED_ARM Free motion of the pendulums with the arm clamped.
%   y = RDP_LOCKED_ARM(p, t, y0)  -> both pendulums free,
%       y0 = [theta1; theta2; dtheta1; dtheta2], y = [theta1 theta2]
%   y = RDP_LOCKED_ARM(p, t, y0, true) -> 1st pendulum also clamped,
%       only the 2nd swings, y0 = [theta2; dtheta2]
%   Returns [] if the integration fails (used inside optimizers).
if nargin < 4, lock_pend1 = false; end
o = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

if lock_pend1
    J = p.m3*p.d2^2 + p.I3_zz;   % inertia about the joint
    f = @(~, z) [z(2);
        (-p.g*p.m3*p.d2*sin(z(1)) - p.b_theta2*z(2) - p.c_theta2*tanh(z(2)/p.v_eps)) / J];
    [~, z] = ode45(f, t, y0, o);
    y = [zeros(size(z,1),1) z(:,1)];
else
    [~, z] = ode45(@(~, z) two_link(z, p), t, y0, o);
    y = z(:,1:2);
end
if size(y,1) ~= numel(t), y = []; end
end

function dz = two_link(z, p)
% Lower-right 2x2 block of the full model with phi = dphi = 0
yfull = [0; z(1); z(2); 0; z(3); z(4)];
M = rdp_mass_matrix(yfull, p);  M = M(2:3, 2:3);
C = rdp_coriolis(yfull, p);     C = C(2:3, 2:3);
s1 = sin(z(1)); s12 = sin(z(1) + z(2));
G = p.g * [(p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12; p.m3*p.d2*s12];
dq = z(3:4);
F = [p.b_theta1; p.b_theta2].*dq + [p.c_theta1; p.c_theta2].*tanh(dq/p.v_eps);
dz = [dq; M \ (-C*dq - G - F)];
end
