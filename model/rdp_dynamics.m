function dydt = rdp_dynamics(t, y, p, u_fun)
%RDP_DYNAMICS State derivative of the rotary double pendulum.
%   dydt = RDP_DYNAMICS(t, y, p, u_fun)
%     y     = [phi; theta1; theta2; dphi; dtheta1; dtheta2]
%     p     = parameter struct (rdp_params)
%     u_fun = function handle u(t) returning the control voltage [V]
%             (e.g. a griddedInterpolant of the measured input)
%
%   Equations of motion (Lagrange, 2nd kind):
%       M(q) ddq + C(q,dq) dq + G(q) = [tau; 0; 0] - D dq - F_c(dq)

q2 = y(2); q3 = y(3);
dq = y(4:6);

% --- Motor: control voltage -> torque on the arm ------------------------
u_in  = max(min(u_fun(t), p.U_ctrl_max), -p.U_ctrl_max);
U_app = u_in / p.U_ctrl_max * p.U_pwm_amp;
I_a   = (U_app - p.Ke_out * dq(1)) / p.R;   % armature current (L neglected)
tau   = p.Kt_out * I_a;

M = rdp_mass_matrix(y, p);
C = rdp_coriolis(y, p);

s1 = sin(q2); s12 = sin(q2 + q3);
G = [0;
     p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
     p.g *   p.m3*p.d2*s12];

% Viscous + (smoothed) Coulomb friction
F = [p.b_phi; p.b_theta1; p.b_theta2] .* dq ...
  + [0; p.c_theta1; p.c_theta2] .* tanh(dq / p.v_eps);

ddq  = M \ ([tau; 0; 0] - F - C*dq - G);
dydt = [dq; ddq];
end
