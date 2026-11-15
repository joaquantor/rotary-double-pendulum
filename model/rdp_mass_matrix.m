function M = rdp_mass_matrix(y, p)
%RDP_MASS_MATRIX Inertia matrix M(q) of the rotary double pendulum.
%   Derived in derivation/derive_equations.m (Denavit-Hartenberg kinematics,
%   geometric Jacobians, T = 1/2 dq' M(q) dq).
q2 = y(2); q3 = y(3);
c1 = cos(q2); s1 = sin(q2);
c2 = cos(q3);
c12 = cos(q2 + q3); s12 = sin(q2 + q3);

% Arm inertia grows as the pendulums swing out sideways
M11 = p.I1_zz + p.m1*p.d_ram^2 ...
    + p.m2*(p.L0^2 + (p.d1*s1)^2) ...
    + p.m3*(p.L0^2 + (p.L1*s1 + p.d2*s12)^2) ...
    + p.I2_yy*s1^2 + p.I2_xx*c1^2 + p.I3_yy*s12^2 + p.I3_xx*c12^2;

M22 = p.m2*p.d1^2 + p.m3*(p.L1^2 + p.d2^2 + 2*p.L1*p.d2*c2) + p.I2_zz + p.I3_zz;
M33 = p.m3*p.d2^2 + p.I3_zz;
M23 = p.m3*p.d2^2 + p.m3*p.L1*p.d2*c2 + p.I3_zz;
M12 = p.m2*p.L0*p.d1*c1 + p.m3*p.L0*(p.L1*c1 + p.d2*c12);
M13 = p.m3*p.L0*p.d2*c12;

M = [M11 M12 M13;
     M12 M22 M23;
     M13 M23 M33];
end
