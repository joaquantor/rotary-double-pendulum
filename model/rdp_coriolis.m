function C = rdp_coriolis(y, p)
%RDP_CORIOLIS Coriolis/centrifugal matrix C(q,dq) (Christoffel symbols).
%   C_kj = sum_i 1/2 (dM_kj/dq_i + dM_ki/dq_j - dM_ij/dq_k) dq_i
%   Checked symbolically against rdp_mass_matrix (see
%   derivation/verify_coriolis.m).
q2 = y(2); q3 = y(3);
dq1 = y(4); dq2 = y(5); dq3 = y(6);

s1 = sin(q2); c1 = cos(q2); s2 = sin(q3);
s12 = sin(q2 + q3); c12 = cos(q2 + q3);

% lambda_i = dM11/dq_i  (change of the arm inertia with the pendulum angles)
H = p.L1*s1 + p.d2*s12;            % sideways offset of the 2nd pendulum CoM
lambda1 = p.m2*p.d1^2*sin(2*q2) + 2*p.m3*H*(p.L1*c1 + p.d2*c12) ...
        + (p.I2_yy - p.I2_xx)*sin(2*q2) + (p.I3_yy - p.I3_xx)*sin(2*(q2+q3));
lambda2 = 2*p.m3*H*p.d2*c12 + (p.I3_yy - p.I3_xx)*sin(2*(q2+q3));

h = -p.m3*p.L1*p.d2*s2;            % dM23/dq3, coupling of the two pendulums

% Terms from dM12/dq and dM13/dq (centripetal reaction of the pendulums on the arm)
k12 = -p.L0*((p.m2*p.d1 + p.m3*p.L1)*s1*dq2 + p.m3*p.d2*s12*(dq2 + dq3));
k13 = -p.L0*p.m3*p.d2*s12*(dq2 + dq3);

C = [0.5*(lambda1*dq2 + lambda2*dq3), 0.5*lambda1*dq1 + k12, 0.5*lambda2*dq1 + k13;
     -0.5*lambda1*dq1,                h*dq3,                 h*(dq2 + dq3);
     -0.5*lambda2*dq1,               -h*dq2,                 0];
end
