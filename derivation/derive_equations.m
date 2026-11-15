%% Equations of motion of the rotary double pendulum (Lagrange method)
% Method: Denavit-Hartenberg kinematics + geometric Jacobians +
% Christoffel symbols.  Requires the Symbolic Math Toolbox.
%
% q1 = phi (arm), q2 = theta1 (1st pendulum), q3 = theta2 (2nd pendulum,
% relative to the 1st).  Result:  M(q) ddq + C(q,dq) dq + G(q) = Q
clear; clc;

%% 1. Symbolic variables
syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real
syms m1 m2 m3 real          % masses (arm, pendulum 1, pendulum 2)
syms L0 L1 L2 real          % link lengths
syms d_ram d1 d2 real       % centre-of-mass offsets
syms g tau real             % gravity, motor torque
syms I1_xx I1_yy I1_zz I2_xx I2_yy I2_zz I3_xx I3_yy I3_zz real

I1 = diag([I1_xx, I1_yy, I1_zz]);
I2 = diag([I2_xx, I2_yy, I2_zz]);
I3 = diag([I3_xx, I3_yy, I3_zz]);
q = [q1; q2; q3];  dq = [dq1; dq2; dq3];  ddq = [ddq1; ddq2; ddq3];

%% 2. Kinematics - standard Denavit-Hartenberg
dh = @(theta, d, a, alpha) [
    cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
    sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
    0           sin(alpha)             cos(alpha)            d;
    0           0                      0                     1];

T01 = dh(q1, 0, L0, pi/2);   % arm: rotation about z0, offset L0, axis tilted 90 deg
T12 = dh(q2, 0, L1, 0);      % 1st pendulum
T23 = dh(q3, 0, L2, 0);      % 2nd pendulum
T02 = simplify(T01 * T12);
T03 = simplify(T02 * T23);

R01 = T01(1:3,1:3); R02 = T02(1:3,1:3); R03 = T03(1:3,1:3);
o0 = [0;0;0]; o1 = T01(1:3,4); o2 = T02(1:3,4);
z0 = [0;0;1]; z1 = T01(1:3,3); z2 = T02(1:3,3);

% Centres of mass (on the local x axis at distance d)
P_C1 = dh(q1, 0, d_ram, pi/2) * [0;0;0;1];        P_C1 = simplify(P_C1(1:3));
P_C2 = T01 * dh(q2, 0, d1, 0) * [0;0;0;1];        P_C2 = simplify(P_C2(1:3));
P_C3 = T02 * dh(q3, 0, d2, 0) * [0;0;0;1];        P_C3 = simplify(P_C3(1:3));

%% 3. Geometric Jacobians  (J_v = z x (P_C - o),  J_w = z)
Jv = @(z, pc, o) cross(z, pc - o);
Jv1 = [Jv(z0,P_C1,o0), [0;0;0], [0;0;0]];            Jw1 = [z0, [0;0;0], [0;0;0]];
Jv2 = [Jv(z0,P_C2,o0), Jv(z1,P_C2,o1), [0;0;0]];     Jw2 = [z0, z1, [0;0;0]];
Jv3 = [Jv(z0,P_C3,o0), Jv(z1,P_C3,o1), Jv(z2,P_C3,o2)]; Jw3 = [z0, z1, z2];

%% 4. Inertia matrix  M = sum( m Jv'Jv + Jw' R I R' Jw )
disp('Computing M(q)...');
M = simplify( m1*(Jv1.'*Jv1) + Jw1.'*(R01*I1*R01.')*Jw1 ...
            + m2*(Jv2.'*Jv2) + Jw2.'*(R02*I2*R02.')*Jw2 ...
            + m3*(Jv3.'*Jv3) + Jw3.'*(R03*I3*R03.')*Jw3 );

%% 5. Coriolis matrix from Christoffel symbols
disp('Computing C(q,dq)...');
C = sym(zeros(3));
for k = 1:3
    for j = 1:3
        for i = 1:3
            C(k,j) = C(k,j) + 0.5*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) ...
                                 - diff(M(i,j),q(k))) * dq(i);
        end
    end
end
C = simplify(C);

%% 6. Gravity vector  G = dV/dq
disp('Computing G(q)...');
V = g*(m1*P_C1(3) + m2*P_C2(3) + m3*P_C3(3));
G = simplify(jacobian(V, q).');

%% 7. Equations of motion
eqn = M*ddq + C*dq + G == [tau; 0; 0];
disp('M(q) ='); disp(M);
disp('C(q,dq) ='); disp(C);
disp('G(q) ='); disp(G);
for k = 1:3, fprintf('\nEquation %d:\n', k); disp(eqn(k)); end
