%% Verify the hand-written Coriolis matrix against M(q)
% Builds M(q) symbolically from model/rdp_mass_matrix.m, computes C(q,dq)
% from Christoffel symbols and compares it with model/rdp_coriolis.m.
% Also checks the skew-symmetry of dM/dt - 2C (passivity property).
% Requires the Symbolic Math Toolbox.
clear; clc;

syms q1 q2 q3 dq1 dq2 dq3 real
names = {'m1','m2','m3','L0','L1','d_ram','d1','d2','I1_zz', ...
         'I2_xx','I2_yy','I2_zz','I3_xx','I3_yy','I3_zz'};
for k = 1:numel(names), p.(names{k}) = sym(names{k}, 'real'); end

q = [q1; q2; q3]; dq = [dq1; dq2; dq3];
y = [q; dq];

M = rdp_mass_matrix(y, p);
C_ref = sym(zeros(3));
for k = 1:3
    for j = 1:3
        for i = 1:3
            C_ref(k,j) = C_ref(k,j) + 0.5*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) ...
                                         - diff(M(i,j),q(k))) * dq(i);
        end
    end
end

C = rdp_coriolis(y, p);
D = simplify(expand(C - C_ref), 'Steps', 100);
dM = sym(zeros(3));
for i = 1:3, dM = dM + diff(M, q(i))*dq(i); end
N = dM - 2*C;
S = simplify(expand(N + N.'), 'Steps', 100);

fprintf('C equals the Christoffel-symbol C:   %d\n', isequal(D, sym(zeros(3))));
fprintf('dM/dt - 2C is skew-symmetric:        %d\n', isequal(S, sym(zeros(3))));
