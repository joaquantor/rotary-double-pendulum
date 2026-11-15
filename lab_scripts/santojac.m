%% Odvození pohybových rovnic rotačního dvojitého kyvadla (Lagrangeova metoda)
% Metoda: Denavit-Hartenberg + Jakobiány + Christoffelovy symboly
clear; clc;

%% 1. Definice symbolických proměnných
% q1 = phi (rotace ramene), q2 = theta1 (1. kyvadlo), q3 = theta2 (2. kyvadlo)
syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real
syms m1 m2 m3 real      % Hmotnosti (1=rameno, 2=kyv1, 3=kyv2)
syms L0 L1 L2 real      % Délky (L0=rameno, L1=kyv1, L2=kyv2)
syms d_ram d1 d2 real   % Vzdálenosti těžišť od osy (d_ram, d1, d2)
syms g real             % Gravitace
syms tau real           % Krouticí moment motoru

% Tenzory setrvačnosti (uvažujeme diagonální v lokálních souřadnicích)
% Ixx, Iyy, Izz pro každé těleso
syms I1_xx I1_yy I1_zz real
syms I2_xx I2_yy I2_zz real
syms I3_xx I3_yy I3_zz real

I1 = diag([I1_xx, I1_yy, I1_zz]);
I2 = diag([I2_xx, I2_yy, I2_zz]);
I3 = diag([I3_xx, I3_yy, I3_zz]);

% Vektory souřadnic
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

%% 2. Kinematika - Denavit-Hartenberg (DH)
% Funkce pro DH matici
dh = @(theta, d, a, alpha) [
    cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
    sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
    0           sin(alpha)             cos(alpha)            d;
    0           0                      0                     1];

% Transformační matice (Standard DH)
% 1. Rameno (rotace kolem Z, posun po X o L0, otočení osy o 90°)
T01 = dh(q1, 0, L0, pi/2); 

% 2. První kyvadlo (rotace kolem Z', posun po X' o L1)
T12 = dh(q2, 0, L1, 0);

% 3. Druhé kyvadlo (rotace kolem Z'', posun po X'' o L2)
T23 = dh(q3, 0, L2, 0);

% Celkové transformace do globálního rámu (0)
T02 = simplify(T01 * T12);
T03 = simplify(T02 * T23);

% Rotační matice (R) - horní levá 3x3
R01 = T01(1:3, 1:3);
R02 = T02(1:3, 1:3);
R03 = T03(1:3, 1:3);

% Polohy počátků souřadných systémů
o0 = [0; 0; 0];
o1 = T01(1:3, 4);
o2 = T02(1:3, 4);

% Osy rotace (osy Z lokálních systémů vyjádřené v nule)
z0 = [0; 0; 1];         % Osa motoru
z1 = T01(1:3, 3);       % Osa 1. kyvadla
z2 = T02(1:3, 3);       % Osa 2. kyvadla

% Polohy těžišť (Centers of Mass - CoM) v globálním systému
% Předpoklad: Těžiště leží na ose x lokálního systému ve vzdálenosti d
P_C1 = simplify(dh(q1, 0, d_ram, pi/2) * [0;0;0;1]); P_C1 = P_C1(1:3);
P_C2 = simplify(T01 * dh(q2, 0, d1, 0) * [0;0;0;1]); P_C2 = P_C2(1:3);
P_C3 = simplify(T02 * dh(q3, 0, d2, 0) * [0;0;0;1]); P_C3 = P_C3(1:3);

%% 3. Sestavení Jakobiánů
% Funkce pro rotační Jakobián (osa rotace) a translační (z x r)
calc_Jv = @(z, p_com, o_joint) cross(z, p_com - o_joint);

% --- Těleso 1 (Rameno) ---
Jv1 = [calc_Jv(z0, P_C1, o0), [0;0;0], [0;0;0]];
Jw1 = [z0, [0;0;0], [0;0;0]];

% --- Těleso 2 (1. kyvadlo) ---
Jv2 = [calc_Jv(z0, P_C2, o0), calc_Jv(z1, P_C2, o1), [0;0;0]];
Jw2 = [z0, z1, [0;0;0]];

% --- Těleso 3 (2. kyvadlo) ---
Jv3 = [calc_Jv(z0, P_C3, o0), calc_Jv(z1, P_C3, o1), calc_Jv(z2, P_C3, o2)];
Jw3 = [z0, z1, z2];

%% 4. Matice setrvačnosti M(q)
% M = sum( m*Jv'*Jv + Jw'*R*I*R'*Jw )
disp('Počítám matici setrvačnosti M(q)...');

M1 = m1*(Jv1.'*Jv1) + Jw1.' * (R01*I1*R01.') * Jw1;
M2 = m2*(Jv2.'*Jv2) + Jw2.' * (R02*I2*R02.') * Jw2;
M3 = m3*(Jv3.'*Jv3) + Jw3.' * (R03*I3*R03.') * Jw3;

M = simplify(M1 + M2 + M3);

%% 5. Matice Coriolisových a odstředivých sil C(q, dq)
% Výpočet pomocí Christoffelových symbolů
disp('Počítám matici C(q, dq) pomocí Christoffelových symbolů...');

n = 3; % Počet stupňů volnosti
C = sym(zeros(n,n));

for k = 1:n
    for j = 1:n
        % c_kj = sum( c_ijk * dq_i )
        temp_sum = 0;
        for i = 1:n
            % Christoffelův symbol 1. druhu: c_ijk
            c_ijk = 0.5 * (diff(M(k,j), q(i)) + diff(M(k,i), q(j)) - diff(M(i,j), q(k)));
            temp_sum = temp_sum + c_ijk * dq(i);
        end
        C(k,j) = temp_sum;
    end
end
C = simplify(C);

%% 6. Vektor gravitačních sil G(q)
% V = m * g * h (z-ová souřadnice těžiště)
disp('Počítám gravitační vektor G(q)...');

% Z-ové souřadnice (3. prvek vektoru)
z_C1 = P_C1(3);
z_C2 = P_C2(3);
z_C3 = P_C3(3);

Ep = m1*g*z_C1 + m2*g*z_C2 + m3*g*z_C3;
G = jacobian(Ep, q).'; % Transpozice na sloupcový vektor
G = simplify(G);

%% 7. Sestavení finálních rovnic pohybu
% M(q)*ddq + C(q,dq)*dq + G(q) = Q
disp('Sestavuji finální rovnice...');

% Vektor vnějších sil (moment pouze na prvním kloubu)
Q_vec = [tau; 0; 0]; 

% Levá strana rovnice
LHS = M*ddq + C*dq + G;

% Výpis výsledků
disp('--------------------------------------------------');
disp('Matice setrvačnosti M(q):');
disp(M);
disp('--------------------------------------------------');
disp('Matice Coriolisových sil C(q, dq):');
disp(C);
disp('--------------------------------------------------');
disp('Vektor gravitace G(q):');
disp(G);
disp('--------------------------------------------------');
disp('Výsledné Lagrangeovy rovnice pohybu (LHS = Q):');
eqn = LHS == Q_vec;

% Pro zobrazení konkrétní rovnice pro q1 (phi):
disp('Rovnice pro phi (q1):');
disp(eqn(1));

% Pro zobrazení konkrétní rovnice pro theta1 (q2):
disp('Rovnice pro theta1 (q2):');
disp(eqn(2));

% Pro zobrazení konkrétní rovnice pro theta2 (q3):
disp('Rovnice pro theta2 (q3):');
disp(eqn(3));