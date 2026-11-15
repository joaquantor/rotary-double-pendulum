%% Porovnání simulace s naměřenými daty (MATLAB)
clear; clc; close all;

%% 1. NAČTENÍ DAT A PARAMETRŮ
load('experimenty.mat');

% Vybereme experiment "kyvadlo2" (Hnané obě kyvadla i základna)
% Struktura v .mat souboru má obvykle formát: time a signals.values
data_exp = kyvadlo2; 

t_meas = data_exp.time - data_exp.time(1);
vals = double(data_exp.signals.values);

% Data jsou SUROVÉ PULZY enkodérů -> přepočet na radiány
ENC_PPR_PEND = 1024;          % kyvadla: 1024 PPR -> 4096 hran/ot
ENC_PPR_MOT  = 500;           % motor:   500 PPR  -> 2000 hran/ot
GEAR_RATIO   = 128 / 14;
k_phi   = (2*pi) / (4 * ENC_PPR_MOT * GEAR_RATIO);
k_theta = (2*pi) / (4 * ENC_PPR_PEND);

% Rozdělení signálů (Napětí, Základna, Kyv1, Kyv2)
u_meas  = vals(:, 1);              % Vstupní napětí [V]
q1_meas = vals(:, 2) * k_phi;      % Základna (phi) [rad]
q2_meas = vals(:, 3) * k_theta;    % Kyvadlo 1 (theta1) [rad]
% POZOR: enkodér 2. kyvadla měří ABSOLUTNÍ úhel (ověřeno z tvaru módů
% v kyvadlo2fixB) -> relativní theta2 = abs - theta1
q3_meas = vals(:, 4) * k_theta - q2_meas;  % Kyvadlo 2 (theta2, relativní) [rad]

% --- PARAMETRY SYSTÉMU (Váš původní kód) ---
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;       
p.m2 = 0.0356 + 0.0059;  
p.m3 = 0.0435;  
p.L0 = 0.201; p.L1 = 0.163; p.L2 = 0.150;   
p.d_ram = -0.0638; p.d1 = 0.078; p.d2 = 0.075;    
p.gear_ratio = 128 / 14; 
p.I1_zz = 0.0051; 
% Původní I_zz (0.000318/0.000326) byly momenty k ose kloubu, ne k těžišti
% -> s m*d^2 se počítaly 2x. Nové hodnoty: fit na kyvadlo2fix + kyvadlo2fixB.
% I_yy (také kolmo k ose tyče) bereme stejné jako I_zz.
p.I2_xx = 0.00005; p.I2_zz = 1.283e-4; p.I2_yy = p.I2_zz;
p.I3_xx = 0.00005; p.I3_zz = 8.91e-5;  p.I3_yy = p.I3_zz;
p.g = 9.81;     
p.U_pwm_amp = 14.4; p.U_ctrl_max = 10.0; 
p.R = 2.20; p.Kt = 0.052; p.Ke = 0.0516; 
p.Ke_out = p.Ke * p.gear_ratio; 
p.Kt_out = p.Kt * p.gear_ratio;
p.b_phi = 0.1714;      % tlumení ramene [Nms/rad] (fit na phi z kyvadlo2)
p.b_theta1 = 0;        p.b_theta2 = 1.705e-4; % viskózní tření kloubů [Nms/rad]
p.c_theta1 = 5.70e-4;  p.c_theta2 = 4.06e-4;  % Coulombovo tření kloubů [Nm] (fit)
p.v_eps = 0.02;        % vyhlazení sign() pomocí tanh(dq/v_eps) [rad/s]

%% 2. NASTAVENÍ SIMULACE
% Počáteční podmínky převezmeme z prvního naměřeného vzorku
% y = [phi; theta1; theta2; dphi; dtheta1; dtheta2]
y0 = [q1_meas(1); q2_meas(1); q3_meas(1); 0; 0; 0];

t_span = [t_meas(1) t_meas(end)];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);

% Interpolant vstupu vytvoříme jednou (rychlejší než interp1 v každém kroku).
% Vstup z Simulinku je po částech konstantní (ZOH) -> 'previous'.
u_fun = griddedInterpolant(t_meas, u_meas, 'previous', 'nearest');
ode_fun = @(t, y) dynamics_measured(t, y, p, u_fun);

fprintf('Spouštím simulaci...\n');
tic;
[t_sim, y_sim] = ode45(ode_fun, t_span, y0, options);
toc;

%% 3. VYKRESLENÍ POROVNÁNÍ
figure('Name', 'Validace modelu vs Experiment', 'Color', 'w');

% 1. Napětí (Kontrola vstupu)
subplot(4,1,1);
plot(t_meas, u_meas, 'g'); hold on;
ylabel('U [V]'); title('Vstupní napětí (z experimentu)'); grid on;
legend('Experiment');

% 2. Rameno (Base)
subplot(4,1,2);
plot(t_meas, q1_meas, 'b--'); hold on;
plot(t_sim, y_sim(:,1), 'r', 'LineWidth', 1.5);
ylabel('\phi [rad]'); title('Základna (Rameno)'); grid on;
legend('Měření', 'Simulace');

% 3. První kyvadlo
subplot(4,1,3);
plot(t_meas, q2_meas, 'b--'); hold on;
plot(t_sim, y_sim(:,2), 'r', 'LineWidth', 1.5);
ylabel('\theta_1 [rad]'); title('1. Kyvadlo'); grid on;
legend('Měření', 'Simulace');

% 4. Druhé kyvadlo
subplot(4,1,4);
plot(t_meas, q3_meas, 'b--'); hold on;
plot(t_sim, y_sim(:,3), 'r', 'LineWidth', 1.5);
ylabel('\theta_2 [rad]'); title('2. Kyvadlo'); grid on;
legend('Měření', 'Simulace');
xlabel('Čas [s]');
linkaxes(findobj(gcf, 'Type', 'axes'), 'x');

%% 4. POMOCNÉ FUNKCE (Dynamika)

function dydt = dynamics_measured(t, y, p, u_fun)
    % 1. Vstup (řídicí napětí) v aktuálním čase simulace
    u_in = u_fun(t);
    
    % Omezení napětí (saturace)
    u_in = max(min(u_in, p.U_ctrl_max), -p.U_ctrl_max);
    
    % Rozbalení stavů (ZDE JE ZMĚNA OPROTI PŮVODNÍMU KÓDU)
    % Nyní simulujeme i q1 (základnu), nefixujeme ji na 0
    q1 = y(1); q2 = y(2); q3 = y(3);       
    dq1 = y(4); dq2 = y(5); dq3 = y(6);    
    dq = [dq1; dq2; dq3];
    
    % --- Model Motoru ---
    U_app = (u_in / p.U_ctrl_max) * p.U_pwm_amp;
    % Back-EMF závisí na rychlosti ramene (dq1)
    Current = (U_app - p.Ke_out * dq1) / p.R; 
    Tau = p.Kt_out * Current;
    
    % --- MATICE M(q) ---
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2+q3); s12 = sin(q2+q3);
    
    term2 = p.m2 * (p.L0^2 + (p.d1*s1)^2);
    term3 = p.m3 * (p.L0^2 + (p.L1*s1 + p.d2*s12)^2);
    rot_inertia = p.I2_yy*s1^2 + p.I2_xx*c1^2 + p.I3_yy*s12^2 + p.I3_xx*c12^2;
    
    M11 = p.I1_zz + p.m1*p.d_ram^2 + term2 + term3 + rot_inertia;
    
    M22 = p.m2*p.d1^2 + p.m3*(p.L1^2 + p.d2^2 + 2*p.L1*p.d2*c2) + p.I2_zz + p.I3_zz;
    M33 = p.m3*p.d2^2 + p.I3_zz;
    M23 = p.m3*p.d2^2 + p.m3*p.L1*p.d2*c2 + p.I3_zz;
    M32 = M23;
    
    M12 = p.m2*p.L0*p.d1*c1 + p.m3*p.L0*(p.L1*c1 + p.d2*c12); 
    M13 = p.m3*p.L0*p.d2*c12;
    M21 = M12; M31 = M13;
    
    M = [M11, M12, M13;
         M21, M22, M23;
         M31, M32, M33];
     
    % --- CORIOLIS C(q,dq) ---
    C = compute_coriolis_num(M, y, p);
    
    % --- GRAVITACE G(q) ---
    G = zeros(3,1);
    G(1) = 0;
    G(2) = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G(3) = p.g * p.m3*p.d2*s12;
    
    % Vnější síly (zahrnují moment motoru na první souřadnici)
    Q = [Tau; 0; 0] - [p.b_phi*dq1; p.b_theta1*dq2; p.b_theta2*dq3] ...
        - [0; p.c_theta1*tanh(dq2/p.v_eps); p.c_theta2*tanh(dq3/p.v_eps)];
    
    % Řešení pro zrychlení
    ddq = M \ (Q - C*dq - G);
    dydt = [dq; ddq];
end

function C = compute_coriolis_num(~, y, p)
    % Rozbalení pro výpočet C (pozor na indexy, y obsahuje i pozice q)
    q2 = y(2); q3 = y(3);
    dq1 = y(4); dq2 = y(5); dq3 = y(6);
    
    s1 = sin(q2); c1 = cos(q2);
    s2 = sin(q3); 
    s12 = sin(q2+q3); c12 = cos(q2+q3);
    
    H = p.L1*s1 + p.d2*s12;
    dH_dq2 = p.L1*c1 + p.d2*c12;
    dH_dq3 = p.d2*c12;
    
    dM11_m2_dq2 = p.m2 * p.d1^2 * sin(2*q2);
    dM11_m3_dq2 = 2 * p.m3 * H * dH_dq2;
    dM11_m3_dq3 = 2 * p.m3 * H * dH_dq3;
    
    dM11_rot_dq2 = (p.I2_yy - p.I2_xx) * sin(2*q2) + ...
                   (p.I3_yy - p.I3_xx) * sin(2*(q2+q3));
    dM11_rot_dq3 = (p.I3_yy - p.I3_xx) * sin(2*(q2+q3));
    
    lambda1 = dM11_m2_dq2 + dM11_m3_dq2 + dM11_rot_dq2;
    lambda2 = dM11_m3_dq3 + dM11_rot_dq3;
    h = -p.m3 * p.L1 * p.d2 * s2; 
    
    C = zeros(3,3);
    C(1,1) = 0.5 * (lambda1 * dq2 + lambda2 * dq3);
    % + členy z derivací M12(q2,q3) a M13(q2+q3) (dostředivé síly na rameno)
    C(1,2) = 0.5 * lambda1 * dq1 ...
             - p.L0 * ((p.m2*p.d1 + p.m3*p.L1)*s1*dq2 + p.m3*p.d2*s12*(dq2 + dq3));
    C(1,3) = 0.5 * lambda2 * dq1 - p.L0 * p.m3 * p.d2 * s12 * (dq2 + dq3);
    C(2,1) = -0.5 * lambda1 * dq1;
    C(2,2) = h * dq3;
    C(2,3) = h * (dq2 + dq3);
    C(3,1) = -0.5 * lambda2 * dq1;
    C(3,2) = -h * dq2;
    C(3,3) = 0;
end

