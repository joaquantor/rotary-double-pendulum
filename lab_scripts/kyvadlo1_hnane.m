%% Validace: Kyvadlo 1 Hnané (Motor + Rameno + 1. Kyvadlo)
% Vstup: Naměřené napětí (U)
% Výstup: Porovnání pohybu ramene (phi) a kyvadla (theta1)
clear; clc; close all;

%% 0. UŽIVATELSKÉ NASTAVENÍ (LADĚNÍ)
% Parametry, které nejvíce ovlivňují shodu s realitou:
p.b_phi = 0.135;       % Tření v rameni (ovlivňuje, jak moc se rameno roztočí)
p.b_theta1 = 0.00005;  % Tření v kyvadle (ovlivňuje útlum kmitů kyvadla)

% Interval vykreslení
t_plot = [0 5]; 

%% 1. NAČTENÍ DAT (kyvadlo1)
try
    load('experimenty.mat');
    
    if ~exist('kyvadlo1', 'var')
        error('Proměnná "kyvadlo1" chybí. Zkontrolujte experimenty.mat.');
    end
    
    data_struct = kyvadlo1;
    fprintf('Načtena struktura "kyvadlo1".\n');
    
    % --- KONFIGURACE ENKODERŮ ---
    ENC_PPR_PEND = 1024; ENC_PPR_MOT = 500; GEAR_RATIO = 128/14;
    k_phi = (2*pi)/(4*ENC_PPR_MOT*GEAR_RATIO);
    k_theta = (2*pi)/(4*ENC_PPR_PEND);
    
    % --- ZPRACOVÁNÍ SIGNÁLŮ ---
    % Struktura kyvadlo1 dle experviz.m: [Napětí, Rameno, Kyvadlo1]
    raw_vals = data_struct.signals.values;
    raw_time = data_struct.time;
    
    exp.t = raw_time(:) - raw_time(1);
    
    % Sloupec 1: Napětí
    exp.u = raw_vals(:, 1); 
    
    % Sloupec 2: Rameno (Phi)
    exp.phi = raw_vals(:, 2) * k_phi;
    
    % Sloupec 3: Kyvadlo 1 (Theta1)
    exp.theta1 = raw_vals(:, 3) * k_theta;
    
    % 2. kyvadlo zde není
    exp.theta2 = zeros(size(exp.t));

    fprintf('Data připravena. Max napětí: %.2f V\n', max(abs(exp.u)));

catch ME
    error(['Chyba dat: ' ME.message]);
end

%% 2. PARAMETRY SYSTÉMU
% Motor a převodovka
p.U_pwm_amp = 14.4; p.U_ctrl_max = 10.0; 
p.R = 2.20; p.Kt = 0.052; p.Ke = 0.0516; 
p.gear_ratio = GEAR_RATIO;
p.Ke_out = p.Ke * p.gear_ratio; p.Kt_out = p.Kt * p.gear_ratio;

% Rameno (Mechanika)
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057; 
p.d_ram = -0.059;   
p.I1_zz = 0.0051; 
% Tření ramene (načteno z horní sekce)
fprintf('Tření ramene b_phi = %g\n', p.b_phi);

% 1. Kyvadlo
p.m2 = 0.0556 + 0.0059;  
p.d1 = 0.078;    
p.L0 = 0.201; 
p.L1 = 0.163; 
p.I2_xx = 0.00005; p.I2_yy = 0.000318; p.I2_zz = 0.000318;
% Tření kyvadla (načteno z horní sekce)
fprintf('Tření kyvadla b_theta1 = %g\n', p.b_theta1);

% 2. Kyvadlo (Virtuálně odstraněno - epsilon)
eps_m = 1e-7;
p.m3 = eps_m; p.d2 = 0; p.L2 = 0; 
p.I3_xx = eps_m; p.I3_yy = eps_m; p.I3_zz = eps_m;
p.b_theta2 = 0;

p.g = 9.81;

%% 3. SIMULACE (Full Dynamics)
% Počáteční podmínky [phi, theta1, theta2, dphi, dtheta1, dtheta2]
y0 = [exp.phi(1); exp.theta1(1); 0; 0; 0; 0]; 

t_span = [exp.t(1) exp.t(end)];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);

fprintf('Spouštím simulaci (Dynamics Validated)...\n');
[t_sim, y_sim] = ode45(@(t, y) dynamics_validated(t, y, p, exp.t, exp.u), t_span, y0, options);

%% 4. VIZUALIZACE
figure('Name', 'Validace: Kyvadlo 1 Hnané', 'Color', 'w');

% 1. Vstupní napětí
ax1 = subplot(3,1,1);
plot(exp.t, exp.u, 'b', 'LineWidth', 1); hold on;
plot(t_sim, interp1(exp.t, exp.u, t_sim, 'nearest'), 'r:');
ylabel('Napětí [V]'); title('Vstupní buzení'); 
legend('Experiment', 'Simulace Input', 'Location', 'northeastoutside'); grid on;

% 2. Rameno (Phi) - Zde vidíme kvalitu modelu motoru a tření
ax2 = subplot(3,1,2);
plot(exp.t, exp.phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,1), 'r--', 'LineWidth', 1.5);
ylabel('\phi [rad]'); title('Pohyb Ramene'); 
legend('Experiment', 'Simulace', 'Location', 'northeastoutside'); grid on;

% 3. Kyvadlo (Theta1) - Zde vidíme dynamiku vazby
ax3 = subplot(3,1,3);
plot(exp.t, exp.theta1, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,2), 'r--', 'LineWidth', 1.5);
ylabel('\theta_1 [rad]'); xlabel('Čas [s]');
title('Pohyb 1. Kyvadla'); 
legend('Experiment', 'Simulace', 'Location', 'northeastoutside'); grid on;

linkaxes([ax1, ax2, ax3], 'x');
xlim(t_plot);

%% 5. FUNKCE DYNAMIKY (Kompletní model s motorem)
function dydt = dynamics_validated(t, y, p, t_exp, u_exp)
    q1 = y(1); q2 = y(2); q3 = y(3);       
    dq1 = y(4); dq2 = y(5); dq3 = y(6);    
    dq = [dq1; dq2; dq3];
    
    % Interpolace vstupu (Napětí)
    u_in_raw = interp1(t_exp, u_exp, t, 'previous', 'extrap');
    u_in = max(min(u_in_raw, p.U_ctrl_max), -p.U_ctrl_max);
    
    % Model DC Motoru
    U_app = (u_in / p.U_ctrl_max) * p.U_pwm_amp;
    Current = (U_app - p.Ke_out * dq1) / p.R;
    Tau = p.Kt_out * Current;
    
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2+q3); s12 = sin(q2+q3);
    
    % Matice M
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
    
    M = [M11, M12, M13; M21, M22, M23; M31, M32, M33];
     
    C = compute_coriolis_num(M, y, p);
    
    G = zeros(3,1);
    G(2) = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G(3) = p.g * p.m3*p.d2*s12;
    
    % Vnější síly
    Q = [Tau; 0; 0] - [p.b_phi*dq1; p.b_theta1*dq2; p.b_theta2*dq3];
    
    ddq = M \ (Q - C*dq - G);
    dydt = [dq; ddq];
end

function C = compute_coriolis_num(~, y, p)
    q2 = y(2); q3 = y(3);
    dq1 = y(4); dq2 = y(5); dq3 = y(6);
    s1 = sin(q2); c1 = cos(q2); s2 = sin(q3); s12 = sin(q2+q3); c12 = cos(q2+q3);
    
    H = p.L1*s1 + p.d2*s12;
    dH_dq2 = p.L1*c1 + p.d2*c12;
    dH_dq3 = p.d2*c12;
    
    dM11_m2_dq2 = p.m2 * p.d1^2 * sin(2*q2);
    dM11_m3_dq2 = 2 * p.m3 * H * dH_dq2;
    dM11_m3_dq3 = 2 * p.m3 * H * dH_dq3;
    
    dM11_rot_dq2 = (p.I2_yy - p.I2_xx) * sin(2*q2) + (p.I3_yy - p.I3_xx) * sin(2*(q2+q3));
    dM11_rot_dq3 = (p.I3_yy - p.I3_xx) * sin(2*(q2+q3));
    
    lambda1 = dM11_m2_dq2 + dM11_m3_dq2 + dM11_rot_dq2;
    lambda2 = dM11_m3_dq3 + dM11_rot_dq3;
    h = -p.m3 * p.L1 * p.d2 * s2; 
    
    C = zeros(3,3);
    C(1,1) = 0.5 * (lambda1 * dq2 + lambda2 * dq3);
    C(1,2) = 0.5 * lambda1 * dq1;
    C(1,3) = 0.5 * lambda2 * dq1;
    C(2,1) = -0.5 * lambda1 * dq1;
    C(2,2) = h * dq3;
    C(2,3) = h * (dq2 + dq3);
    C(3,1) = -0.5 * lambda2 * dq1;
    C(3,2) = -h * dq2;
    C(3,3) = 0;
end