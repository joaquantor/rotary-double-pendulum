%% Porovnání Simulace vs. Experiment (Validace modelu)
% Verze: Legendy mimo graf + nastavitelný časový úsek
clear; clc; close all;

%% 0. UŽIVATELSKÉ NASTAVENÍ
% Zde nastavte časový interval pro grafy [od do] v sekundách
t_plot = [0 37]; 

%% 1. NAČTENÍ DAT A PŘEPOČET ENKODERŮ
try
    load('experimenty.mat');
    
    if ~exist('kyvadlo2', 'var')
        error('Proměnná "kyvadlo2" nebyla nalezena v experimenty.mat');
    end
    
    % --- KONFIGURACE ENKODERŮ ---
    % Kyvadla: 1024 pulzů/ot (PPR) -> 4096 hran/ot
    % Motor:   500 pulzů/ot (PPR)  -> 2000 hran/ot
    ENC_PPR_PEND = 1024;  
    ENC_PPR_MOT  = 500;   
    GEAR_RATIO   = 128 / 14; 
    
    k_phi    = (2*pi) / (4 * ENC_PPR_MOT * GEAR_RATIO); 
    k_theta  = (2*pi) / (4 * ENC_PPR_PEND);             
    
    % --- ZPRACOVÁNÍ SIGNÁLŮ ---
    raw_vals = kyvadlo2.signals.values;
    raw_time = kyvadlo2.time;
    
    exp.t = raw_time - raw_time(1); % Čas od 0
    exp.u = raw_vals(:, 1);         % Napětí [V]
    
    exp.phi    = raw_vals(:, 2) * k_phi;
    exp.theta1 = raw_vals(:, 3) * k_theta;
    exp.theta2 = raw_vals(:, 4) * k_theta;

    fprintf('Data načtena. Celková délka záznamu: %.2f s\n', exp.t(end));

catch ME
    error(['Chyba při zpracování dat: ' ME.message]);
end

%% 2. PARAMETRY SYSTÉMU
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;       
p.m2 = 0.0356 + 0.0059;  
p.m3 = 0.0435;  
p.L0 = 0.201;   
p.L1 = 0.163;   
p.L2 = 0.150;   
p.d_ram = -0.059;  
p.d1 = 0.078;    
p.d2 = 0.075;    
p.gear_ratio = GEAR_RATIO; 
p.I1_zz = 0.0051; 
p.I2_xx = 0.00005; p.I2_yy = 0.000318; p.I2_zz = 0.000318;    
p.I3_xx = 0.00005; p.I3_yy = 0.000326; p.I3_zz = 0.000326;    
p.g = 9.81;     
p.U_pwm_amp = 14.4;  
p.U_ctrl_max = 10.0; 
p.R = 2.20;          
p.Kt = 0.052;        
p.Ke = 0.0516; 
p.Ke_out = p.Ke * p.gear_ratio; 
p.Kt_out = p.Kt * p.gear_ratio;
p.b_phi = 0.135;      
p.b_theta1 = 0.0003;  
p.b_theta2 = 0.0003;  

%% 3. SIMULACE
y0 = [exp.phi(1); exp.theta1(1); exp.theta2(1); 0; 0; 0]; 

% Oříznutí simulace podle t_plot (pokud chceme ušetřit čas výpočtu)
sim_end_time = min(t_plot(2), exp.t(end));
if sim_end_time < t_plot(1), sim_end_time = exp.t(end); end

t_span = [exp.t(1) exp.t(end)]; 
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);

fprintf('Spouštím simulaci...\n');
[t_sim, y_sim] = ode45(@(t, y) dynamics_validated(t, y, p, exp.t, exp.u), t_span, y0, options);
fprintf('Simulace dokončena.\n');

%% 4. VIZUALIZACE S LEGENDAMI MIMO
f = figure('Name', 'Validace Modelu', 'Color', 'w');

% 1. Napětí
ax1 = subplot(4,1,1);
plot(exp.t, exp.u, 'b', 'LineWidth', 1.0); hold on;
u_sim_check = interp1(exp.t, exp.u, t_sim, 'nearest');
plot(t_sim, u_sim_check, 'r:', 'LineWidth', 1.5);
ylabel('U [V]'); title('Vstupní napětí'); 
grid on;
% Legenda mimo graf
legend({'Experiment', 'Model Input'}, 'Location', 'northeastoutside');

% 2. Rameno
ax2 = subplot(4,1,2);
plot(exp.t, exp.phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,1), 'r--', 'LineWidth', 1.5);
ylabel('\phi [rad]'); title('Rameno'); 
grid on;
legend({'Experiment', 'Simulace'}, 'Location', 'northeastoutside');

% 3. Kyvadlo 1
ax3 = subplot(4,1,3);
plot(exp.t, exp.theta1, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,2), 'r--', 'LineWidth', 1.5);
ylabel('\theta_1 [rad]'); title('1. Kyvadlo (Vnitřní)'); 
grid on;
legend({'Experiment', 'Simulace'}, 'Location', 'northeastoutside');

% 4. Kyvadlo 2
ax4 = subplot(4,1,4);
plot(exp.t, exp.theta2, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,3), 'r--', 'LineWidth', 1.5);
ylabel('\theta_2 [rad]'); title('2. Kyvadlo (Vnější)'); 
xlabel('Čas [s]'); grid on;
legend({'Experiment', 'Simulace'}, 'Location', 'northeastoutside');

% Propojení os a nastavení limitů
linkaxes([ax1, ax2, ax3, ax4], 'x');
xlim(t_plot);

%% 5. FUNKCE DYNAMIKY
function dydt = dynamics_validated(t, y, p, t_exp, u_exp)
    q1 = y(1); q2 = y(2); q3 = y(3);       
    dq1 = y(4); dq2 = y(5); dq3 = y(6);    
    dq = [dq1; dq2; dq3];
    
    u_in_raw = interp1(t_exp, u_exp, t, 'previous', 'extrap');
    u_in = max(min(u_in_raw, p.U_ctrl_max), -p.U_ctrl_max);
    
    U_app = (u_in / p.U_ctrl_max) * p.U_pwm_amp;
    Current = (U_app - p.Ke_out * dq1) / p.R;
    Tau = p.Kt_out * Current;
    
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
    
    M = [M11, M12, M13; M21, M22, M23; M31, M32, M33];
     
    C = compute_coriolis_num(M, y, p);
    
    G = zeros(3,1);
    G(2) = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G(3) = p.g * p.m3*p.d2*s12;
    
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