%% Validace: Dvojité kyvadlo na fixní základně (kyvadlo2fixB)
% CÍL: Naladit tření b_theta1 a b_theta2 na datech volného kmitání.
% VSTUP: Naměřený pohyb základny (měl by být cca 0)
% VÝSTUP: Pohyb obou kyvadel (Theta1, Theta2)
clear; clc; close all;

%% 0. UŽIVATELSKÉ NASTAVENÍ (LADĚNÍ)
% Zde laďte parametry tření pro obě kyvadla:
p.b_theta1 = 0.00005;  % Tření v 1. kloubu (vnitřní)
p.b_theta2 = 0.00005;  % Tření v 2. kloubu (vnější)

t_plot = [0 20];       % Délka zobrazení (dvojkyvadlo se uklidňuje déle)
smooth_window = 40;    % Filtrace derivací vstupu

%% 1. NAČTENÍ DAT
try
    load('experimenty.mat');
    if ~exist('kyvadlo2fixB', 'var'), error('Chybí data kyvadlo2fixB'); end
    data = kyvadlo2fixB;
    
    % Enkodéry
    ENC_PPR_PEND = 1024; ENC_PPR_MOT = 500; GEAR = 128/14;
    k_phi = (2*pi)/(4*ENC_PPR_MOT*GEAR);
    k_theta = (2*pi)/(4*ENC_PPR_PEND);
    
    vals = data.signals.values;
    time = data.time;
    exp.t = time(:) - time(1);
    
    % Dle experviz.m: col2=Base, col3=Pend1, col4=Pend2
    exp.phi    = vals(:,2) * k_phi;
    exp.theta1 = vals(:,3) * k_theta;
    exp.theta2 = vals(:,4) * k_theta;

    % Kinematika základny (prochvění "fixní" základny může budit kyvadla)
    dt = mean(diff(exp.t));
    raw_dphi = gradient(exp.phi, dt);
    exp.dphi = smoothdata(raw_dphi, 'gaussian', smooth_window);
    raw_ddphi = gradient(exp.dphi, dt);
    exp.ddphi = smoothdata(raw_ddphi, 'gaussian', smooth_window);
    
    fprintf('Data dvojitého kyvadla načtena.\n');

catch ME
    error(['Chyba: ' ME.message]);
end

%% 2. PARAMETRY SYSTÉMU (Plná sestava)
% Rameno
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057; 
p.d_ram = -0.059;   
p.I1_zz = 0.0051; 

% 1. Kyvadlo
p.m2 = 0.0356 + 0.0059;  
p.d1 = 0.078;    
p.L0 = 0.201; 
p.L1 = 0.163; 
p.I2_xx = 0.00005; p.I2_yy = 0.000318; p.I2_zz = 0.000318;

% 2. Kyvadlo (ZDE MUSÍ BÝT AKTIVNÍ)
p.m3 = 0.0435;       % Hmotnost 2. kyvadla
p.d2 = 0.073;        % Těžiště (dle PDF 0.073, kód měl 0.075)
p.L2 = 0.150;        % Délka
p.I3_xx = 0.00005;   
p.I3_yy = 0.000326; 
p.I3_zz = 0.000326;

p.g = 9.81;

%% 3. SIMULACE (Sub-dynamika pro Theta1, Theta2)
% Stav: y = [theta1; theta2; dtheta1; dtheta2]
y0 = [exp.theta1(1); exp.theta2(1); 0; 0]; 

t_span = [exp.t(1) exp.t(end)];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);

[t_sim, y_sim] = ode45(@(t,y) dynamics_double_pendulum(t, y, p, exp), t_span, y0, options);

%% 4. VIZUALIZACE
figure('Name', 'Ladění Dvojitého Kyvadla (fixB)', 'Color', 'w');

% 1. Pohyb základny (Kontrola, zda je fixní)
subplot(3,1,1);
plot(exp.t, exp.phi, 'k');
ylabel('\phi [rad]'); title('Základna (Měla by být fixní)');
grid on; xlim(t_plot);

% 2. Kyvadlo 1
subplot(3,1,2);
plot(exp.t, exp.theta1, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,1), 'r--', 'LineWidth', 1.5);
ylabel('\theta_1 [rad]'); title('1. Kyvadlo (Vnitřní)');
legend('Experiment', 'Simulace'); grid on; xlim(t_plot);

% 3. Kyvadlo 2
subplot(3,1,3);
plot(exp.t, exp.theta2, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,2), 'r--', 'LineWidth', 1.5);
ylabel('\theta_2 [rad]'); title(sprintf('2. Kyvadlo (b_{\\theta2}=%g)', p.b_theta2));
xlabel('Čas [s]');
legend('Experiment', 'Simulace'); grid on; xlim(t_plot);

%% 5. DYNAMIKA (Kinematické buzení Dvojkyvadla)
function dydt = dynamics_double_pendulum(t, y, p, data)
    % y = [theta1; theta2; dtheta1; dtheta2]
    q2 = y(1); q3 = y(2);       % theta1, theta2
    dq2 = y(3); dq3 = y(4);     % dtheta1, dtheta2
    
    % Interpolace vstupu (Základna)
    q1_in  = interp1(data.t, data.phi,   t, 'linear', 0); % phi
    dq1_in = interp1(data.t, data.dphi,  t, 'linear', 0); % dphi
    ddq1_in= interp1(data.t, data.ddphi, t, 'linear', 0); % ddphi
    
    % Sestavení stavového vektoru pro výpočet matic (plný vektor 3x1)
    % q_full = [phi; theta1; theta2]
    q_all = [q1_in; q2; q3];
    dq_all = [dq1_in; dq2; dq3];
    
    % --- MATICE M(q) ---
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2+q3); s12 = sin(q2+q3);
    
    M22 = p.m2*p.d1^2 + p.m3*(p.L1^2 + p.d2^2 + 2*p.L1*p.d2*c2) + p.I2_zz + p.I3_zz;
    M33 = p.m3*p.d2^2 + p.I3_zz;
    M23 = p.m3*p.d2^2 + p.m3*p.L1*p.d2*c2 + p.I3_zz;
    M32 = M23;
    
    % Vazební členy od základny
    M21 = p.m2*p.L0*p.d1*c1 + p.m3*p.L0*(p.L1*c1 + p.d2*c12); 
    M31 = p.m3*p.L0*p.d2*c12;
    
    % Submatice pro kyvadla
    M_sub = [M22, M23; 
             M32, M33];
         
    M_coupling = [M21; M31];
    
    % --- CORIOLIS C(q,dq) ---
    % Potřebujeme spočítat Coriolisovy síly. Využijeme plnou funkci, 
    % ale vybereme jen řádky 2 a 3.
    
    % Pomocný vektor pro full dynamics function [q2, q3, dq1, dq2, dq3]
    y_full_calc = [0; q2; q3; dq1_in; dq2; dq3]; % 0 na pozici q1(phi) není v C důležité
    C_full = compute_coriolis_num_internal(y_full_calc, p);
    
    % Coriolisův vektor (C * dq)
    C_vec_full = C_full * dq_all;
    C_vec_sub = C_vec_full(2:3);
    
    % --- GRAVITACE G(q) ---
    G2 = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G3 = p.g * p.m3*p.d2*s12;
    G_sub = [G2; G3];
    
    % --- TŘENÍ ---
    Friction_sub = [p.b_theta1 * dq2; 
                    p.b_theta2 * dq3];
                
    % --- ŘEŠENÍ SOUSTAVY ---
    % M_sub * ddq_sub + C_vec_sub + G_sub + M_coupling * ddq1_in = -Friction
    % M_sub * ddq_sub = -Friction - C - G - M_couple * ddq1
    
    RHS = -Friction_sub - C_vec_sub - G_sub - (M_coupling * ddq1_in);
    
    ddq_sub = M_sub \ RHS;
    
    dydt = [dq2; dq3; ddq_sub];
end

% Kopie funkce pro Coriolis z hlavního skriptu (nutná pro běh)
function C = compute_coriolis_num_internal(y_full, p)
    % y_full = [unused; q2; q3; dq1; dq2; dq3]
    q2 = y_full(2); q3 = y_full(3);
    dq1 = y_full(4); dq2 = y_full(5); dq3 = y_full(6);
    
    s1 = sin(q2); c1 = cos(q2); s2 = sin(q3); s12 = sin(q2+q3); c12 = cos(q2+q3);
    H = p.L1*s1 + p.d2*s12;
    dH_dq2 = p.L1*c1 + p.d2*c12; dH_dq3 = p.d2*c12;
    
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