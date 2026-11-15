%% Ladění 1. Kyvadla: Kinematické buzení
% VSTUP: Naměřený pohyb ramene (phi, dphi, ddphi)
% VÝSTUP: Odezva kyvadla (theta1)
% CÍL: Naladit b_theta1 tak, aby simulace seděla na experiment.
clear; clc; close all;

%% 0. UŽIVATELSKÉ NASTAVENÍ
% Zde měňte tření, dokud se simulace (červená) nepřekryje s daty (modrá)
p.b_theta1 = 0.00005;  

% Parametr vyhlazení derivace (čím vyšší, tím hladší průběh zrychlení ramene)
smooth_window = 40; 
t_plot = [0 15]; % Časový úsek

%% 1. NAČTENÍ DAT A PŘÍPRAVA KINEMATIKY
try
    load('experimenty.mat');
    if ~exist('kyvadlo1vol', 'var'), error('Chybí data kyvadlo1vol'); end
    data = kyvadlo1vol;

    % Konfigurace enkodérů
    ENC_PPR_PEND = 1024; ENC_PPR_MOT = 500; GEAR = 128/14;
    k_phi = (2*pi)/(4*ENC_PPR_MOT*GEAR);
    k_theta = (2*pi)/(4*ENC_PPR_PEND);

    % Extrakce
    vals = data.signals.values;
    time = data.time;
    
    exp.t = time(:) - time(1);
    exp.phi = vals(:,2) * k_phi;      % Poloha ramene (VSTUP)
    exp.theta1 = vals(:,3) * k_theta; % Poloha kyvadla (PRO POROVNÁNÍ)

    % --- VÝPOČET DERIVACÍ RAMENE (NUMERICKY) ---
    % Nutné pro výpočet setrvačných sil působících na kyvadlo
    dt = mean(diff(exp.t));
    
    % Rychlost ramene
    raw_dphi = gradient(exp.phi, dt);
    exp.dphi = smoothdata(raw_dphi, 'gaussian', smooth_window);
    
    % Zrychlení ramene (Klíčové!)
    raw_ddphi = gradient(exp.dphi, dt);
    exp.ddphi = smoothdata(raw_ddphi, 'gaussian', smooth_window);
    
    fprintf('Kinematika ramene vypočtena.\n');

catch ME
    error(['Chyba: ' ME.message]);
end

%% 2. PARAMETRY SYSTÉMU
% Parametry 2. kyvadla nulujeme, aby nepletly rovnice
p.m3 = 0; p.d2 = 0; p.L2 = 0; p.I3_xx=0; p.I3_yy=0; p.I3_zz=0; 
p.b_theta2 = 0;

% Parametry 1. kyvadla a ramene
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057; 
p.d_ram = -0.059;   
p.I1_zz = 0.0051; 

p.m2 = 0.0356 + 0.0059;  
p.d1 = 0.078;    
p.L0 = 0.201; 
p.L1 = 0.163; 
p.I2_xx = 0.00005; p.I2_yy = 0.000318; p.I2_zz = 0.000318;
p.g = 9.81;

%% 3. SIMULACE (Pouze rovnice pro Theta1)
% Stavový vektor: y = [theta1; dtheta1]
y0 = [exp.theta1(1); 0]; 

t_span = [exp.t(1) exp.t(end)];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);

% Voláme speciální funkci dynamiky s předepsaným pohybem
[t_sim, y_sim] = ode45(@(t,y) dynamics_kinematic(t, y, p, exp), t_span, y0, options);

%% 4. VIZUALIZACE
figure('Name', 'Ladění b_theta1 (Kinematické)', 'Color', 'w');

% Graf 1: Vstupní pohyb ramene (pro kontrolu)
subplot(3,1,1);
plot(exp.t, exp.phi, 'k', 'LineWidth', 1); hold on;
yyaxis right;
plot(exp.t, exp.ddphi, 'm:', 'LineWidth', 0.5);
ylabel('Zrychlení ramene [rad/s^2]');
yyaxis left;
ylabel('Poloha ramene [rad]');
title('Vstup simulace (Naměřený pohyb ramene)');
legend('Poloha \phi', 'Zrychlení \epsilon'); grid on; xlim(t_plot);

% Graf 2: Výsledek - Kyvadlo
subplot(3,1, [2 3]);
plot(exp.t, exp.theta1, 'b', 'LineWidth', 2); hold on;
plot(t_sim, y_sim(:,1), 'r', 'LineWidth', 1.5);
ylabel('\theta_1 [rad]'); xlabel('Čas [s]');
title(sprintf('VALIDACE: Odezva kyvadla (b_{\\theta1} = %g)', p.b_theta1));
legend('Experiment', 'Simulace (Model)', 'Location', 'northeast'); grid on;
xlim(t_plot);

%% 5. DYNAMIKA (Kinematické buzení)
function dydt = dynamics_kinematic(t, y, p, data)
    % y(1) = theta1, y(2) = dtheta1
    theta1 = y(1);
    dtheta1 = y(2);
    
    % 1. Zjistíme, co v čase 't' dělá rameno (Interpolace dat)
    % Linear interp je rychlá a pro hustá data dostatečná
    ddphi_in = interp1(data.t, data.ddphi, t, 'linear', 0);
    dphi_in  = interp1(data.t, data.dphi,  t, 'linear', 0);
    
    % 2. Sestavení rovnice pro Theta1
    % Rovnice z PDF (pro 2. řádek matice M):
    % M21*ddPhi + M22*ddTheta1 - 0.5*lam1*dPhi^2 + G2 = -b*dTheta1
    
    s1 = sin(theta1); c1 = cos(theta1);
    
    % Setrvačnost kyvadla (M22)
    M22 = p.m2 * p.d1^2 + p.I2_zz;
    
    % Vazební člen od zrychlení ramene (M21)
    M21 = p.m2 * p.L0 * p.d1 * c1;
    
    % Derivace M11 (pro odstředivou sílu) - Lambda1
    % lambda1 = d(M11)/dtheta1 = m2*d1^2*sin(2*theta1) + (Iyy-Ixx)*sin(2*theta1)
    dM11_dq2 = (p.m2 * p.d1^2 + p.I2_yy - p.I2_xx) * sin(2*theta1);
    
    % Gravitace
    G2 = p.g * p.m2 * p.d1 * s1;
    
    % 3. Výpočet úhlového zrychlení kyvadla
    % ddTheta1 = (-Tření - Gravitace - Setrvačná_od_ramene + Odstředivá) / M22
    
    Torque_friction = -p.b_theta1 * dtheta1;
    Torque_inertial = -M21 * ddphi_in;          % "Trhnutí" ramenem
    Torque_centrifugal = 0.5 * dM11_dq2 * dphi_in^2; % "Vynášení" odstředivou silou
    
    ddtheta1 = (Torque_friction - G2 + Torque_inertial + Torque_centrifugal) / M22;
    
    dydt = [dtheta1; ddtheta1];
end