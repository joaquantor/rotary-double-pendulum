%% Simulace Dvojitého Kyvadla buzeného naměřenými daty
clear; clc; close all;

%% 1. NAČTENÍ A PŘÍPRAVA DAT (Vstup pohybu ramene)
% Načtení dat (předpokládá se soubor z předchozího kroku)
data_table = readtable('kyvadlo2vol_rameno_data.csv');
t_meas = data_table.Time;
phi_meas = data_table.Arm_Angle_Rad;

% Filtrace a výpočet derivací (Rychlost a Zrychlení ramene)
% Používáme vyhlazení, protože numerická derivace šumu vytváří obrovské chyby
dt = mean(diff(t_meas));
windowSize = 20; % Velikost okna pro klouzavý průměr/filtr
b_fir = (1/windowSize)*ones(1,windowSize);
a_fir = 1;

% 1. Poloha (Phi)
phi_smooth = filtfilt(b_fir, a_fir, phi_meas);

% 2. Rychlost (dPhi) - centrální diference
dphi_raw = gradient(phi_smooth, dt);
dphi_smooth = filtfilt(b_fir, a_fir, dphi_raw);

% 3. Zrychlení (ddPhi)
ddphi_raw = gradient(dphi_smooth, dt);
ddphi_smooth = filtfilt(b_fir, a_fir, ddphi_raw);

% Vytvoření interpolantů pro použití v ode45 (aby funkce fungovala v libovolném čase t)
F_phi   = griddedInterpolant(t_meas, phi_smooth, 'linear', 'nearest');
F_dphi  = griddedInterpolant(t_meas, dphi_smooth, 'linear', 'nearest');
F_ddphi = griddedInterpolant(t_meas, ddphi_smooth, 'linear', 'nearest');

%% 2. PARAMETRY SYSTÉMU
% (Stejné jako ve vašem původním skriptu)
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;       
p.m2 = 0.0356 + 0.0059;  
p.m3 = 0.0435;  
p.L0 = 0.201;   
p.L1 = 0.163;   
p.L2 = 0.150;   
p.d_ram = 0;  
p.d1 = 0.081;    
p.d2 = 0.075;    
p.I1_zz = 0.0051; 
p.I2_xx = 0.00005; p.I2_yy = 0.000318; p.I2_zz = 0.000318;    
p.I3_xx = 0.00005; p.I3_yy = 0.000326; p.I3_zz = 0.000326;    
p.g = 9.81;     
% Tření
p.b_theta1 = 0.0007;  
p.b_theta2 = 0.0007;  

% Předání interpolantů do struktury parametrů
p.F_phi = F_phi;
p.F_dphi = F_dphi;
p.F_ddphi = F_ddphi;

%% 3. NASTAVENÍ A SPUŠTĚNÍ SIMULACE
% Stavy nyní jsou pouze: [theta1; theta2; dtheta1; dtheta2]
init_pos = [0; 0]; % Počáteční poloha kyvadel (svisle dolů)
init_vel = [0; 0]; % Počáteční rychlost kyvadel
y0 = [init_pos; init_vel];

t_span = [t_meas(1) t_meas(end)]; 
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

disp('Spouštím simulaci...');
[t, y] = ode45(@(t, y) dynamics_driven(t, y, p), t_span, y0, options);

%% 4. VIZUALIZACE A POROVNÁNÍ
% Rekonstrukce vstupu pro grafy
phi_sim = F_phi(t);

figure('Name', 'Výsledek simulace - Buzeno daty');
subplot(3,1,1);
plot(t_meas, phi_meas, 'Color', [0.7 0.7 0.7], 'LineWidth', 2); hold on;
plot(t, phi_sim, 'k--');
ylabel('\phi [rad]'); title('Vstup: Natočení ramene (Data vs Interpolace)');
legend('Naměřeno', 'Použito v simulaci'); grid on;

subplot(3,1,2);
plot(t, y(:,1)); ylabel('\theta_1 [rad]'); title('Výstup: 1. Kyvadlo'); grid on;

subplot(3,1,3);
plot(t, y(:,2)); ylabel('\theta_2 [rad]'); title('Výstup: 2. Kyvadlo'); grid on;

% Animace (volitelné)
% animate_system_driven(t, y, p);

%% 5. DYNAMIKA (Redukovaný model)
function dydt = dynamics_driven(t, y_red, p)
    % y_red obsahuje jen [theta1; theta2; dtheta1; dtheta2]
    q2 = y_red(1); q3 = y_red(2);
    dq2 = y_red(3); dq3 = y_red(4);
    
    % Získání vynuceného pohybu ramene v čase t
    q1 = p.F_phi(t);
    dq1 = p.F_dphi(t);
    ddq1 = p.F_ddphi(t); % Nutné pro setrvačné síly
    
    % Sestavíme "plný" stavový vektor, abychom mohli recyklovat vaši logiku matic
    % y_full = [phi; theta1; theta2; dphi; dtheta1; dtheta2]
    y_full = [q1; q2; q3; dq1; dq2; dq3];
    
    % --- MATICE M(q) ---
    % (Kopie vašeho kódu pro výpočet prvků M)
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2+q3); s12 = sin(q2+q3);
    
    % Prvky M, které potřebujeme (jen řádky 2 a 3)
    % M21, M22, M23
    M22 = p.m2*p.d1^2 + p.m3*(p.L1^2 + p.d2^2 + 2*p.L1*p.d2*c2) + p.I2_zz + p.I3_zz;
    M33 = p.m3*p.d2^2 + p.I3_zz;
    M23 = p.m3*p.d2^2 + p.m3*p.L1*p.d2*c2 + p.I3_zz;
    M32 = M23;
    
    M12 = p.m2*p.L0*p.d1*c1 + p.m3*p.L0*(p.L1*c1 + p.d2*c12); 
    M13 = p.m3*p.L0*p.d2*c12;
    M21 = M12;
    M31 = M13;
    
    % Sestavení sub-matic
    % Rovnice systému: 
    % [M11 M12 M13] [ddq1]   + C1... = tau
    % [M21 M22 M23] [ddq2]   + C2... = 0
    % [M31 M32 M33] [ddq3]   + C3... = 0
    
    % My řešíme spodní dvě řádky pro [ddq2; ddq3]:
    % [M22 M23] * [ddq2] = - [M21] * ddq1 - C_terms - G_terms
    % [M32 M33]   [ddq3]   - [M31]
    
    M_sub = [M22, M23; M32, M33];
    M_coupling = [M21; M31];
    
    % --- CORIOLIS C(q,dq) ---
    % Použijeme vaši funkci, která vrací plnou matici 3x3
    % Musíme jí ale předat "plný" vektor stavů y_full
    C_full = compute_coriolis_num_local(y_full, p);
    
    % Coriolisovy síly pro řádky 2 a 3
    dq_full = [dq1; dq2; dq3];
    C_vec_full = C_full * dq_full;
    C_sub = C_vec_full(2:3); % Jen 2. a 3. prvek
    
    % --- GRAVITACE G(q) ---
    G2 = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G3 = p.g * p.m3*p.d2*s12;
    G_sub = [G2; G3];
    
    % --- TŘENÍ ---
    Fric_sub = [p.b_theta1*dq2; p.b_theta2*dq3];
    
    % --- VÝPOČET ZRYCHLENÍ ---
    % M_sub * ddq_sub = - M_coupling * ddq1 - C_sub - G_sub - Fric
    
    RHS = -M_coupling * ddq1 - C_sub - G_sub - Fric_sub;
    
    ddq_sub = M_sub \ RHS;
    
    % Výstup: [dq2; dq3; ddq2; ddq3]
    dydt = [dq2; dq3; ddq_sub];
end

% Lokální kopie vaší funkce pro Coriolis (beze změny, jen přejmenovaná pro pořádek)
function C = compute_coriolis_num_local(y_full, p)
    q2 = y_full(2); q3 = y_full(3);
    dq1 = y_full(4); dq2 = y_full(5); dq3 = y_full(6);
    
    s1 = sin(q2); c1 = cos(q2);
    s2 = sin(q3); % c2 = cos(q3);
    s12 = sin(q2+q3); c12 = cos(q2+q3);
    
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
    h_val = -p.m3 * p.L1 * p.d2 * s2; 
    
    C = zeros(3,3);
    C(1,1) = 0.5 * (lambda1 * dq2 + lambda2 * dq3);
    C(1,2) = 0.5 * lambda1 * dq1;
    C(1,3) = 0.5 * lambda2 * dq1;
    C(2,1) = -0.5 * lambda1 * dq1;
    C(2,2) = h_val * dq3;
    C(2,3) = h_val * (dq2 + dq3);
    C(3,1) = -0.5 * lambda2 * dq1;
    C(3,2) = -h_val * dq2;
    C(3,3) = 0;
end