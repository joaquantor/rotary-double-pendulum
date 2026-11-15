%% Testování Ramene Bez Zátěže (Identifikace motoru)
clear; clc; close all;

%% 1. NAČTENÍ DAT "BEZ ZÁTĚŽE"
try
    load('experimenty.mat');
    
    % Hledáme strukturu 'bezzateze'
    if exist('sezavazim', 'var')
        raw_struct = bezzateze;
        fprintf('Načtena struktura "sezavazim".\n');
    else
        error('Proměnná "sezavazim" nebyla v souboru experimenty.mat nalezena.');
    end
    
    % --- KONFIGURACE ENKODERU (Pouze motor) ---
    ENC_PPR_MOT  = 500;        % HEDS-55xx na motoru
    GEAR_RATIO   = 128 / 14;   % Převodovka
    
    % Přepočet tiků na radiány (rameno za převodovkou)
    k_phi = (2*pi) / (4 * ENC_PPR_MOT * GEAR_RATIO); 
    
    % Zpracování signálů
    raw_vals = raw_struct.signals.values;
    raw_time = raw_struct.time;
    
    % --- OPRAVA: Vynucení sloupcových vektorů pomocí (:), aby nedošlo k chybě dimenzí ---
    exp.t = raw_time(:); 
    exp.t = exp.t - exp.t(1); % Čas od 0
    
    u_raw = raw_vals(:, 1);
    phi_raw = raw_vals(:, 2);
    
    exp.u = u_raw(:);          % Napětí [V]
    exp.phi = phi_raw(:) * k_phi; % Poloha ramene [rad]
    
    % Odhad rychlosti z dat (numerická derivace)
    % Nyní máme jistotu, že exp.phi i exp.t jsou sloupce
    diff_phi = diff(exp.phi);
    diff_t = diff(exp.t);
    
    % Ošetření dělení nulou (pokud by se v čase objevily duplicity)
    diff_t(diff_t == 0) = 1e-6; 
    
    velocity_calc = diff_phi ./ diff_t;
    
    % Doplnění nuly na začátek, aby seděl počet prvků
    exp.dphi = [0; velocity_calc]; 
    
    % Filtrace rychlosti (hrubá derivace šumí)
    exp.dphi = smoothdata(exp.dphi, 'gaussian', 20);

    fprintf('Data připravena. Délka: %.2f s\n', exp.t(end));

catch ME
    error(['Chyba ve zpracování dat: ' ME.message]);
end

%% 2. PARAMETRY SYSTÉMU (Pouze Rameno)
% Parametry motoru a elektroniky
p.U_pwm_amp = 14.4;  
p.U_ctrl_max = 10.0; 
p.R = 2.20;          
p.Kt = 0.052;        
p.Ke = 0.0516; 
p.gear_ratio = GEAR_RATIO;
p.Ke_out = p.Ke * p.gear_ratio; 
p.Kt_out = p.Kt * p.gear_ratio;

% --- MECHANICKÉ PARAMETRY (Rameno) ---
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057; 
p.d_ram = -0.059;   
p.I1_zz = 0.0051;   
p.b_phi = 0.135;     

p.J_total = p.I1_zz; 

%% 3. SIMULACE (1 DOF Model)
y0 = [exp.phi(1); 0]; % [Poloha; Rychlost]
t_span = [exp.t(1) exp.t(end)];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);

fprintf('Spouštím simulaci ramene bez zátěže...\n');
[t_sim, y_sim] = ode45(@(t, y) dynamics_arm_only(t, y, p, exp.t, exp.u), t_span, y0, options);

%% 4. VIZUALIZACE
figure('Name', 'Test Ramene Bez Zátěže', 'Color', 'w');

% 1. Napětí
ax1 = subplot(3,1,1);
plot(exp.t, exp.u, 'b'); hold on;
% Interpolace pro vykreslení simulovaného vstupu (kontrola)
u_sim_plot = interp1(exp.t, exp.u, t_sim, 'previous');
plot(t_sim, u_sim_plot, 'r:');
ylabel('Napětí [V]'); title('Vstupní signál');
legend({'Měření', 'Simulace'}, 'Location', 'northeastoutside'); grid on;

% 2. Poloha
ax2 = subplot(3,1,2);
plot(exp.t, exp.phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim(:,1), 'r--', 'LineWidth', 1.5);
ylabel('\phi [rad]'); title('Poloha Ramene');
legend({'Měření', 'Simulace'}, 'Location', 'northeastoutside'); grid on;

% 3. Rychlost
ax3 = subplot(3,1,3);
plot(exp.t, exp.dphi, 'b', 'LineWidth', 1); hold on;
plot(t_sim, y_sim(:,2), 'r', 'LineWidth', 1.5);
ylabel('\omega [rad/s]'); title('Rychlost Ramene');
xlabel('Čas [s]');
legend({'Měření (deriv)', 'Simulace'}, 'Location', 'northeastoutside'); grid on;

linkaxes([ax1, ax2, ax3], 'x');
xlim([0 exp.t(end)]);

%% 5. DYNAMIKA (Pouze Rameno)
function dydt = dynamics_arm_only(t, y, p, t_exp, u_exp)
    % y(1) = phi
    % y(2) = dphi
    dphi = y(2);
    
    % Vstup (Napětí)
    u_in_raw = interp1(t_exp, u_exp, t, 'previous', 'extrap');
    u_in = max(min(u_in_raw, p.U_ctrl_max), -p.U_ctrl_max);
    
    % Rovnice motoru
    U_app = (u_in / p.U_ctrl_max) * p.U_pwm_amp;
    Current = (U_app - p.Ke_out * dphi) / p.R;
    Tau_motor = p.Kt_out * Current;
    
    % Tření
    Tau_friction = p.b_phi * dphi;
    
    % Zrychlení: J * alpha = Tau_motor - Tření
    ddphi = (Tau_motor - Tau_friction) / p.J_total;
    
    dydt = [dphi; ddphi];
end