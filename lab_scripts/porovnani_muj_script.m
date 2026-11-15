clear; clc; close all;

%% 1. DEFINICE SYMBOLICKÝCH PROMĚNNÝCH
% Čas a gravitace
syms t g_acc real

% Geometrické parametry
syms L1 L2 L3 Lc1 Lc2 Lc3 d1 d2 real

% Hmotnostní parametry
syms m1 m2 m3 real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real
syms I3xx I3yy I3zz real

% Zobecněné souřadnice a derivace
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real

% Vektory stavu
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

%% 2. KINEMATIKA
disp('Počítám kinematiku...');

% Pomocné anonymní funkce
rotx = @(th) [1, 0, 0, 0; 0, cos(th), -sin(th), 0; 0, sin(th), cos(th), 0; 0, 0, 0, 1];
roty = @(th) [cos(th), 0, sin(th), 0; 0, 1, 0, 0; -sin(th), 0, cos(th), 0; 0, 0, 0, 1];
rotz = @(th) [cos(th), -sin(th), 0, 0; sin(th), cos(th), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
tx = @(d) [1 0 0 d; 0 1 0 0; 0 0 1 0; 0 0 0 1];
tz = @(d) [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];

% --- Transformace ---

% Těleso 1 (Rameno)
T01 = rotz(q1)*tx(L1)*tz(d1)*roty(-pi/2);

% Těžiště 1 (p1)
T_com1 = rotz(q1) * tx(-0.059) * tz(d1); % Dle tvého nastavení
p1_h = T_com1 * [0; 0; 0; 1];

% Těleso 2 (1. Kyvadlo)
T12 = rotz(q2)*tx(L2)*tz(d2);
T02 = T01 * T12; 

% Těžiště 2 (p2)
p2_h = T02 * [-Lc2; 0; -d2; 1];

% Těleso 3 (2. Kyvadlo)
T23 = rotz(q3)*tx(L3);
T03 = T02 * T23;

% Těžiště 3 (p3)
p3_h = T03 * [-Lc3; 0; 0; 1];

% --- Extrakce finálních vektorů ---
p1 = p1_h(1:3);
p2 = p2_h(1:3);
p3 = p3_h(1:3);

R1 = T01(1:3, 1:3); 
R2 = T02(1:3, 1:3); 
R3 = T03(1:3, 1:3); 

%% 3. VÝPOČET RYCHLOSTÍ
disp('Počítám Jakobiány...');

Jv1 = jacobian(p1, q); v1 = simplify(Jv1 * dq);
Jv2 = jacobian(p2, q); v2 = simplify(Jv2 * dq);
Jv3 = jacobian(p3, q); v3 = simplify(Jv3 * dq);

dR1 = diff(R1, q1)*dq1 + diff(R1, q2)*dq2 + diff(R1, q3)*dq3;
S1 = dR1 * R1.'; w1 = [S1(3,2); S1(1,3); S1(2,1)];

dR2 = diff(R2, q1)*dq1 + diff(R2, q2)*dq2 + diff(R2, q3)*dq3;
S2 = dR2 * R2.'; w2 = [S2(3,2); S2(1,3); S2(2,1)];

dR3 = diff(R3, q1)*dq1 + diff(R3, q2)*dq2 + diff(R3, q3)*dq3;
S3 = dR3 * R3.'; w3 = [S3(3,2); S3(1,3); S3(2,1)];

%% 4. ENERGIE A LAGRANGIÁN
disp('Sestavuji Lagrangián...');

I1_loc = diag([I1xx, I1yy, I1zz]);
I2_loc = diag([I2xx, I2yy, I2zz]);
I3_loc = diag([I3xx, I3yy, I3zz]);

T_trans = 0.5*m1*(v1.'*v1) + 0.5*m2*(v2.'*v2) + 0.5*m3*(v3.'*v3);

w1_loc = R1.' * w1;
w2_loc = R2.' * w2;
w3_loc = R3.' * w3;

T_rot = 0.5*(w1_loc.' * I1_loc * w1_loc) + ...
        0.5*(w2_loc.' * I2_loc * w2_loc) + ...
        0.5*(w3_loc.' * I3_loc * w3_loc);

T_total = simplify(T_trans + T_rot);
V_total = m1*g_acc*p1(3) + m2*g_acc*p2(3) + m3*g_acc*p3(3);

Lagrangian = simplify(T_total - V_total);

%% 5. DYNAMICKÉ MATICE (M, C, G)
% Poznámka: Pokud už máš funkce vygenerované, můžeš tuto sekci přeskočit pro zrychlení
if ~exist('get_M.m', 'file')
    disp('Odvozuji rovnice pohybu (M, C, G)...');
    M_mat = jacobian(jacobian(T_total, dq), dq);
    M_mat = simplify(M_mat);

    G_vec = jacobian(V_total, q).';
    G_vec = simplify(G_vec);

    n = 3;
    C_mat = sym(zeros(n,n));
    for k = 1:n
        for j = 1:n
            c_kj = 0;
            for i = 1:n
                term = 0.5 * (diff(M_mat(k,j), q(i)) + diff(M_mat(k,i), q(j)) - diff(M_mat(i,j), q(k)));
                c_kj = c_kj + term * dq(i);
            end
            C_mat(k,j) = c_kj;
        end
    end
    C_mat = simplify(C_mat);

    disp('Generuji soubory .m...');
    vars_list = {q, dq, m1, m2, m3, ...
                 I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz, ...
                 g_acc, L1, L2, L3, Lc1, Lc2, Lc3, d1, d2};

    matlabFunction(M_mat, 'File', 'get_M', 'Vars', vars_list);
    matlabFunction(C_mat, 'File', 'get_C', 'Vars', vars_list);
    matlabFunction(G_vec, 'File', 'get_G', 'Vars', vars_list);
else
    disp('Používám existující funkce get_M, get_C, get_G.');
end

%% 7. PŘÍPRAVA SIMULACE A DATA Z EXPERIMENTU
disp('Připravuji simulaci a načítám data...');

% --- A. DEFINICE PARAMETRŮ MECHANIKY ---
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;   
p.m2 = 0.0356 + 0.0059;  
p.m3 = 0.0435;
p.g = 9.81;

p.L1 = 0.201; p.d1 = 0.02; p.Lc1 = 0.01; 
p.L2 = 0.163;  p.d2 = 0.01; p.Lc2 = 0.081; 
p.L3 = 0.150; p.Lc3 = 0.075; 

p.I1 = [0.000 0.000 0.0051]; 
p.I2 = [0.00005 0.000318 0.000318];
p.I3 = [0.00005 0.000326 0.000326];

% --- B. DEFINICE PARAMETRŮ MOTORU A TŘENÍ ---
p.N  = 128/14;   
p.R  = 4.5;      
p.Kt = 0.05;     
p.Ke = 0.05;     

% Ladicí parametry pro shodu s experimentem
p.b = [0.0; 0.0005; 0.0005]; % Viskózní tření
p.fc = [0.0; 0.001; 0.001];  % Coulombovo tření


% --- C. NAČTENÍ DAT A PŘÍPRAVA VSTUPU ---
try
    load('experimenty.mat');
    disp('Data experimenty.mat úspěšně načtena.');
    
    % Použijeme dataset 'kyvadlo2'
    % Struktura: Col 1: Napětí, Col 2: Základna, Col 3: Kyv1, Col 4: Kyv2
    exp_data_struct = kyvadlo2; 
    
    t_exp = exp_data_struct.time;
    u_exp = exp_data_struct.signals.values(:,1);
    q_exp_raw = exp_data_struct.signals.values(:,2:4);
    
    % --- KONTROLA A PŘEVOD JEDNOTEK (Senzory vs Model) ---
    % Model pracuje v RADIÁNECH.
    % Pokud jsou data v grafu > 2*pi (cca 6.28), pravděpodobně jsou ve stupních nebo ticích.
    
    max_val = max(abs(q_exp_raw(:)));
    
    if max_val > 2000
        disp('Detekovány jednotky: IMPULZY ENKODÉRU (Ticks). Provádím převod na rad.');
        % Zde musíte doplnit rozlišení vašeho enkodéru (např. 2048 nebo 4096 na otáčku)
        enc_resolution = 2048 * 4; % PŘÍKLAD (ticks * kvadratura)
        q_exp_measured = (q_exp_raw / enc_resolution) * 2*pi;
        
    elseif max_val > 15
        disp('Detekovány jednotky: STUPNĚ. Provádím převod na rad.');
        q_exp_measured = deg2rad(q_exp_raw);
        
    else
        disp('Detekovány jednotky: RADIÁNY (nebo malý pohyb). Ponechávám.');
        q_exp_measured = q_exp_raw;
    end

    % 1. Interpolace vstupního napětí
    u_interpolant = griddedInterpolant(t_exp, u_exp, 'linear');
    p.u_func = @(t) u_interpolant(t);
    
    % 2. Počáteční podmínky
    q0 = q_exp_measured(1,:)';
    % Odhad počáteční rychlosti (numerická derivace prvních pár vzorků)
    dt_sample = t_exp(2) - t_exp(1);
    dq0 = (q_exp_measured(2,:)' - q_exp_measured(1,:)') / dt_sample;
    
    x0 = [q0; dq0];
    % ... (kód načítání dat) ...
    
    % --- KOREKCE SOUŘADNÝCH SYSTÉMŮ (OFFSETY) ---
    % q_sim = q_measured + offset
    % Zkuste tyto varianty pro q2 a q3:
    % A) Pokud model 0 = vodorovně a senzor 0 = svisle dolů -> offset = -pi/2
    % B) Pokud model 0 = vzpřímeně nahoru a senzor 0 = svisle dolů -> offset = pi (nebo -pi)
    
    % [Rameno, Kyvadlo1, Kyvadlo2]
    % Zkuste postupně: [0, -pi/2, 0] nebo [0, -pi/2, -pi/2] nebo [0, pi, pi]
    q_offset = [0; -pi; 0]; 
    q_offset2 = [0; 0; 0]; 
    
    % Aplikace offsetu na počáteční podmínku pro simulaci
    % (Simulátor dostane "svou" správnou polohu, aby gravitace fungovala)
    q0_sim = q_exp_measured(1,:)' + q_offset;
    q_exp_measured = q_exp_measured + q_offset2';
    disp('Otáčím znaménko u 3. kloubu (poslední rameno)...');
    %q_exp_measured(:, 3) = -q_exp_measured(:,2)+q_exp_measured(:, 3);
    dq0 = [0; 0; 0]; % Nebo vypočtené z derivace
    x0 = [q0_sim; dq0];
    
    % 3. Čas simulace
    tspan = [t_exp(1) t_exp(end)];
    
    has_experiment = true;
    
catch ME
    warning('Chyba při načítání experimentu: %s', ME.message);
    % Fallback pro testování bez dat
    p.u_func = @(t) 0; % Nulový vstup
    x0 = [0; 0; 0; 0; 0; 0];
    tspan = [0 t_exp(end)];
    has_experiment = false;
end

% --- D. SPUŠTĚNÍ INTEGRÁTORU ---
disp('Běží simulace...');
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
[t_out, x_out] = ode45(@(t,x) sys_dynamics(t, x, p), tspan, x0);

%% 8. VYKRESLENÍ A POROVNÁNÍ
if has_experiment
    figure('Name', 'Validace modelu (s korekcí offsetu) - 5s detail', 'Color', 'w');
    
    % Převedení simulace zpět do "senzorových" souřadnic pro vizualizaci
    q_sim_visual = x_out(:, 1:3) - (q_offset'); 
    
    % Nastavení limitu času
    time_limit = [0 5]; 
    
    % --- 1. Poloha q1 (Rameno) ---
    ax1 = subplot(4,1,1); hold on; % Ukládám si "handle" osy do ax1
    plot(t_exp, q_exp_measured(:,1), 'k-', 'LineWidth', 1.2); 
    plot(t_out, q_sim_visual(:,1), 'b-', 'LineWidth', 1.5);
    ylabel('phi [rad]'); title('Úhel ramene');
    legend('Měření', 'Simulace'); grid on;
    xlim(time_limit); % Omezení osy X
    
    % --- 2. Poloha q2 (Kyvadlo 1) ---
    ax2 = subplot(4,1,2); hold on;
    plot(t_exp, q_exp_measured(:,2), 'k-', 'LineWidth', 1.2);
    plot(t_out, q_sim_visual(:,2), 'r-', 'LineWidth', 1.5);
    ylabel('\theta_1 [rad]'); title('Úhel 1. kyvadla');
    grid on;
    xlim(time_limit); % Omezení osy X
    
    % --- 3. Poloha q3 (Kyvadlo 2) ---
    ax3 = subplot(4,1,3); hold on;
    plot(t_exp, q_exp_measured(:,3), 'k-', 'LineWidth', 1.2);
    plot(t_out, q_sim_visual(:,3), 'g-', 'LineWidth', 1.5);
    ylabel('\theta_2 [rad]'); title('Úhel 2. kyvadla');
    grid on;
    xlim(time_limit); % Omezení osy X
    
    % --- 4. Vstupní napětí ---
    ax4 = subplot(4,1,4); hold on;
    plot(t_exp, u_exp, 'm-', 'LineWidth', 1);
    ylabel('Napětí U [V]'); 
    xlabel('Čas [s]'); grid on;
    xlim(time_limit); % Omezení osy X
    
    % --- BONUS: Propojení os ---
    % Toto zajistí, že když v grafu zazoomujete myší, pohnou se všechny grafy najednou
    linkaxes([ax1, ax2, ax3, ax4], 'x');
end
%% --- FUNKCE DYNAMIKY ---
function dxdt = sys_dynamics(t, x, p)
    q_val = x(1:3);
    dq_val = x(4:6);
    
    % Rozbalení parametrů
    m1=p.m1; m2=p.m2; m3=p.m3;
    I1xx=p.I1(1); I1yy=p.I1(2); I1zz=p.I1(3);
    I2xx=p.I2(1); I2yy=p.I2(2); I2zz=p.I2(3);
    I3xx=p.I3(1); I3yy=p.I3(2); I3zz=p.I3(3);
    g_acc=p.g;
    L1=p.L1; L2=p.L2; L3=p.L3;
    Lc1=p.Lc1; Lc2=p.Lc2; Lc3=p.Lc3;
    d1=p.d1; d2=p.d2;
    
    % Matice
    M = get_M(q_val, dq_val, m1, m2, m3, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz, g_acc, L1, L2, L3, Lc1, Lc2, Lc3, d1, d2);
    C = get_C(q_val, dq_val, m1, m2, m3, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz, g_acc, L1, L2, L3, Lc1, Lc2, Lc3, d1, d2);
    G = get_G(q_val, dq_val, m1, m2, m3, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz, g_acc, L1, L2, L3, Lc1, Lc2, Lc3, d1, d2);
    
    % Síly
    tau_friction = (p.b .* dq_val) + (p.fc .* tanh(100 * dq_val));
    
    u_voltage = p.u_func(t); 
    
    % Saturace napětí (dle experimentu +/- 14.4V, pokud je potřeba)
    u_voltage = max(min(u_voltage, 14.4), -14.4);

    dq1 = dq_val(1);
    current = (u_voltage - p.N * p.Ke * dq1) / p.R; 
    tau_actuator_scalar = p.N * p.Kt * current;     
    
    tau_motor = [tau_actuator_scalar; 0; 0];
    
    % Řešení
    RHS = tau_motor - tau_friction - C*dq_val - G;
    ddq_val = M \ RHS;
    
    dxdt = [dq_val; ddq_val];
end