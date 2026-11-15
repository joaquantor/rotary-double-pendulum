%% Simulace Rotačního Dvojitého Kyvadla
clear; clc; close all;

%% 1. PARAMETRY SYSTÉMU (Zde doplňte vaše vypočtené hodnoty)
% -------------------------------------------------------------------------
% Hmotnosti [kg]
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;       % Hmotnost rotujícího ramene (včetně protizávaží, velkého ozubeného kola)
p.m2 = 0.0356 + 0.0059;  % Hmotnost 1. kyvadla (viz zadání: 35.6 g) + hmotnost propojovacích nýtů
p.m3 = 0.0435;  % Hmotnost 2. kyvadla

% Délky a geometrie [m]
p.L0 = 0.201;   % Délka ramene od osy motoru k ose 1. kyvadla
p.L1 = 0.163;   % Délka 1. kyvadla (od osy rotace po "konec" ramene)
p.L2 = 0.150;   % Délka 2. kyvadla (od osy rotace po "konec" ramene)
p.d_ram = 0;    % Vzdálenost těžiště ramene od osy motoru (řekli jsme že pro správnou funkčnost by bylo dobré, aby soustava byla vyvážená)
p.d1 = 0.081;   % Vzdálenost těžiště 1. kyvadla od jeho osy rotace
p.d2 = 0.075;   % Vzdálenost těžiště 2. kyvadla od jeho osy rotace

% Momenty setrvačnosti [kg*m^2]
p.I1_zz = 0.0051; % Moment setrvačnosti ramene k ose Z (motoru)
p.I2_xx = 0;    % Setrvačnost 1. kyvadla k podélné ose
p.I2_yy = 0.000318;    % Setrvačnost 1. kyvadla k příčné ose
p.I2_zz = 0.000318;    % Setrvačnost 1. kyvadla k ose rotace (hlavní)
p.I3_xx = 0;    % Setrvačnost 2. kyvadla k podélné ose
p.I3_yy = 0.000326;    % Setrvačnost 2. kyvadla k příčné ose
p.I3_zz = 0.000326;    % Setrvačnost 2. kyvadla k ose rotace (hlavní)

% Ostatní konstanty
p.g = 9.81;     % Gravitační zrychlení [m/s^2]

% --- NOVÉ: NASTAVENÍ IMPULZNÍHO BUZENÍ MOTOREM ---
% Definice periodického "kopání" do motoru
p.tau_amp = 0.15;     % Amplituda impulzu [Nm] (síla kopnutí)
p.tau_period = 2.0;   % Perioda: impulz přijde každé 2 sekundy
p.tau_width = 0.1;    % Šířka: impulz trvá 0.1 sekundy

% Tření
p.b_phi = 0.1;     % Viskozní tření v ose motoru (zvýšeno, aby se systém po impulzu uklidnil)
p.b_theta1 = 0.0001;  % Viskozní tření v 1. kloubu
p.b_theta2 = 0.0001;  % Viskozní tření v 2. kloubu
% -------------------------------------------------------------------------

%% 2. NASTAVENÍ A SPUŠTĚNÍ SIMULACE
% Počáteční podmínky: [phi; theta1; theta2; dphi; dtheta1; dtheta2]
init_pos = [0; -pi/2 + 0.1; 0]; % Mírně vychýleno ze spodní polohy
init_vel = [0; 0; 0];           % Nulové počáteční rychlosti
y0 = [init_pos; init_vel];

t_span = [0 10]; % Čas simulace [s]

% Volání řešiče ODE45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t, y] = ode45(@(t, y) dynamics(t, y, p), t_span, y0, options);

%% 3. VIZUALIZACE VÝSLEDKŮ
figure('Name', 'Průběhy veličin');

% 1. Poloha ramene
subplot(4,1,1); plot(t, y(:,1)); 
title('Poloha ramene \phi'); ylabel('rad'); grid on;

% 2. Poloha 1. kyvadla
subplot(4,1,2); plot(t, y(:,2)); 
title('Poloha 1. kyvadla \theta_1'); ylabel('rad'); grid on;

% 3. Poloha 2. kyvadla
subplot(4,1,3); plot(t, y(:,3)); 
title('Poloha 2. kyvadla \theta_2'); ylabel('rad'); grid on;

% 4. Průběh momentu (NOVÉ: rekonstrukce pro graf)
% Protože ODE45 má proměnný krok, musíme zpětně dopočítat, kdy byl moment aktivní
tau_rec = zeros(size(t));
for i = 1:length(t)
    if mod(t(i), p.tau_period) < p.tau_width
        tau_rec(i) = p.tau_amp;
    else
        tau_rec(i) = 0;
    end
end
subplot(4,1,4); area(t, tau_rec, 'FaceColor', 'r', 'FaceAlpha', 0.3);
title('Vstupní moment motoru \tau (Impulzy)'); ylabel('Nm'); xlabel('Čas [s]'); grid on;

% Jednoduchá animace (pohled shora/zboku)
animate_system(t, y, p);


%% 4. FUNKCE DYNAMIKY (Lagrangeovy rovnice)
function dydt = dynamics(t, y, p)
    % Rozbalení stavového vektoru
    q1 = y(1); q2 = y(2); q3 = y(3);       % phi, theta1, theta2
    dq1 = y(4); dq2 = y(5); dq3 = y(6);    % dphi, dtheta1, dtheta2
    
    dq = [dq1; dq2; dq3];
    
    % --- NOVÉ: VÝPOČET IMPULZNÍHO MOMENTU ---
    % Funkce mod(t, perioda) vrací pilovitý průběh času (0 -> perioda -> 0).
    % Pokud je tento "lokální čas" menší než šířka impulzu, motor zabere.
    if mod(t, p.tau_period) < p.tau_width
        tau_current = p.tau_amp;
    else
        tau_current = 0;
    end
    % ----------------------------------
    
    % Zkrácené goniometrické funkce
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2 + q3); s12 = sin(q2 + q3);
    
    % --- MATICE SETRVAČNOSTI M(q) ---
    M11 = p.I1_zz + p.m1 * p.d_ram^2 + ...
          p.m2 * (p.L0 + p.d1 * c1)^2 + ...
          p.m3 * (p.L0 + p.L1 * c1 + p.d2 * c12)^2 + ...
          (p.I2_xx * s1^2 + p.I2_yy * c1^2) + ...
          (p.I3_xx * s12^2 + p.I3_yy * c12^2);
    
    M22 = p.I2_zz + p.m2 * p.d1^2 + p.I3_zz + p.m3 * (p.L1^2 + p.d2^2) + ...
          2 * p.m3 * p.L1 * p.d2 * c2;
      
    M33 = p.I3_zz + p.m3 * p.d2^2;
    M23 = p.I3_zz + p.m3 * p.d2^2 + p.m3 * p.L1 * p.d2 * c2;
    
    M = [M11, 0, 0;
         0,   M22, M23;
         0,   M23, M33];
     
    % --- MATICE CORIOLISOVÝCH A ODSTŘEDIVÝCH SIL C(q, dq) ---
    h = -p.m3 * p.L1 * p.d2 * s2;
    
    dM11_dq2 = -2 * p.m2 * (p.L0 + p.d1*c1) * p.d1*s1 ...
               -2 * p.m3 * (p.L0 + p.L1*c1 + p.d2*c12) * (p.L1*s1 + p.d2*s12) ...
               + (p.I2_xx - p.I2_yy)*sin(2*q2) + (p.I3_xx - p.I3_yy)*sin(2*(q2+q3));
           
    dM11_dq3 = -2 * p.m3 * (p.L0 + p.L1*c1 + p.d2*c12) * p.d2*s12 ...
               + (p.I3_xx - p.I3_yy)*sin(2*(q2+q3));
    
    C = zeros(3,3);
    C(1,1) = 0.5 * dM11_dq2 * dq2 + 0.5 * dM11_dq3 * dq3;
    C(1,2) = 0.5 * dM11_dq2 * dq1;
    C(1,3) = 0.5 * dM11_dq3 * dq1;
    
    C(2,1) = -0.5 * dM11_dq2 * dq1;
    C(2,2) = h * dq3;
    C(2,3) = h * (dq2 + dq3);
    
    C(3,1) = -0.5 * dM11_dq3 * dq1;
    C(3,2) = -h * dq2;
    C(3,3) = 0;
    
    % --- VEKTOR GRAVITACE G(q) ---
    G = zeros(3,1);
    G(1) = 0; % Gravitace nepůsobí na rotaci ramene (pohybuje se vodorovně)
    G(2) = p.g * ((p.m2 * p.d1 + p.m3 * p.L1) * c1 + p.m3 * p.d2 * c12);
    G(3) = p.g * p.m3 * p.d2 * c12;
    
    % --- VEKTOR VNĚJŠÍCH SIL Q ---
    % Přidání tření: Q = Q_motor - Tření
    % Zde používáme proměnnou tau_current vypočtenou na začátku funkce
    F_fric = [p.b_phi * dq1; p.b_theta1 * dq2; p.b_theta2 * dq3];
    Q = [tau_current; 0; 0] - F_fric;
    
    % --- VÝPOČET ZRYCHLENÍ ---
    ddq = M \ (Q - C*dq - G);
    dydt = [dq; ddq];
end

function animate_system(t, y, p)
    figure('Name', 'Animace'); hold on; grid on; axis equal;
    % Fixní limity os pro přehlednost
    limit = p.L0 + p.L1 + p.L2 + 0.1;
    xlim([-limit limit]); ylim([-limit limit]); zlim([-limit limit]);
    view(3); xlabel('X'); ylabel('Y'); zlabel('Z');
    
    arm_line = plot3([0 0], [0 0], [0 0], 'k-', 'LineWidth', 3);
    pend1_line = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2);
    pend2_line = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2);
    joint1 = plot3(0,0,0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    
    title('Simulace (Zavřete okno pro ukončení)');
    
    step = max(1, floor(length(t)/200)); 
    
    for i = 1:step:length(t)
        if ~ishandle(arm_line), break; end 
        
        q1 = y(i,1); q2 = y(i,2); q3 = y(i,3);
        
        % Kinematika bodů (stejná jako v předchozí verzi)
        x1 = p.L0 * cos(q1);
        y1 = p.L0 * sin(q1);
        z1 = 0;
        
        % O2 (Konec 1. kyvadla)
        x_p1_loc = p.L1 * cos(q2); 
        z_p1_loc = p.L1 * sin(q2);
        
        x2 = x1 + x_p1_loc * cos(q1);
        y2 = y1 + x_p1_loc * sin(q1);
        z2 = z1 + z_p1_loc;
        
        % O3 (Konec 2. kyvadla)
        x_p2_loc = p.L2 * cos(q2+q3);
        z_p2_loc = p.L2 * sin(q2+q3);
        
        x3 = x2 + x_p2_loc * cos(q1);
        y3 = y2 + x_p2_loc * sin(q1);
        z3 = z2 + z_p2_loc;
        
        set(arm_line, 'XData', [0 x1], 'YData', [0 y1], 'ZData', [0 z1]);
        set(pend1_line, 'XData', [x1 x2], 'YData', [y1 y2], 'ZData', [z1 z2]);
        set(pend2_line, 'XData', [x2 x3], 'YData', [y2 y3], 'ZData', [z2 z3]);
        set(joint1, 'XData', x1, 'YData', y1, 'ZData', z1);
        
        title(sprintf('Čas: %.2f s', t(i)));
        drawnow limitrate;
    end
end