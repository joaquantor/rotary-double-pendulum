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
p.d_ram = 0;  % Vzdálenost těžiště ramene od osy motoru
p.d1 = 0.081;    % Vzdálenost těžiště 1. kyvadla od jeho osy rotace
p.d2 = 0.075;    % Vzdálenost těžiště 2. kyvadla od jeho osy rotace

% Momenty setrvačnosti [kg*m^2]
% I_xx: rotace kolem podélné osy profilu
% I_yy, I_zz: rotace kolem příčných os (hlavní setrvačnost při kyvu)
p.I1_zz = 0.0051; % Moment setrvačnosti ramene k ose Z (motoru)

p.I2_xx = 0;    % Setrvačnost 1. kyvadla k podélné ose
p.I2_yy = 0.000318;    % Setrvačnost 1. kyvadla k příčné ose
p.I2_zz = 0.000318;    % Setrvačnost 1. kyvadla k ose rotace (hlavní)

p.I3_xx = 0;    % Setrvačnost 2. kyvadla k podélné ose
p.I3_yy = 0.000326;    % Setrvačnost 2. kyvadla k příčné ose
p.I3_zz = 0.000326;    % Setrvačnost 2. kyvadla k ose rotace (hlavní)

% Ostatní konstanty
p.g = 9.81;     % Gravitační zrychlení [m/s^2]
p.tau = 0;      % Moment motoru (zatím 0 pro volný pohyb)

% Tření
p.b_phi = 0.0;     % Viskozní tření v ose motoru a ose ratační základny
p.b_theta1 = 0.0;  % Viskozní tření v 1. kloubu
p.b_theta2 = 0.0;  % Viskozní tření v 2. kloubu
% -------------------------------------------------------------------------

%% 2. NASTAVENÍ A SPUŠTĚNÍ SIMULACE
% Počáteční podmínky: [phi; theta1; theta2; dphi; dtheta1; dtheta2]
% Úhly v radiánech. theta = 0 je vodorovná poloha (nataženo).
% Svisle dolů je -pi/2 (dle naší DH konvence).
init_pos = [0; -pi/2 + 1; 0]; % Mírně vychýleno ze spodní polohy
init_vel = [0; 0; 0];              % Nulové počáteční rychlosti
y0 = [init_pos; init_vel];

t_span = [0 10]; % Čas simulace [s]

% Volání řešiče ODE45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t, y] = ode45(@(t, y) dynamics(t, y, p), t_span, y0, options);

%% 3. VIZUALIZACE VÝSLEDKŮ
figure;
subplot(3,1,1); plot(t, y(:,1)); title('Poloha ramene \phi'); ylabel('rad'); grid on;
subplot(3,1,2); plot(t, y(:,2)); title('Poloha 1. kyvadla \theta_1'); ylabel('rad'); grid on;
subplot(3,1,3); plot(t, y(:,3)); title('Poloha 2. kyvadla \theta_2'); ylabel('rad'); grid on; xlabel('Čas [s]');

% Jednoduchá animace (pohled shora/zboku)
animate_system(t, y, p);


%% 4. FUNKCE DYNAMIKY (Lagrangeovy rovnice)
function dydt = dynamics(t, y, p)
    % Rozbalení stavového vektoru
    q1 = y(1); q2 = y(2); q3 = y(3);       % phi, theta1, theta2
    dq1 = y(4); dq2 = y(5); dq3 = y(6);    % dphi, dtheta1, dtheta2
    
    dq = [dq1; dq2; dq3];
    
    % Zkrácené goniometrické funkce
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2 + q3); s12 = sin(q2 + q3);
    
    % --- MATICE SETRVAČNOSTI M(q) ---
    % Prvek M11 (Setrvačnost ramene + proměnná setrvačnost kyvadel vůči Z)
    % Pozor: Zahrnujeme Steinerovu větu a rotaci tenzorů
    M11 = p.I1_zz + p.m1 * p.d_ram^2 + ...
          p.m2 * (p.L0 + p.d1 * c1)^2 + ...
          p.m3 * (p.L0 + p.L1 * c1 + p.d2 * c12)^2 + ...
          (p.I2_xx * s1^2 + p.I2_yy * c1^2) + ...
          (p.I3_xx * s12^2 + p.I3_yy * c12^2);
    
    % Prvek M22 (Setrvačnost kyvadel v rovině kyvu)
    M22 = p.I2_zz + p.m2 * p.d1^2 + p.I3_zz + p.m3 * (p.L1^2 + p.d2^2) + ...
          2 * p.m3 * p.L1 * p.d2 * c2;
      
    % Prvek M33 (Setrvačnost 2. kyvadla)
    M33 = p.I3_zz + p.m3 * p.d2^2;
    
    % Prvek M23 = M32 (Vazba mezi kyvadly)
    M23 = p.I3_zz + p.m3 * p.d2^2 + p.m3 * p.L1 * p.d2 * c2;
    
    % Sestavení matice M
    M = [M11, 0, 0;
         0,   M22, M23;
         0,   M23, M33];
     
    % --- MATICE CORIOLISOVÝCH A ODSTŘEDIVÝCH SIL C(q, dq) ---
    % Pomocný člen h (derivace vazby podle q3)
    h = -p.m3 * p.L1 * p.d2 * s2;
    
    % Derivace M11 pro Coriolisovy členy v 1. rovnici
    % dM11/dq2 (theta1)
    dM11_dq2 = -2 * p.m2 * (p.L0 + p.d1*c1) * p.d1*s1 ...
               -2 * p.m3 * (p.L0 + p.L1*c1 + p.d2*c12) * (p.L1*s1 + p.d2*s12) ...
               + (p.I2_xx - p.I2_yy)*sin(2*q2) + (p.I3_xx - p.I3_yy)*sin(2*(q2+q3));
           
    % dM11/dq3 (theta2)
    dM11_dq3 = -2 * p.m3 * (p.L0 + p.L1*c1 + p.d2*c12) * p.d2*s12 ...
               + (p.I3_xx - p.I3_yy)*sin(2*(q2+q3));
    
    % Prvky matice C
    C11 = 0.5 * (dM11_dq2 * dq2 + dM11_dq3 * dq3); % Zde je chyba v definici C, oprava dle odvození:
    % Správně dle Christoffelových symbolů pro 1. řádek:
    % C11*dq1 + C12*dq2 + C13*dq3 = 0.5*dM11_dq2*dq1*dq2 ...
    % V maticovém násobení C*dq:
    
    C = zeros(3,3);
    % 1. řádek
    C(1,1) = 0.5 * dM11_dq2 * dq2 + 0.5 * dM11_dq3 * dq3;
    C(1,2) = 0.5 * dM11_dq2 * dq1;
    C(1,3) = 0.5 * dM11_dq3 * dq1;
    
    % 2. řádek
    C(2,1) = -0.5 * dM11_dq2 * dq1;
    C(2,2) = h * dq3;
    C(2,3) = h * (dq2 + dq3);
    
    % 3. řádek
    C(3,1) = -0.5 * dM11_dq3 * dq1;
    C(3,2) = -h * dq2;
    C(3,3) = 0;
    
    % --- VEKTOR GRAVITACE G(q) ---
    % Potenciální energie V = m*g*z
    % G = [dV/dphi; dV/dtheta1; dV/dtheta2]
    
    G = zeros(3,1);
    G(1) = 0; % Gravitace nepůsobí na rotaci ramene (pohybuje se vodorovně)
    G(2) = p.g * ((p.m2 * p.d1 + p.m3 * p.L1) * c1 + p.m3 * p.d2 * c12);
    G(3) = p.g * p.m3 * p.d2 * c12;
    
    % --- VEKTOR VNĚJŠÍCH SIL Q ---
    % Přidání tření: Q = Q_motor - Tření
    F_fric = [p.b_phi * dq1; p.b_theta1 * dq2; p.b_theta2 * dq3];
    Q = [p.tau; 0; 0] - F_fric;
    
    % --- VÝPOČET ZRYCHLENÍ ---
    % M * ddq + C * dq + G = Q  =>  ddq = M \ (Q - C*dq - G)
    
    ddq = M \ (Q - C*dq - G);
    
    dydt = [dq; ddq];
end

function animate_system(t, y, p)
    % Velmi jednoduchá animace pro kontrolu
    figure; hold on; grid on; axis equal;
    xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-0.5 0.5]);
    view(3);
    
    arm_line = plot3([0 0], [0 0], [0 0], 'k-', 'LineWidth', 2);
    pend1_line = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2);
    pend2_line = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2);
    
    for i = 1:1:length(t) % Přeskakování snímků pro rychlost
        q1 = y(i,1); q2 = y(i,2); q3 = y(i,3);
        
        % Kinematika bodů
        % Konec ramene (O1)
        x1 = p.L0 * cos(q1);
        y1 = p.L0 * sin(q1);
        z1 = 0;
        
        % Konec 1. kyvadla (O2)
        % V globálním systému (transformace T01 * T12)
        % Jednodušší: projekce
        % V lokálním rámu ramene je kyvadlo v rovině XZ (otočené o 90 deg?)
        % Dle DH: x2 směřuje podél kyvadla.
        
        % Vektor 1. kyvadla v souřadnicích ramene: [L1*cos(q2); 0; L1*sin(q2)]?
        % Ne, v DH T12: x os je osa kyvadla.
        % Zjednodušeně pro vizualizaci:
        % Kyvadlo visí "dolů" podél osy Z, když q2 = -pi/2
        
        % Poloha konce 1. kyvadla
        vec1 = [p.L1 * cos(q2); 0; p.L1 * sin(q2)]; % V rámu ramene (otočeném o q1)
        % Rotace o q1 kolem Z
        R01 = [cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
        % Ale pozor na DH konvenci, tam byla osa Z vodorovně pro kyvadlo.
        % Pro vizualizaci použijeme prostou geometrii:
        % Rameno je v rovině XY. Kyvadlo se kýve v rovině kolmé na rameno (obsahující Z a rameno).
        
        % Konec 1. kyvadla
        r_pend1 = [p.L1 * cos(q2); 0; p.L1 * sin(q2)]; % Lokálně v rovině kyvu
        % Transformace do globálu:
        % Osa x kyvu směřuje radiálně (podél ramene). Osa z kyvu je vodorovná tečná.
        % Takže rovina kyvu je svislá rovina procházející ramenem.
        % Souřadnice konce 1. kyvadla vůči konci ramene:
        dx2 = p.L1 * cos(q2) * cos(q1);
        dy2 = p.L1 * cos(q2) * sin(q1);
        dz2 = p.L1 * sin(q2);
        
        x2 = x1 + dx2;
        y2 = y1 + dy2;
        z2 = z1 + dz2;
        
        % Konec 2. kyvadla
        dx3 = p.L2 * cos(q2+q3) * cos(q1);
        dy3 = p.L2 * cos(q2+q3) * sin(q1);
        dz3 = p.L2 * sin(q2+q3);
        
        x3 = x2 + dx3;
        y3 = y2 + dy3;
        z3 = z2 + dz3;
        
        set(arm_line, 'XData', [0 x1], 'YData', [0 y1], 'ZData', [0 z1]);
        set(pend1_line, 'XData', [x1 x2], 'YData', [y1 y2], 'ZData', [z1 z2]);
        set(pend2_line, 'XData', [x2 x3], 'YData', [y2 y3], 'ZData', [z2 z3]);
        
        title(sprintf('Time: %.2f s', t(i)));
        drawnow;
    end
end