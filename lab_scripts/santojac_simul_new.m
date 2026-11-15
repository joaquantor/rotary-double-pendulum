%% Simulace Rotačního Dvojitého Kyvadla
% clear; clc; close all;

%% 1. PARAMETRY SYSTÉMU
% -------------------------------------------------------------------------
% Hmotnosti [kg]
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057;       
p.m2 = 0.0356 + 0.0059;  
p.m3 = 0.0435;  

% Délky [m]
p.L0 = 0.201;   
p.L1 = 0.163;   
p.L2 = 0.150;   
p.d_ram = -0.025;  
p.d1 = 0.078;    
p.d2 = 0.075;    

% Převodovka a setrvačnost
p.gear_ratio = 128 / 14; 
p.I1_zz = 0.0051; 

% Setrvačnosti kyvadel (lokální)
p.I2_xx = 0.00005;
p.I2_yy = 0.000318;    
p.I2_zz = 0.000318;    
p.I3_xx = 0.00005;    
p.I3_yy = 0.000326;    
p.I3_zz = 0.000326;    

% Tihove zrychleni
p.g = 9.81;     

% --- MOTOR A ŘÍZENÍ ---
p.U_pwm_amp = 14.4;  
p.U_ctrl_max = 10.0; 
p.R = 2.20;          
p.Kt = 0.052;        
p.Ke = 0.0516; % 5.4V/1000rpm
% Konstanty přepočtené na výstup (za převodovkou)
p.Ke_out = p.Ke * p.gear_ratio; 
p.Kt_out = p.Kt * p.gear_ratio;

% Testovací signál
p.ctrl_signal_amp = 8.0;  % Skok na 8V
p.ctrl_period = 2.0;      
p.ctrl_width = 0.05;       

% Tření
p.b_phi = 0.03;      
p.b_theta1 = 0.0003;  
p.b_theta2 = 0.0003;  
% -------------------------------------------------------------------------

%% 2. NASTAVENÍ A SPUŠTĚNÍ
% 0 definujeme jako svisle dolů.
init_pos = [0; 3.14; 1]; % Mírná výchylka, aby se to rozhoupalo
init_vel = [0; 0; 0];           
y0 = [init_pos; init_vel];

t_span = [0 30]; 

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t, y] = ode45(@(t, y) dynamics(t, y, p), t_span, y0, options);

%% 3. VIZUALIZACE
%% 3. VIZUALIZACE
figure('Name', 'Výsledek simulace "Obě kyvadla - hnané"');

% 1. Graf: Rameno
subplot(4,1,1); 
plot(t, y(:,1)); 
title('Rameno $\varphi$', 'Interpreter', 'latex'); 
grid on;
xlabel('t [s]'); % Popisek časové osy
ylabel('$\varphi$ [rad]', 'Interpreter', 'latex');

% 2. Graf: První kyvadlo
subplot(4,1,2); 
plot(t, y(:,2)); 
title('1. Kyvadlo \theta_1'); 
grid on;
xlabel('t [s]'); % Popisek časové osy
ylabel('$\theta_1$ [rad]', 'Interpreter', 'latex');

% 3. Graf: Druhé kyvadlo
subplot(4,1,3); 
plot(t, y(:,3)); 
title('2. Kyvadlo \theta_2'); 
grid on;
xlabel('t [s]'); % Popisek časové osy
ylabel('$\theta_2$ [rad]', 'Interpreter', 'latex');

% Výpočet elektrických veličin pro 4. graf
u_rec = zeros(size(t)); i_rec = zeros(size(t));
for k = 1:length(t)
    if mod(t(k), p.ctrl_period) < p.ctrl_width
        u = p.ctrl_signal_amp;
    else
        u = 0;
    end
    u_mot = (u/p.U_ctrl_max) * p.U_pwm_amp;
    i_act = (u_mot - p.Ke_out * y(k,4)) / p.R;
    u_rec(k) = u_mot; i_rec(k) = i_act;
end

% 4. Graf: Motor (Napětí a Proud)
subplot(4,1,4); 
yyaxis left; 
area(t, u_rec, 'FaceAlpha', 0.2); 
ylabel('U [V]');
yyaxis right; 
plot(t, i_rec); 
ylabel('I [A]'); 
title('Motor');
xlabel('t [s]'); % Popisek časové osy

animate_system(t, y, p);

%% 4. DYNAMIKA (Tangenciální vazba)
function dydt = dynamics(t, y, p)
    q1 = 0; q2 = y(2); q3 = y(3);       
    dq1 = 0; dq2 = y(5); dq3 = y(6);    
    dq = [dq1; dq2; dq3];
    
    % --- VSTUP (Motor) ---
    if mod(t, p.ctrl_period) < p.ctrl_width
        % MODE Obe kyvadla hnane
        % u_in = p.ctrl_signal_amp;
        % MODE Obe kyvadla volne
        u_in = 0;
    else
        u_in = 0;
    end
    u_in = max(min(u_in, p.U_ctrl_max), -p.U_ctrl_max);
    U_app = (u_in / p.U_ctrl_max) * p.U_pwm_amp;
    Current = (U_app - p.Ke_out * dq1) / p.R;
    Tau = p.Kt_out * Current;
    
    % --- MATICE M(q) PRO TANGENCIÁLNÍ MODEL ---
    c1 = cos(q2); s1 = sin(q2);
    c2 = cos(q3); s2 = sin(q3);
    c12 = cos(q2+q3); s12 = sin(q2+q3);
    
    % M11: Setrvačnost kolem osy motoru (zahrnuje roztažení kyvadel do boku)
    % r^2 = L0^2 + (vychylka_do_boku)^2
    % Vychylka 1: d1*s1
    % Vychylka 2: L1*s1 + d2*s12
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
    M21 = M12;
    M31 = M13;
    
    M = [M11, M12, M13;
         M21, M22, M23;
         M31, M32, M33];
     
    % --- CORIOLIS C(q,dq) ---
    % prostě vypočte coriolisovu matici numericky bez dlouheho vypisovani
    C = compute_coriolis_num(M, y, p);

    % --- GRAVITACE G(q) ---
    G = zeros(3,1);
    G(1) = 0;
    G(2) = p.g * ((p.m2*p.d1 + p.m3*p.L1)*s1 + p.m3*p.d2*s12);
    G(3) = p.g * p.m3*p.d2*s12;
    
    % Vnější síly
    Q = [Tau; 0; 0] - [p.b_phi*dq1; p.b_theta1*dq2; p.b_theta2*dq3];
    
    ddq = M \ (Q - C*dq - G);
    dydt = [dq; ddq];
end

% Pomocná funkce pro numerický výpočet Coriolisovy matice
function C = compute_coriolis_num(~, y, p)

    
    % Rozbalení stavů
    % y = [phi; theta1; theta2; dphi; dtheta1; dtheta2]
    q2 = y(2);      % theta1
    q3 = y(3);      % theta2
    dq1 = y(4);     % dphi
    dq2 = y(5);     % dtheta1
    dq3 = y(6);     % dtheta2
    
    % Goniometrické funkce
    s1 = sin(q2);   % sin(theta1)
    c1 = cos(q2);   % cos(theta1)
    s2 = sin(q3);   % sin(theta2)
    % c2 = cos(q3); % cos(theta2) - není přímo potřeba pro C, ale h využívá sin(q3)
    s12 = sin(q2+q3);
    c12 = cos(q2+q3);
    
    % --- 1. Pomocné proměnné pro derivace M11 ---
    % M11 závisí na q2 a q3 (přes "vyložení" do boku).
    % lambda1 = d(M11)/dq2
    % lambda2 = d(M11)/dq3
    
    % Horizontální výchylka H (bez L0) a její derivace
    H = p.L1*s1 + p.d2*s12;
    dH_dq2 = p.L1*c1 + p.d2*c12;
    dH_dq3 = p.d2*c12;
    
    % Derivace členu od m2 (term2 v dynamics): m2*(d1*s1)^2
    % d/dq2 [m2 * d1^2 * sin^2(q2)] = m2 * d1^2 * sin(2*q2)
    dM11_m2_dq2 = p.m2 * p.d1^2 * sin(2*q2);
    
    % Derivace členu od m3 (term3 v dynamics): m3*(L0^2 + H^2) -> zajímá nás jen H^2
    % d/dq2 [m3 * H^2] = 2 * m3 * H * dH/dq2
    dM11_m3_dq2 = 2 * p.m3 * H * dH_dq2;
    % d/dq3 [m3 * H^2] = 2 * m3 * H * dH/dq3
    dM11_m3_dq3 = 2 * p.m3 * H * dH_dq3;
    
    % Derivace rotační setrvačnosti (rot_inertia v dynamics)
    % I_eff = Iyy*sin^2 + Ixx*cos^2
    % d/dq [Iyy*sin^2(q) + Ixx*cos^2(q)] = (Iyy - Ixx) * sin(2*q)
    dM11_rot_dq2 = (p.I2_yy - p.I2_xx) * sin(2*q2) + ...
                   (p.I3_yy - p.I3_xx) * sin(2*(q2+q3)); % Chain rule d(q2+q3)/dq2 = 1
               
    dM11_rot_dq3 = (p.I3_yy - p.I3_xx) * sin(2*(q2+q3)); % Chain rule d(q2+q3)/dq3 = 1
    
    % Sumace pro lambdy
    lambda1 = dM11_m2_dq2 + dM11_m3_dq2 + dM11_rot_dq2;
    lambda2 = dM11_m3_dq3 + dM11_rot_dq3;
    
    % --- 2. Pomocná proměnná h (dynamika dvojkyvadla) ---
    % h = d(M23)/dq3 = d(M32)/dq3
    % M23 část závislá na q3 je: m3 * L1 * d2 * cos(q3)
    % Derivace cos(q3) je -sin(q3)
    h = -p.m3 * p.L1 * p.d2 * s2; 
    
    % --- 3. Sestavení matice C ---
    % Podle odvozeného vzorce:
    % C11 = 0.5*(lam1*dq2 + lam2*dq3)   C12 = 0.5*lam1*dq1   C13 = 0.5*lam2*dq1
    % C21 = -0.5*lam1*dq1               C22 = h*dq3          C23 = h*(dq2+dq3)
    % C31 = -0.5*lam2*dq1               C32 = -h*dq2         C33 = 0
    
    C = zeros(3,3);
    
    % 1. řádek (Rovnice pro phi - Coriolis od kyvadel)
    C(1,1) = 0.5 * (lambda1 * dq2 + lambda2 * dq3);
    C(1,2) = 0.5 * lambda1 * dq1;
    C(1,3) = 0.5 * lambda2 * dq1;
    
    % 2. řádek (Rovnice pro theta1 - Odstředivá od phi + vazba s theta2)
    C(2,1) = -0.5 * lambda1 * dq1;
    C(2,2) = h * dq3;
    C(2,3) = h * (dq2 + dq3);
    
    % 3. řádek (Rovnice pro theta2 - Odstředivá od phi + vazba s theta1)
    C(3,1) = -0.5 * lambda2 * dq1;
    C(3,2) = -h * dq2;
    C(3,3) = 0;
end

function animate_system(t, y, p)

    f = figure('Name', 'Animace Kyvadla (Tangenciální)', 'Color', 'w');
    ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    view(ax, 45, 30); 
    limit = p.L0 + p.L1 + p.L2 + 0.1;
    xlim(ax, [-limit limit]); ylim(ax, [-limit limit]); zlim(ax, [-limit limit]);
    
    plot3(ax, [0 0], [0 0], [-0.1 0], 'k-', 'LineWidth', 4); 
    h_arm = plot3(ax, [0 0], [0 0], [0 0], 'k-o', 'LineWidth', 3, 'MarkerFaceColor', 'k');
    h_pen1 = plot3(ax, [0 0], [0 0], [0 0], 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    h_pen2 = plot3(ax, [0 0], [0 0], [0 0], 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    
    tic; sim_time = 0;
    while ishandle(f) && sim_time <= t(end)
        sim_time = toc; if sim_time > t(end), break; end
        cur_q = interp1(t, y, sim_time);
        q1 = cur_q(1); q2 = cur_q(2); q3 = cur_q(3);
        
        x1 = p.L0*cos(q1); y1 = p.L0*sin(q1); z1 = 0;
        
        y_loc_1 = p.L1*sin(q2); z_loc_1 = -p.L1*cos(q2); % 0=dolů
        x2 = x1 - y_loc_1*sin(q1); y2 = y1 + y_loc_1*cos(q1); z2 = z1 + z_loc_1;
        
        y_loc_2 = p.L1*sin(q2) + p.L2*sin(q2+q3);
        z_loc_2 = -p.L1*cos(q2) - p.L2*cos(q2+q3);
        x3 = x1 - y_loc_2*sin(q1); y3 = y1 + y_loc_2*cos(q1); z3 = z1 + z_loc_2;
        
        set(h_arm, 'XData', [0 x1], 'YData', [0 y1], 'ZData', [0 z1]);
        set(h_pen1, 'XData', [x1 x2], 'YData', [y1 y2], 'ZData', [z1 z2]);
        set(h_pen2, 'XData', [x2 x3], 'YData', [y2 y3], 'ZData', [z2 z3]);
        drawnow;
    end
end