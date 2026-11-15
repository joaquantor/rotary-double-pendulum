function p = rdp_params()
%RDP_PARAMS Physical parameters of the rotary double pendulum.
%   Geometry and masses were measured on the lab rig; the inertias of the
%   pendulum links and all friction coefficients were identified from
%   measured data (see identification/).
%
%   Generalized coordinates:  q = [phi; theta1; theta2]
%     phi    - rotation of the horizontal arm (motor axis, vertical)
%     theta1 - 1st pendulum, measured from the downward vertical
%     theta2 - 2nd pendulum, RELATIVE to the 1st one

% --- Masses [kg] ---------------------------------------------------------
p.m1 = 0.4362 + 0.110 + 2*0.020 + 0.2057; % arm: counterweight + bar + 2 encoders + gear
p.m2 = 0.0356 + 0.0059;                   % 1st pendulum + connecting rivet
p.m3 = 0.0435;                            % 2nd pendulum

% --- Geometry [m] --------------------------------------------------------
p.L0 = 0.201;      % motor axis -> 1st pendulum joint
p.L1 = 0.163;      % 1st pendulum joint -> 2nd pendulum joint
p.L2 = 0.150;      % length of the 2nd pendulum
p.d_ram = -0.0638; % arm centre of mass offset (identified)
p.d1 = 0.078;      % centre of mass of the 1st pendulum from its joint
p.d2 = 0.075;      % centre of mass of the 2nd pendulum from its joint

% --- Inertias [kg m^2] (principal, about the centre of mass) -------------
% x = along the link, y/z = perpendicular to it (z = joint axis).
p.I1_zz = 0.0051;
p.I2_xx = 0.00005; p.I2_zz = 1.283e-4; p.I2_yy = p.I2_zz; % identified
p.I3_xx = 0.00005; p.I3_zz = 8.91e-5;  p.I3_yy = p.I3_zz; % identified

p.g = 9.81;

% --- DC motor + gearbox + PWM power stage -------------------------------
p.gear_ratio = 128 / 14;
p.U_pwm_amp  = 14.4;   % PWM supply amplitude [V]
p.U_ctrl_max = 10.0;   % control-signal range [V]
p.R  = 2.20;           % armature resistance [Ohm]
p.Kt = 0.052;          % torque constant [Nm/A]
p.Ke = 0.0516;         % back-EMF constant [V s/rad]
p.Kt_out = p.Kt * p.gear_ratio;   % referred to the gearbox output
p.Ke_out = p.Ke * p.gear_ratio;

% --- Friction (identified) ----------------------------------------------
p.b_phi    = 0.1714;   % arm viscous friction [N m s/rad]
p.b_theta1 = 0;        % joint 1 viscous friction [N m s/rad]
p.b_theta2 = 1.705e-4; % joint 2 viscous friction [N m s/rad]
p.c_theta1 = 5.70e-4;  % joint 1 Coulomb friction [N m]
p.c_theta2 = 4.06e-4;  % joint 2 Coulomb friction [N m]
p.v_eps    = 0.02;     % smoothing of sign(): tanh(dq/v_eps) [rad/s]
end
