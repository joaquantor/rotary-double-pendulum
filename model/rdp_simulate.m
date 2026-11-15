function [t_sim, y_sim] = rdp_simulate(p, e, opts)
%RDP_SIMULATE Simulate the full system driven by the measured input voltage.
%   [t_sim, y_sim] = RDP_SIMULATE(p, e) uses the voltage and the initial
%   angles of experiment e (rdp_load_experiment) and returns the solution
%   on the measurement time grid.
if nargin < 3, opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-6); end
% The input comes from a Simulink zero-order hold -> piecewise constant
u_fun = griddedInterpolant(e.t, e.u, 'previous', 'nearest');
y0 = [e.phi(1); e.theta1(1); e.theta2(1); 0; 0; 0];
[t_sim, y_sim] = ode45(@(t, y) rdp_dynamics(t, y, p, u_fun), e.t, y0, opts);
end
