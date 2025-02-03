%% Simulation of B-TENG (Biomimetic Triboelectric Nanogenerator)
% It includes:
%   1. Energy harvesting calculations (Equations 7–11)
%   2. Fluid–structure interaction and nonlinear film dynamics (Equations 1–6)
%   3. Electrical output simulation vs. wind velocity and load resistance.
%
% Additionally, it generates MATLAB-based figures that simulate:
%   - The Bernoulli effect and two interacting (flapping) films.
%   - The working principle: electric potential distribution versus contact area variation.
%   - Output performance plots (voltage, current, power, conversion efficiency).
%

clear; clc; close all;
warning('off', 'all');
%% ------------------------------------------------------------------------
%% Section 1: Energy Harvesting Calculations (Equations 7–11)
%% ------------------------------------------------------------------------
% Equations:
%   (7) Q_v = v * S         --> Volume flow rate
%   (8) Q_m = Q_v * rho     --> Mass flow rate
%   (9) P_o = V^2 / R       --> Instantaneous output power
%  (10) P_i = (Q_m * v^2)/2  --> Instantaneous input power
%  (11) ? = (P_o / P_i)*100  --> Conversion efficiency

% Constants
v_fixed = 8;          % Wind velocity in m/s (nominal value)
S = 0.02;             % Air outlet area in m^2
rho_air = 1.225;      % Air density in kg/m^3
R_load = 1e6;         % Load resistance in ohms
V_oc = 175;           % Open-circuit voltage (V)
I_sc = 43e-3;         % Short-circuit current in A (43 mA)

% Equation 7: Volume Flow Rate
Q_v = v_fixed * S;
fprintf('Equation 7 - Q_v: %.4f m^3/s\n', Q_v);

% Equation 8: Mass Flow Rate
Q_m = Q_v * rho_air;
fprintf('Equation 8 - Q_m: %.4f kg/s\n', Q_m);

% Equation 9: Instantaneous Output Power
P_o = (V_oc^2) / R_load;
fprintf('Equation 9 - P_o: %.4f W\n', P_o);

% Equation 10: Instantaneous Input Power
P_i = (Q_m * v_fixed^2) / 2;
fprintf('Equation 10 - P_i: %.4f W\n', P_i);

% Equation 11: Conversion Efficiency
eta_fixed = (P_o / P_i) * 100;
fprintf('Equation 11 - Conversion Efficiency: %.2f%%\n\n', eta_fixed);

%% ------------------------------------------------------------------------
%% Section 2: Fluid–Structure Interaction & Film Dynamics (Equations 1–6)
%% ------------------------------------------------------------------------
% The film (or flexible PVDF) is modeled with a geometrically nonlinear equation.
% Equations:
%   (3) rs*(d^2X/dt^2) - d/ds(T*dX/ds) + Kb*(d^4X/ds^4) = F
%   (4) T(s) = Ks*(|dX/ds| - 1)
%   (6) Nondimensional parameters: M = rs/(rho_f*L), Ks^, Kb^

% Filament and fluid properties
rho_s = 1.2;        % Linear density (kg/m)
K_s = 1e3;          % Stretching coefficient (N/m)
K_b = 1e-2;         % Flexural rigidity (N*m^2)
L = 0.08;           % Characteristic length (m)
U = 8;              % Characteristic velocity (m/s)
rho_f = rho_air;    % Fluid density

% Equation 6: Nondimensional parameters
M = rho_s / (rho_f * L);
K_s_dimless = K_s / (rho_f * U^2 * L);
K_b_dimless = K_b / (rho_f * U^2 * L^3);
fprintf('Equation 6 - Non-dimensional parameters:\n  M = %.4f\n  K_s_dimless = %.4f\n  K_b_dimless = %.4f\n\n',...
    M, K_s_dimless, K_b_dimless);

% Define the spatial coordinate along the film (Lagrangian coordinate)
s = linspace(0, L, 100);

% An initial sinusoidal displacement profile X(s)
X = sin(2 * pi * s / L);

% Compute spatial derivatives using finite differences
dX_ds = gradient(X, s);
d2X_ds2 = gradient(dX_ds, s);
d4X_ds4 = gradient(gradient(d2X_ds2, s), s);

% Equation 4: Tensile Stress
T = K_s * (abs(dX_ds) - 1);

% Equation 3: Force Balance on the Film
F = rho_s * d2X_ds2 - gradient(T .* dX_ds, s) + K_b * d4X_ds4;
fprintf('Equation 3 - Computed film force along its length.\n');

% Equation 5: Immersed Boundary Method (using a smooth delta function)
sigma = 0.005;  % Standard deviation for the Gaussian (adjust as needed)
delta_smooth = @(x) exp(-x.^2/(2*sigma^2)) / (sigma*sqrt(2*pi));

% Restrict integration to the region where the delta function is significant.
s_center = 0.5 * L;
s_lower = s_center - 3*sigma;
s_upper = s_center + 3*sigma;

f_density = -integral(@(s_val) interp1(s, F, s_val) .* delta_smooth(s_val - s_center), ...
                      s_lower, s_upper, 'ArrayValued', true);
fprintf('Equation 5 - Force density: %.4f N/m\n\n', f_density);

%% ------------------------------------------------------------------------
%% Section 3: Fluid Dynamics (Equations 1–2)
%% ------------------------------------------------------------------------
% Using a simplified 1D form for demonstration:
%   Equation 1: ? (?u/?t + u·?u) = -?p + ??^2u + f_ext
%   Equation 2: ?·u = 0

u = 1;              % fluid velocity (m/s)
p = 101325;         % Atmospheric pressure (Pa)
mu = 1.81e-5;       % Dynamic viscosity (Pa·s)
f_external = 0;     % No external force

du_dt = -gradient(p, s) + mu * gradient(gradient(u, s), s) + f_external;
fprintf('Equation 1 - Average fluid acceleration: %.4f m/s^2\n', mean(du_dt));

div_u = gradient(u, s);
fprintf('Equation 2 - Average velocity divergence: %.4f 1/s\n\n', mean(div_u));

%% ------------------------------------------------------------------------
%% Section 4: Simulation of Electrical Output vs. Wind Velocity & Load
%% ------------------------------------------------------------------------
% The electrical output is influenced by the contact–propagation–separation
% cycle of the films. We assume a piecewise linear behavior:
%   For wind speeds from 4 to 8 m/s, voltage and current increase;
%   Above 8 m/s, they decrease.
%
% For each wind speed, we calculate:
%   - Output voltage V (V)
%   - Output current I (mA)
%   - Output power P_o (W) [Equation 9]
%   - Input power P_i (W) [Equation 10]
%   - Conversion efficiency ? (%) [Equation 11]

v_array = linspace(4, 15, 100);  % Wind velocity range (m/s)
f_array = 5 + 2 * v_array;        % Example motion frequency (Hz)

% Preallocate arrays for electrical outputs
V_array = zeros(size(v_array));  % Output voltage (V)
I_array = zeros(size(v_array));  % Output current (mA)
Po_array = zeros(size(v_array)); % Output power (W)
Pi_array = zeros(size(v_array)); % Input power (W)
eta_array = zeros(size(v_array));% Conversion efficiency (%)

for i = 1:length(v_array)
    v_val = v_array(i);
    if v_val <= 8
        V_array(i) = 115 + (V_oc - 115) * (v_val - 4) / (8 - 4);
        I_array(i) = 15 + (43 - 15) * (v_val - 4) / (8 - 4);
    else
        V_array(i) = V_oc - (V_oc - 115) * (v_val - 8) / (15 - 8);
        I_array(i) = 43 - (43 - 15) * (v_val - 8) / (15 - 8);
    end
    Po_array(i) = (V_array(i)^2) / R_load;
    
    Qv = v_val * S;
    Qm = Qv * rho_air;
    Pi_array(i) = (Qm * v_val^2) / 2;
    
    eta_array(i) = (Po_array(i) / Pi_array(i)) * 100;
end

%% Frequency Analysis Using FFT
fs = 1000;                  % Sampling frequency (Hz)
t_signal = 0:1/fs:1-1/fs;   % 1-second time vector
flap_frequency = 8;         % Example flapping frequency (Hz)
signal = sin(2*pi*flap_frequency*t_signal) + 0.1*randn(size(t_signal));
Y = fft(signal);
P2 = abs(Y/fs);
P1 = P2(1:length(t_signal)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f_fft = fs*(0:(length(t_signal)/2))/length(t_signal);

%% ------------------------------------------------------------------------
%% Section 5: Plotting Main Simulation Results
%% ------------------------------------------------------------------------

% ----- Figure 1: Output Electrical Characteristics vs. Wind Velocity -----
figure;
subplot(3,1,1);
plot(v_array, V_array, 'b', 'LineWidth', 2);
xlabel('Wind Velocity (m/s)'); ylabel('Output Voltage (V)');
title('Output Voltage vs. Wind Velocity'); grid on;

subplot(3,1,2);
plot(v_array, I_array, 'r', 'LineWidth', 2);
xlabel('Wind Velocity (m/s)'); ylabel('Output Current (mA)');
title('Output Current vs. Wind Velocity'); grid on;

subplot(3,1,3);
plot(v_array, f_array, 'g', 'LineWidth', 2);
xlabel('Wind Velocity (m/s)'); ylabel('Motion Frequency (Hz)');
title('Motion Frequency vs. Wind Velocity'); grid on;

% ----- Figure 2: Measured (Simulated) Data: Voltage & Power vs. Wind Speed -----
wind_speeds_exp = [4, 6, 8, 10, 15];
voltages_exp = [115, 150, 175, 160, 140]; % Example measured voltages (V)
power_exp = (voltages_exp.^2) / R_load;      % Power (W)
power_exp_mW = power_exp * 1e3;              % Power (mW)

figure;
subplot(2,1,1);
plot(wind_speeds_exp, voltages_exp, '-o', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Wind Speed (m/s)'); ylabel('Voltage (V)');
title('Measured Voltage vs. Wind Speed'); grid on;

subplot(2,1,2);
plot(wind_speeds_exp, power_exp_mW, '-o', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Wind Speed (m/s)'); ylabel('Power Output (mW)');
title('Measured Power Output vs. Wind Speed'); grid on;

% ----- Figure 3: Output vs. Load Resistance -----
R_vals = logspace(3, 7, 100);  % Load resistance range: 1 k? to 10 M?
R_internal = 1e6;              % Internal resistance (1 M?)
V_out = V_oc .* R_vals ./ (R_vals + R_internal);  % Voltage divider model
I_out = V_oc ./ (R_vals + R_internal);
P_out = (V_out.^2) ./ R_vals;  % Output power (W)

figure;
subplot(2,1,1);
semilogx(R_vals, V_out, 'b', 'LineWidth', 2);
xlabel('Load Resistance (\Omega)'); ylabel('Output Voltage (V)');
title('Output Voltage vs. Load Resistance'); grid on;

subplot(2,1,2);
semilogx(R_vals, P_out*1e3, 'k', 'LineWidth', 2); % mW
xlabel('Load Resistance (\Omega)'); ylabel('Output Power (mW)');
title('Output Power vs. Load Resistance'); grid on;
hold on;
[maxP, idx_max] = max(P_out);
plot(R_vals(idx_max), maxP*1e3, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('Output Power','Maximum Power');

% ----- Figure 4: Conversion Efficiency vs. Wind Velocity -----
figure;
plot(v_array, eta_array, 'm', 'LineWidth', 2);
xlabel('Wind Velocity (m/s)'); ylabel('Conversion Efficiency (%)');
title('Conversion Efficiency vs. Wind Velocity'); grid on;

% ----- Figure 5: Filament Force Distribution along the Film -----
figure;
plot(s, F, 'LineWidth', 2);
xlabel('Lagrangian Coordinate s (m)'); ylabel('Force F (N)');
title('Force Distribution along the Filament'); grid on;

% ----- Figure 6: FFT Spectrum of the Simulated Flapping Signal -----
figure;
plot(f_fft, P1, 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('FFT Spectrum of Flapping Signal (Simulated)'); grid on;

%% ------------------------------------------------------------------------
%% Section 6: MATLAB Simulation of the B-TENG Working Principle
% (Electric Potential Distribution vs. Contact Area Variation)
%% ------------------------------------------------------------------------
% According to the article, the B-TENG output is generated by the contact–propagation–separation
% cycle. The maximum contact area occurs at d/L = 0.25 and primary contact at about 1/3 of the film.
%
% In this simulation:
%   - We assume the effective contact area varies sinusoidally with time.
%   - The electric potential along the film is modeled as a Gaussian profile
%     centered at 1/3 of the film length.
%   - The instantaneous electric potential is given by the product of the spatial profile and
%     a time-dependent modulation function.
%
% Parameters for the working principle simulation:
T_period = 1;      % Period of one cycle (s)
V_max = 175;       % Maximum potential (V)
sigma_profile = 0.1;       % Width (standard deviation) of Gaussian (normalized units)
center = 1/3;      % Primary contact position (normalized)

% Spatial coordinate (normalized from 0 to 1)
x = linspace(0, 1, 200);

% Time vector for one cycle
t_cycle = linspace(0, T_period, 100);
[X_grid, T_grid] = meshgrid(x, t_cycle);

% Modulation: effective contact area variation (sinusoidal, ranging 0 to 1)
modulation = (sin(2*pi*T_grid/T_period) + 1) / 2;

% Gaussian spatial profile for electric potential distribution
potential_distribution = V_max * exp(-((X_grid - center).^2) / (2*sigma_profile^2));

% Effective electric potential: product of modulation and spatial profile
V_elec = modulation .* potential_distribution;

% ----- Figure 7: Snapshot of Electric Potential Distribution at Mid-cycle -----
figure;
plot(x, V_elec(round(end/2),:), 'LineWidth', 2);
xlabel('Normalized Position along Film');
ylabel('Electric Potential (V)');
title('Electric Potential Distribution at t = T/2');
grid on;

% ----- Figure 8: Surface Plot of Electric Potential Distribution over One Cycle -----
figure;
surf(x, t_cycle, V_elec, 'EdgeColor', 'none');
xlabel('Normalized Position along Film');
ylabel('Time (s)');
zlabel('Electric Potential (V)');
title('Electric Potential Distribution over One Cycle');
colorbar;
view(45,30);

%% ------------------------------------------------------------------------
%% Section 7: Additional Simulated Figures from the Article
%% ------------------------------------------------------------------------
% --- Figure 9: Simulated Flapping of Two Interacting Films (Out-of-Phase) ---
% We simulate two film trajectories that are 180° out-of-phase.
t = linspace(0, 2*pi, 200);
film1 = 0.1 * sin(t);
film2 = 0.1 * sin(t + pi);  % Out-of-phase motion

figure;
plot(t, film1, 'b', 'LineWidth', 2); hold on;
plot(t, film2, 'r', 'LineWidth', 2);
xlabel('Time (rad)'); ylabel('Displacement (m)');
title('Simulated Out-of-Phase Flapping of Two Films');
legend('Film 1','Film 2'); grid on;

% --- Figure 10: Simulated Pressure Distribution due to the Bernoulli Effect ---
% We model a simple pressure drop where pressure decreases as velocity increases.
x_pos = linspace(0,1,200);
% Velocity increases sinusoidally along the gap between films.
velocity_profile = 8 + 2*sin(2*pi*x_pos);
pressure_profile = 101325 - 0.5*1.225*(velocity_profile.^2);

figure;
plot(x_pos, pressure_profile, 'k', 'LineWidth', 2);
xlabel('Normalized Position between Films');
ylabel('Pressure (Pa)');
title('Simulated Pressure Distribution (Bernoulli Effect)');
grid on;

%% ------------------------------------------------------------------------
%% End of Simulation
%% ------------------------------------------------------------------------
% This script demonstrates the interplay of fluid mechanics, film dynamics,
% and electrical outputs in the B-TENG device. The generated figures simulate
% key illustrations from the article, including the flapping mode, pressure variation,
% contact area evolution, and electric potential distribution.
