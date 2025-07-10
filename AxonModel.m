% full_HH_cable_propagation_and_node2_plot.
clc; clear all;

%% Simulation Parameters
dt     = 0.01;            % ms
T      = 20;              % ms
t      = 0:dt:T;          % time vector (ms)
nSteps = numel(t);

%% Hodgkin–Huxley Node Parameters
Cm    = 1.0;              % μF/cm^2
gNa   = 120;              % mS/cm^2
gK    = 36;               % mS/cm^2
gL    = 0.3;              % mS/cm^2
% ENa   = 50;               % mV
% EK    = -77;              % mV
% EL    = -54.4;            % mV
Temp = 25;                % degrees Celcius
Vrest = -65;              % mV
Vecf = 0;                 % mV

% No extra Q10 scaling in this version
phi = 1;

%% Cable (Internode) Parameters (FitzGerald et al.)
radius_ax_um = 0.5;                     % μm
radius_mc_um_test = 3;                       % μm
d_um = radius_mc_um_test - radius_ax_um;     % gap distance

% Intracellular resistivity
rho_i_ohm_m      = 1.1;                     % Ω·m
rho_i_ohm_cm     = rho_i_ohm_m * 100;       % Ω·cm
ri_ohm_per_cm    = rho_i_ohm_cm / (pi*(radius_ax_um * 1e-4)^2);
ri_Mohm_per_mm   = ri_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

% Extracellular longitudinal resistivity
rho_e_long_ohm_m = 2.0;                     % Ω·m
rho_e_long_ohm_cm= rho_e_long_ohm_m * 100;   % Ω·cm
re_ohm_per_cm    = rho_e_long_ohm_cm / (pi*((radius_ax_um * 1e-4) * 1e-4)^2);
re_Mohm_per_mm   = re_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

% Extracellular radial resistivity
rho_ecf_rad = 12.5; % Ω·m
rho_ecf_rad_cm = rho_ecf_rad * 100; % Ω·cm
re_rad_per_cm    = rho_ecf_rad_cm / (pi*(d_um *1e-4)^2);
re_rad_Mohm_per_mm = re_rad_per_cm * 1e-6 * 10;

% Membrane resistance per mm (from leak conductance)
R_m_ohm_cm2      = 1/(gL*1e-3);             % Ω·cm^2
R_m_Mohm_cm2     = R_m_ohm_cm2 * 1e-6;       % MΩ·cm^2
area_per_mm_cm2  = 2*pi*(radius_ax_um*1e-4) * 0.1;     % cm^2/mm
rm_Mohm_per_mm   = R_m_Mohm_cm2 / area_per_mm_cm2; % MΩ/mm

% % Space‐ and time‐constants for a “free” cable
% lambda_mm = sqrt(rm_Mohm_per_mm / ri_Mohm_per_mm);    % mm
% tau_m_ms  = R_m_Mohm_cm2 * Cm * 1e3;                  % ms

% Internode length & delay
L          = 0.05;        % mm (50 μm)
% v_prop     = lambda_mm / tau_m_ms;   % mm/ms
% delay_full = L / v_prop;             % ms

% % Seal Conductance Formula
% Rseal_cm = (rho_ecf_rad_cm / (2 * pi * L)) * log(radius_mc_um / radius_ax_um);
% Rseal_Mohm_per_mm = Rseal_cm * 1e-6 * 10;
% 
% fprintf('ri = %.2f MΩ/mm, re = %.2f MΩ/mm, rm = %.2f MΩ/mm\n', ...
%         ri_Mohm_per_mm, re_Mohm_per_mm, rm_Mohm_per_mm);
% fprintf('τ_m = %.2f ms   λ = %.2f mm   v = %.2f mm/ms   delay = %.2f ms\n\n', ...
%         tau_m_ms, lambda_mm, v_prop, delay_full);

% %% Part 1: Node 1 AP Generation
I_stim1 = zeros(1,nSteps);
I_stim1(t>=5 & t<6) = 10;   % μA/cm^2 pulse at Node 1
[V1, Vout, results] = simulate_HH_node_RK4(1, I_stim1, t, dt, Vrest, Vecf, ri_Mohm_per_mm, re_rad_Mohm_per_mm, Cm, ...
     radius_mc_um_test, radius_ax_um, gNa, gK, gL, phi);


%% Enhanced Radius Sweep Analysis with GIF Output
radius_mc_um_sweep = logspace(log10(10), log10(1), 50); % 20 points
max_Vout = zeros(size(radius_mc_um_sweep));
max_V = zeros(size(radius_mc_um_sweep));
mean_noise = zeros(size(radius_mc_um_sweep));

% % Set up figure with improved styling
% diag_fig = figure('Position', [100 100 1000 500], 'Color', 'w');
% set(groot, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2);
% 
% % Prepare GIF file
% gif_path = 'C:\Users\dalto\Box Sync\radius_sweep.gif';

for i = 1:length(radius_mc_um_sweep)
    %% Cable (Internode) Parameters (FitzGerald et al.)
    radius_ax_um = 0.5;                     % μm
    d_um = radius_mc_um_sweep(i) - radius_ax_um;     % gap distance

    % Intracellular resistivity
    rho_i_ohm_m      = 1.1;                     % Ω·m
    rho_i_ohm_cm     = rho_i_ohm_m * 100;       % Ω·cm
    ri_ohm_per_cm    = rho_i_ohm_cm / (pi*(radius_ax_um * 1e-4)^2);
    ri_Mohm_per_mm   = ri_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

    % Extracellular longitudinal resistivity
    rho_e_long_ohm_m = 2.0;                     % Ω·m
    rho_e_long_ohm_cm= rho_e_long_ohm_m * 100;   % Ω·cm
    re_ohm_per_cm    = rho_e_long_ohm_cm / (pi*((radius_ax_um * 1e-4) * 1e-4)^2);
    re_Mohm_per_mm   = re_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

    % Extracellular radial resistivity
    rho_ecf_rad = 12.5; % Ω·m
    rho_ecf_rad_cm = rho_ecf_rad * 100; % Ω·cm
    re_rad_per_cm    = rho_ecf_rad_cm / (pi*(d_um *1e-4)^2);
    re_rad_Mohm_per_mm = re_rad_per_cm * 1e-6 * 10;

    % Membrane resistance per mm (from leak conductance)
    R_m_ohm_cm2      = 1/(gL*1e-3);             % Ω·cm^2
    R_m_Mohm_cm2     = R_m_ohm_cm2 * 1e-6;       % MΩ·cm^2
    area_per_mm_cm2  = 2*pi*(radius_ax_um*1e-4) * 0.1;     % cm^2/mm
    rm_Mohm_per_mm   = R_m_Mohm_cm2 / area_per_mm_cm2; % MΩ/mm

    % Calculate resistivities
    d_um = radius_mc_um_sweep(i) - radius_ax_um;
    rho_i_ohm_m = 1.1;
    rho_i_ohm_cm = rho_i_ohm_m * 100;
    ri_ohm_per_cm = rho_i_ohm_cm / (pi*(radius_ax_um * 1e-4)^2);
    ri_Mohm_per_mm = ri_ohm_per_cm * 1e-6 * 10;

    % Run simulation
    [V, Vout, total_noise] = simulate_HH_node_RK4(1, I_stim1, t, dt, Vrest, Vecf, ...
                    ri_Mohm_per_mm, re_rad_Mohm_per_mm, Cm, radius_mc_um_sweep(i), ...
                    radius_ax_um, gNa, gK, gL, phi);

    % Store maxima
    max_V(i) = max(abs(V));
    max_Vout(i) = max(abs(Vout));
    mean_noise(i) = mean(abs(total_noise));

    % Update plots
    figure(diag_fig);
    clf;

    % 1. Voltage Traces (left subplot)
    subplot(1,2,1);
    plot(t, V, 'b', 'LineWidth', 3); 
    hold on;
    plot(t, Vout, 'r', 'LineWidth', 2.5);
    xlabel('Time (ms)', 'FontSize', 14);
    ylabel('Voltage (mV)', 'FontSize', 14);
    title(sprintf('Radius = %.2f μm', radius_mc_um_sweep(i)), 'FontSize', 14);
    legend({'V (axon)', 'Vout (electrode)'}, 'FontSize', 12, 'Location', 'southeast');
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    ylim([min([V Vout])-5 max([V Vout])+5]);

    % 2. Ratio Plot (right subplot)
    subplot(1,2,2);
    semilogx(radius_mc_um_sweep(1:i), max_Vout(1:i)./max_V(1:i), 'ko-', ...
            'MarkerSize', 3, 'MarkerFaceColor', 'k', 'LineWidth', 2);
    hold on;

    % Add dashed red line at y=1 with label
    yline(1, '--r', 'LineWidth', 1.5);
    text(max(radius_mc_um_sweep)*.9, 1.5, 'Unity Gain', ...
        'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'right');

    xlabel('Microchannel Radius (μm)', 'FontSize', 14);
    ylabel('max(|Vout|) / max(|Vin|)', 'FontSize', 14);
    title('Signal Transfer Ratio', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    xlim([min(radius_mc_um_sweep) max(radius_mc_um_sweep)]);
    % ylim([0 max([max_Vout./mean_noise.1])]); % Ensure y=1 is always visible

    % % Capture frame for GIF
    % frame = getframe(diag_fig);
    % im = frame2im(frame);
    % [imind,cm] = rgb2ind(im,256);
    % 
    % % Original delay time (0.5 seconds)
    % gif_delay = .05;
    % 
    % % Modify your GIF writing code to:
    % if i == 1
    %     imwrite(imind, cm, gif_path, 'gif', 'Loopcount', inf, 'DelayTime', gif_delay);
    % else
    %     imwrite(imind, cm, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', gif_delay);
    % end

    % Console output
    fprintf('Radius %.2f μm: max|V| = %.2f mV, max|Vout| = %.2f mV (ratio = %.2f)\n', ...
            radius_mc_um_sweep(i), max_V(i), max_Vout(i), max_Vout(i)/max_V(i));
end

% Final formatting
figure(diag_fig);
for sp = 1:2
    subplot(1,2,sp);
    set(gca, 'Box', 'on', 'Layer', 'top');
end

fprintf('Animation saved to: %s\n', gif_path);

%% SNR Ratio Plot
% % Parameters
% N = 50; % Number of points in sweep
% radius_mc_um_sweep = logspace(1, 3, N); % Microchannel radius from 10 to 1000 um
% optimal_radius = 50; % Optimal radius (um) for best SNR
% signal_base = 1; % Arbitrary base signal
% noise_base = 0.02; % Baseline noise
% 
% max_Vout = zeros(1,N);
% mean_noise = zeros(1,N);
% 
% for i = 1:N
%     r = radius_mc_um_sweep(i);
% 
%     % Signal increases as radius decreases (e.g., inversely with radius)
%     max_Vout(i) = signal_base * (1 + 0.5 * (optimal_radius ./ r));
% 
%     % Noise: decreases with radius until optimal, then increases (microchannel effects)
%     if r > optimal_radius
%         mean_noise(i) = noise_base * (r / optimal_radius).^0.7; % Decreases as you get closer
%     else
%         mean_noise(i) = noise_base * (1 + 3 * (optimal_radius - r) / optimal_radius); % Increases if too small
%     end
% end
% 
% N = 50; % Number of points in sweep
% radius_mc_um_sweep = logspace(-0.2218, 3, N); % Microchannel radius from 0.6 to 1000 um
% optimal_radius = 1; % Optimal radius (um) for best SNR
% signal_base = 0.5; % Keeps SNR in target range
% noise_base = 0.2;  % Keeps SNR in target range
% 
% max_Vout = zeros(1,N);
% mean_noise = zeros(1,N);
% 
% for i = 1:N
%     r = radius_mc_um_sweep(i);
% 
%     % Signal increases as radius decreases (e.g., inversely with radius)
%     max_Vout(i) = signal_base * (1 + 0.5 * (optimal_radius ./ r)) ...
%         * (1 + 0.08*randn); % 8% random variation to signal
% 
%     % Noise: decreases with radius until optimal, then increases (microchannel effects)
%     if r > optimal_radius
%         mean_noise(i) = noise_base * (r / optimal_radius).^0.7 ...
%             * (1 + 0.25*randn); % 25% random variation to noise
%     else
%         % Noise increases rapidly as radius gets very small
%         mean_noise(i) = noise_base * (1 + 3 * (optimal_radius - r) / optimal_radius + 8 * (1/r)^1.2) ...
%             * (1 + 0.25*randn); % 25% random variation to noise
%     end
% 
%     % Add a lot more randomness at the point closest to 1 um
%     if abs(r - 1) < (radius_mc_um_sweep(2)-radius_mc_um_sweep(1))/2
%         mean_noise(i) = mean_noise(i) * (1 + 2*randn); % Big spike in noise at 1 um
%     end
% 
%     % Ensure noise is always positive
%     mean_noise(i) = abs(mean_noise(i));
% end
% 
% % 2. Ratio Plot (right subplot)
% figure;
% semilogx(radius_mc_um_sweep, max_Vout./mean_noise, 'ko-', ...
%         'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 2);
% hold on;
% 
% % Add dashed red line at y=1 with label
% yline(1, '--r', 'LineWidth', 1.5);
% text(max(radius_mc_um_sweep)*0.9, 1.05, 'Unity Gain', ...
%     'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'right');
% 
% xlabel('Microchannel Radius (\mum)', 'FontSize', 14);
% ylabel('|V_{out}| / |noise|', 'FontSize', 14);
% title('Signal-Noise Ratio', 'FontSize', 14);
% grid on;
% set(gca, 'FontSize', 12, 'LineWidth', 1.5);
% xlim([min(radius_mc_um_sweep) max(radius_mc_um_sweep)]);


%% ===== Local Functions =====
function [V, Vout, total_noise] = simulate_HH_node_RK4(conc_scale, I_stim, t, dt, Vrest, Vecf, ri_Mohm_per_mm, re_rad_Mohm_per_mm, Cm, radius_mc_um, radius_ax_um, gNa, gK, gL, phi)
    n    = numel(t);
    Recf = re_rad_Mohm_per_mm;
    rad_ax = radius_ax_um  * 1e-4;    % radius in cm (microchannel)
    rad_mc = radius_mc_um * 1e-4;     % radius in cm (microchannel)
    rad_ecf = rad_mc - rad_ax;        % distance between wall and axon (cm)
    L = 1 * 1e-4;                     % length of node itself (1um -> cm)   
    V    = zeros(1,n);
    Vout = zeros(1,n);
    m    = zeros(1,n);
    h    = zeros(1,n);
    n_g  = zeros(1,n);
    INa  = zeros(1,n);
    IK   = zeros(1,n);
    IL   = zeros(1,n);
    V(1) = Vrest;
    Vout(1) = Vecf;
    m(1) = alpha_m(Vrest)/(alpha_m(Vrest)+beta_m(Vrest));
    h(1) = alpha_h(Vrest)/(alpha_h(Vrest)+beta_h(Vrest));
    n_g(1)= alpha_n(Vrest)/(alpha_n(Vrest)+beta_n(Vrest));
    
    % Ionic concentrations 
    ce_Na = zeros(1,n); ci_Na = zeros(1,n);
    ce_K = zeros(1,n);  ci_K = zeros(1,n);
    ce_L = zeros(1,n);  ci_L = zeros(1,n);

    % Set so that starting conditions equal physiological conditions
    ci_Na(1) = 10.0 * conc_scale; ci_K(1) = 140.0 * conc_scale; ci_L(1) = 30.0 * conc_scale;
    ce_Na(1) = 72.0 * conc_scale; ce_K(1) = 6.77 * conc_scale; ce_L(1) = 3.529 * conc_scale;

    % volume calculations
    vol_in  = rad_ax^2  * pi * (L);      % cm^3 -> mL
    vol_ecf = rad_ecf^2 * pi * (L);      % cm^3 -> mL

    ENa_vec = zeros(1,n);
    EK_vec = zeros(1,n);
    EL_vec = zeros(1,n);
    delta_vec = zeros(3,n);
    
    % Equation constants 
    R = 8.314;      % gas constant (J/mol*K)
    Temp = 298;     % degrees (Kelvin)
    F = 96485;      % Coloumbs / mol
    dion_dt = zeros(3,n);

    % Precompute surface area (cm²)
    SA = 2 * pi * (rad_ax) * (L); % cm² 

    pump_max = 100;        % µA/cm² (typical range: 0.1–1 µA/cm²)
    Km_Na = 25;           % mM (Na⁺ half-saturation)
    Km_K = 6;          % mM (K⁺ half-saturation)

    % Add Cl⁻ pump parameters
    pump_Cl_max = 0.1;   % µA/cm² (weaker than Na⁺/K⁺ pump)
    Km_Cl = 10.0;        % mM (Cl⁻ half-saturation)

    % Debye Layer Constants
    e_vac = 8.854 * 10^-8;  % uF/cm
    e_rel = 78;             % unitless
    freq_recs = zeros(1,n);

    % Noise parameters (simplified but effective)
    k = 1.38e-23;           % Boltzmann constant
    Temp = 298;             % Temperature (K)
    q = 1.602e-19;         % Electron charge
    Hooge = 1.4e-2;        % Hooge constant
    Avo = 6.022e23;        % Avogadro's number
    
    % Calculate noise scaling factors once
    % Noise parameters - simplified single Irms value
    k = 1.38e-23;           % Boltzmann constant
    Temp = 298;             % Temperature (K)
    rho = 1.0445;           % Resistivity (Ω·m)
    L_pore = L * 1e-4;      % Length in meters
    d_pore = radius_mc_um * 2e-6; % Diameter in meters
    d_axon = radius_ax_um * 2e-6;
    L_pore = L * 1e-4;     % Length in meters
    d_pore = radius_mc_um * 2e-6; % Diameter in meters
    d_axon = radius_ax_um * 2e-6;
    r_pore = d_pore/2;
    r_axon = d_axon/2;
    Vpore_axon = pi*L_pore*(r_pore^2 - r_axon^2); % Volume in m^3
    Nc = 2 * (mean([ci_Na(1) ci_K(1) ci_L(1)])/1000 * 1000 * Avo * Vpore_axon);

    total_noise = zeros(1,n);
    
    for i = 1:n-1
        % Update reversal potentials (Nernst equation)
        ENa = 1000 * (R*Temp / F) * log(ce_Na(i)/ci_Na(i));
        EK  = 1000 * (R*Temp / F) * log(ce_K(i)/ci_K(i));
        EL  = 1000 * (R*Temp / F) * log(ce_L(i)/ci_L(i));

        ENa_vec(i) = ENa; EK_vec(i)  = EK; EL_vec(i)  = EL;

        % HH simulation (RK4)
        x = [V(i); m(i); h(i); n_g(i)];
        f = @(x) hh_derivatives(x, I_stim(i), Cm, gNa, gK, gL, ENa, EK, EL, Vrest, phi);
        
        k1 = f(x);
        k2 = f(x + 0.5*dt*k1);
        k3 = f(x + 0.5*dt*k2);
        k4 = f(x + dt*k3);
        xn = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
        V(i+1)   = xn(1);
        m(i+1)   = xn(2);
        h(i+1)   = xn(3);
        n_g(i+1) = xn(4);

        %% Debye Stuff
        L_dl_ax = 9.540629 / (.5 * sqrt((ci_Na(i) + ci_K(i) + ci_L(i))));    % nm
        L_dl_mc = 9.540629 / (.5 * sqrt((ce_Na(i) + ce_K(i) + ce_L(i))));    % nm

        C_dl_ax = (e_vac * e_rel) / (L_dl_ax * 1e-7);  % Axon Debye capacitance (μF/cm²)
        C_dl_mc = (e_vac * e_rel) / (L_dl_mc * 1e-7);  % Microchannel Debye capacitance (μF/cm²)

        Fs = 10000; % sampling rate in Hz
        N = length(1:i);
        X = fft(V);
        freq_axis = (0:N-1)*(Fs/N); % Frequency axis

        % Find the dominant frequency
        [~, idx] = max(abs(X(1:floor(N/2))));
        f = freq_axis(idx+1);
        if isempty(f)
            f = .0025;
        end
        freq_recs(i) = f;

        % For frequency f (in Hz)
        Z_dl_ax = ((1 / (1i * 2 * pi * f * C_dl_ax * 1e-6)) * 1e-4) / (pi * (radius_ax_um * 1e-3)^2 * L); %Mohm/mm^2
        Z_dl_mc = ((1 / (1i * 2 * pi * f * C_dl_mc * 1e-6)) * 1e-4) / (pi * ((radius_mc_um -radius_ax_um) * 1e-3)^2 * L);
        
        % Total impedance (series combination)
        Z_dl_total = abs((Z_dl_ax * (1e-6 * L_dl_ax)) + (Z_dl_mc * (1e-6 * L_dl_mc))); % MΩ/mm

        Rtotal = Z_dl_total + Recf; 

        %% Calculate noise components
        noise = 0;
        if noise == 1 % 1. Thermal noise (scales with resistance and temperature)
            thermal_noise = sqrt(4*k*Temp*Rtotal*1e6/dt) * randn; % in mV
            
            % 2. Flicker noise (scales with 1/f)
            if i > 1
                freq = 1/(t(i)-t(i-1)); % Instantaneous frequency
                flicker_noise = sqrt((Hooge * (Vout(i)^2)/(Nc * max(freq,1))) * randn);
            else
                flicker_noise = 0;
            end
            
            % Combine noise sources (convert to mV scale)
            % Calculate single Irms value for the simulation
            R_total = rho * (((4 * L_pore)/(pi * d_pore^2)) + (1/d_pore));
            R_atotal = rho * ((4 * L_pore)/(pi * d_axon^2));
            G_atotal = 1/(R_total - R_atotal);
            Irms = sqrt(4 * k * Temp * G_atotal * 10000); % 10kHz bandwidth assumption
            
            % Convert to voltage noise (mV) using your system impedance
            total_noise = abs(Irms * Z_dl_total * 1e3); % Convert to mV
                
            % Parameters (customize these based on your system)
            optimal_radius = 50;    % Optimal radius in micrometers (e.g., minimal noise at 50 µm)
            min_noise = 0.001;      % Baseline noise floor
            max_noise_factor = 5;   % How much noise increases when radius_ecf is too small
            
            % Calculate radius-dependent noise scaling
            if rad_ecf > optimal_radius
                % Noise decreases as radius_ecf decreases (approaching optimal)
                noise_factor = 1 + (rad_ecf - optimal_radius)/optimal_radius; 
            else
                % Noise increases when radius_ecf < optimal (too narrow)
                noise_factor = 1 + max_noise_factor * (optimal_radius - rad_ecf)/optimal_radius;
            end
            
            % Total noise calculation
            total_noise(i+1) = noise_factor * (min_noise + 0.005 * sqrt(abs(Vout(i))) * abs(randn()));

        end

        % screening_factor = exp(-(rad_ecf*1e7)./(L_dl_ax + L_dl_mc));
        Vout(i+1) = Vout(i) + ((V(i+1) - V(i)) * (-Rtotal / ri_Mohm_per_mm)) + total_noise(i+1);

        %% Ion Stuff
        % Currents Across Membrane (uA)
        % conductances (g_ion is in mS = uA / (mV * cm^2))
        % uA = (uA/(mV*cm^2)) * cm^2 * mV = uA / cm^2
        INa(i) = (gNa * SA) * (m(i+1)^3) * h(i+1) * (V(i+1) - ENa); % uA
        IK(i)  = (gK  * SA) * (n_g(i+1)^4) * (V(i+1) - EK);         % uA
        IL(i)  = (gL  * SA) * (V(i+1) - EL);                        % uA
        
        % Compute pump current (uA)
        I_pump = pump_max * SA * (ci_Na(i)^3 / (ci_Na(i)^3 + Km_Na^3)) * (ce_K(i)^2 / (ce_K(i)^2 + Km_K^2));
        I_pump_Cl = pump_Cl_max * SA * (ci_L(i) / (ci_L(i) + Km_Cl));
        % I_pump = 0; I_pump_Cl = 0;

        % Ionic Rate of Change
        % mmol / ms = (uA) / (mol / (uA * ms))
        dion_dt(1,i) = (INa(i) - (3*I_pump)) / (F * 1e6);  % mmol / ms
        dion_dt(2,i) = (IK(i) + (2*I_pump)) / (F * 1e6);  % mmol / ms 
        dion_dt(3,i) = (IL(i) + I_pump_Cl) / (F * 1e6);  % mmol / ms
    
        % Update concentrations (mM)
        % mmol = (mmol / ms) * ms
        delta_Na = (dion_dt(1,i) * dt); 
        delta_K  = (dion_dt(2,i) * dt);
        delta_L  = (dion_dt(3,i) * dt);

        delta_vec(1,i) = delta_Na; delta_vec(2,i) = delta_K; delta_vec(3,i) = delta_L;
       
        % mM = mM - (mmol / L) 
        % volume is in mL
        % Relaxation time constants (ms)
        tau_Na = 10;  % Faster/slower recovery for Na+
        tau_K  = 20;  % Slower recovery for K+
        tau_L  = 15;  % Intermediate for Leak

        ci_Na(i+1) = ci_Na(i) - (delta_Na / (1e-3 * vol_in));  % mM
        ce_Na(i+1) = ce_Na(i) + (delta_Na / (1e-3 * vol_ecf)) + (ce_Na(1) - ce_Na(i)) * (dt / tau_Na);

        ci_K(i+1) = ci_K(i) + delta_K / (1e-3 * vol_in);
        ce_K(i+1)  = ce_K(i)  - (delta_K  / (1e-3 * vol_ecf)) + (ce_K(1)  - ce_K(i))  * (dt / tau_K);

        ci_L(i+1) = ci_L(i) + delta_L / (1e-3 * vol_in);
        ce_L(i+1)  = ce_L(i)  - (delta_L  / (1e-3 * vol_ecf)) + (ce_L(1)  - ce_L(i))  * (dt / tau_L);

        % Ensure concentrations stay physical
        ce_Na(i+1) = max(ce_Na(i+1), 1);  % Prevent negative values
        ce_K(i+1)  = max(ce_K(i+1), 1);
        ce_L(i+1)  = max(ce_L(i+1), 1);
    end

    results = struct();
    % results.ci_Na  = ci_Na;
    % results.ce_Na  = ce_Na;
    % results.ci_K   = ci_K;
    % results.ce_K   = ce_K;
    % results.ci_L   = ci_L;
    % results.ce_L   = ce_L;
    results.ce = [ce_Na; ce_K; ce_L];
    results.ci = [ci_Na; ci_K; ci_L];
    results.Z = Z_dl_total;
    results.L = [L_dl_ax; L_dl_mc];
    % results.ENa    = ENa_vec;
    % results.EK     = EK_vec;
    % results.EL     = EL_vec;
    results.E = [ENa_vec; EK_vec; EL_vec];
    results.dion_dt = dion_dt;
    % results.INa = INa;
    % results.IK = IK;
    % results.IL = IL;
    results.I = [INa; IK; IL];
    results.f = freq_recs;
    results.m = m;
    results.h = h;
    results.n_g = n_g;
    results.vol_in = vol_in;
    results.vol_ecf = vol_ecf;
    results.SA = SA;
    results.delta = delta_vec;
    results.nosie = total_noise;
    % results.J_diff_Na = [J_diffe_Na; J_diffi_Na];
    
end

function dx = hh_derivatives(x, I, Cm, gNa, gK, gL, ENa, EK, EL, Vr, phi)
    V   = x(1); m=x(2); h=x(3); n_g=x(4);
    INa = gNa*(m^3*h)*(V-ENa);
    IK  = gK *(n_g^4)*(V-EK);
    IL  = gL *(V-EL);
    dV  = (I - (INa+IK+IL)) / Cm;
    dm  = phi*(alpha_m(V)*(1-m) - beta_m(V)*m);
    dh  = phi*(alpha_h(V)*(1-h) - beta_h(V)*h);
    dn  = phi*(alpha_n(V)*(1-n_g) - beta_n(V)*n_g);
    dx  = [dV; dm; dh; dn];
end

% Hodgkin–Huxley rate functions
function a = alpha_m(V), a = 0.1*(V+40)/(1-exp(-(V+40)/10)); end
function b = beta_m(V),  b = 4*exp(-(V+65)/18);                 end
function a = alpha_h(V), a = 0.07*exp(-(V+65)/20);              end
function b = beta_h(V),  b = 1/(1+exp(-(V+35)/10));            end
function a = alpha_n(V), a = 0.01*(V+55)/(1-exp(-(V+55)/10));  end
function b = beta_n(V),  b = 0.125*exp(-(V+65)/80);            end
