%% 
clc; clear; close all;

% full_HH_cable_propagation_and_node2_plot.


%% Simulation Parameters
dt     = 0.01;            % ms
T      = 40;              % ms
t      = 0:dt:T-dt;          % time vector (ms)
nSteps = numel(t);

%% Hodgkin–Huxley Node Parameters
Cm    = 1.0;              % μF/cm^2
gNa   = 120;              % mS/cm^2
gK    = 36;               % mS/cm^2
gL    = 0.3;              % mS/cm^2
ENa   = 50;               % mV
EK    = -77;              % mV
EL    = -54.4;            % mV
Vrest = -65;              % mV
Vecf = 0;                 % mV

% No extra Q10 scaling in this version
phi = 1;

%% Cable (Internode) Parameters (FitzGerald et al.)
radius_ax_um = .5;                    % μm
radius_ax_cm = radius_ax_um * 1e-4;     % cm

radius_mc_um = 1;                       % μm
radius_mc_cm = radius_mc_um * 1e-4;     % cm

% Intracellular resistivity
rho_i_ohm_m      = 1.1;                     % Ω·m
rho_i_ohm_cm     = rho_i_ohm_m * 100;       % Ω·cm
ri_ohm_per_cm    = rho_i_ohm_cm / (pi*radius_ax_cm^2);
ri_Mohm_per_mm   = ri_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

% Extracellular longitudinal resistivity
rho_e_long_ohm_m = 2.0;                     % Ω·m
rho_e_long_ohm_cm= rho_e_long_ohm_m * 100;   % Ω·cm
re_ohm_per_cm    = rho_e_long_ohm_cm / (pi*radius_ax_cm^2);
re_Mohm_per_mm   = re_ohm_per_cm * 1e-6 * 10;  % MΩ/mm

% Extracellular radial resistivity
rho_ecf_rad = 12.5; % Ω·m
rho_ecf_rad_cm = rho_ecf_rad * 100; % Ω·cm
re_rad_per_cm    = rho_ecf_rad_cm / (pi*(radius_mc_cm - radius_ax_cm)^2);

% Membrane resistance per mm (from leak conductance)
R_m_ohm_cm2      = 1/(gL*1e-3);             % Ω·cm^2
R_m_Mohm_cm2     = R_m_ohm_cm2 * 1e-6;       % MΩ·cm^2
area_per_mm_cm2  = 2*pi*radius_ax_cm * 0.1;     % cm^2/mm
rm_Mohm_per_mm   = R_m_Mohm_cm2 / area_per_mm_cm2; % MΩ/mm

% Space‐ and time‐constants for a “free” cable
lambda_mm = sqrt(rm_Mohm_per_mm / ri_Mohm_per_mm);    % mm
tau_m_ms  = R_m_Mohm_cm2 * Cm * 1e3;                  % ms

% Internode length & delay
L          = 0.05;        % mm (50 μm)
v_prop     = lambda_mm / tau_m_ms;   % mm/ms
delay_full = L / v_prop;             % ms

fprintf('ri = %.2f MΩ/mm, re = %.2f MΩ/mm, rm = %.2f MΩ/mm\n', ...
        ri_Mohm_per_mm, re_Mohm_per_mm, rm_Mohm_per_mm);
fprintf('τ_m = %.2f ms   λ = %.2f mm   v = %.2f mm/ms   delay = %.2f ms\n\n', ...
        tau_m_ms, lambda_mm, v_prop, delay_full);

%% Part 1: Node 1 AP Generation
I_stim1 = zeros(1,nSteps);
I_stim1(t>=5 & t<6) = 10;   % μA/cm^2 pulse at Node 1
[V1 Vout] = simulate_HH_node_RK4(I_stim1, t, dt, Vrest, Vecf, ri_ohm_per_cm, re_rad_per_cm, Cm, ...
     gNa, gK, gL, ENa, EK, EL, phi);

%% Part 2: Passive Decay Along Internode at L/3, L/2, 3L/4
distances = [L/3, L/2, 3*L/4];
V_decay   = zeros(numel(distances), nSteps);
for j = 1:numel(distances)
    d         = distances(j);
    idx_delay = round((d/v_prop)/dt);
    alpha     = exp(-d/lambda_mm);
    for i = 1:nSteps
        if i > idx_delay
            Vp = Vrest + (V1(i-idx_delay) - Vrest)*alpha;
        else
            Vp = Vrest;
        end
        V_decay(j,i) = Vp;
    end
end

%% Part 3: Passive Voltage at Node 2 → Extracellular Drive
idx_full       = round(delay_full/dt);
V_passive_full = Vrest*ones(1,nSteps);
for i = idx_full+1:nSteps
    V_passive_full(i) = Vrest + ...
        (V1(i-idx_full) - Vrest)*exp(-L/lambda_mm);
end
V_extra = V_passive_full - Vrest;  % extracellular voltage at Node 2

%% Part 4: Node 2 AP Regeneration
I_stim2 = gL * V_extra;      % only depolarizing component
I_stim2(I_stim2<0) = 0;
[V2, ~] = simulate_HH_node_RK4(I_stim2, t, dt, Vrest, Vecf, ri_ohm_per_cm, re_rad_per_cm, Cm, ...
     gNa, gK, gL, ENa, EK, EL, phi);

%% Part 5: Plot All Traces (Node1 → L/3 → L/2 → 3L/4 → Node2)
figure;
hold on;
plot(t, V1,           'k-',  'LineWidth',1.5);
plot(t, V_decay(1,:), 'b--', 'LineWidth',1.2);
plot(t, V_decay(2,:), 'g--', 'LineWidth',1.2);
plot(t, V_decay(3,:), 'm--', 'LineWidth',1.2);
plot(t, V2,           'r-',  'LineWidth',1.5);
plot(t, Vout,         'c-',  'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Node 1','L/3','L/2','3L/4','Node 2','Location','Best');
title('AP Propagation @37°C: Active & Passive + Regeneration');
grid on; hold off;

%% Part 6: Plot Intracellular vs Extracellular at Node 2
figure;
hold on;
plot(t, V2,      'k-',  'LineWidth',1.5);
plot(t, V_extra, 'b--', 'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Voltage at Node 2 (mV)');
legend('Intracellular (V_{intra})','Extracellular (V_{extra})','Location','Best');
title('Node 2: Regenerated AP vs. Extracellular Drive');
grid on; hold off;


%% ===== Local Functions =====

function [V Vout] = simulate_HH_node_RK4(I_stim, t, dt, Vrest, Vecf, ri_ohm_per_cm, re_rad_per_cm, Cm, gNa, gK, gL, ENa, EK, EL, phi)
    n    = numel(t);
    V    = zeros(1,n);
    Vout = zeros(1,n);
    m    = zeros(1,n);
    h    = zeros(1,n);
    n_g  = zeros(1,n);
    V(1) = Vrest;
    Vout(1) = Vecf;
    m(1) = alpha_m(Vrest)/(alpha_m(Vrest)+beta_m(Vrest));
    h(1) = alpha_h(Vrest)/(alpha_h(Vrest)+beta_h(Vrest));
    n_g(1)= alpha_n(Vrest)/(alpha_n(Vrest)+beta_n(Vrest));

    for i=1:n-1
        x = [V(i); m(i); h(i); n_g(i)];
        f = @(x) hh_derivatives(x, I_stim(i), Cm, gNa, gK, gL, ENa, EK, EL, Vrest, phi);
        k1= f(x);
        k2= f(x + 0.5*dt*k1);
        k3= f(x + 0.5*dt*k2);
        k4= f(x +    dt*k3);
        xn = x + (dt/6)*(k1+2*k2+2*k3+k4);
        V(i+1)   = xn(1);
        Vout(i+1) = Vout(i) + (V(i+1) - V(i)) * -ri_ohm_per_cm / (ri_ohm_per_cm + re_rad_per_cm);
        m(i+1)   = xn(2);
        h(i+1)   = xn(3);
        n_g(i+1) = xn(4);
    end
end

function dx = hh_derivatives(x, I, Cm, gNa, gK, gL, ENa, EK, EL, Vr, phi)
    V   = x(1); m=x(2); h=x(3); n_g=x(4);
    INa = gNa*(m^3*h)*(V-ENa);
    IK  = gK *(n_g^4)*(V-EK);
    IL  = gL      *(V-EL);
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


% NOISE STUFF

% Flicker Noise Inputs
Hooge = 5e-4;         % Hooge constant
%I = 1e-5;               % Current, A
%I_dyn = G_atotal * (max(Vout_temp) - min(Vout_temp));  % approximate current amplitude
Avo = 6.022e23;         % Avogadros number
M = 1;                  % Molarity, NaCl

% White Noise Inputs
k = 1.38e-23;           % Boltzmann constant
Temp = 298;                % Temperature 
rho = .089;            % Resistivity of NaCl at 1M 298K
L_pore = .5e-6;               % Length of Nanopore, m 
d_pore = radius_mc_um *2e-6;             % Diameter of Nanopore, m
d_axon = radius_ax_um * 2e-6;        % Diameter of axon, m
q = 1.602e-19;         % abs value, elementary charge of electron, C

% Dielectric Noise Inputs
sigma = 11.2;          % Conductivity of NaCl, S/m (assuming 1M and 298K)
plate_L = 1e-10;           %half the distance between capacitor plates, m assumed this
plate_area = 1e-7;       %capacitor plate area, m2, assumed, circular plate
E_SiN = 7;                  % dielectric constant of nanopore
E_Si = 11.7;                 %dielectric constant of coating, assuming coating on one side
Er = 78.3;              % permittivity of water at 298K
E0 = 8.854e-12;        % permitivity at vacuum

% Amplifier Noise Inputs
y = 2 / 3;             % FET constant
mu_n = .97e-8;         % ionic mobility of NaCl


% Frequency range
f1_low = -1;
f4_high = 3.8;

f1 = logspace(f1_low, 1, nSteps);  % From 0.1 Hz to 10^4 Hz (10 kHz)
f2 = logspace(1, 2.55, nSteps);   % from 10 kHz to 50 kHz
f3 = logspace(2.55, 3.41, nSteps);  % from 50 to 200 kHz
f4 = logspace(3.41, f4_high, nSteps); % from 200 to 1000 kHz
ftotal = logspace(f1_low, f4_high, nSteps);  % total Hz range

% Volume of the Nanopore, factoring in axon
r_pore = d_pore / 2;             % radius of nanopore
r_axon = d_axon / 2;     % radius of axon
Vpore = pi * r_pore ^ 2 * L_pore;   %volume of nanopore
V_axon = pi * r_axon ^ 2 * L_pore;      % volume of axon
Vpore_axon = Vpore - V_axon;      % vol of pore - axon

% Number of Charge Carriers 
n = 2 * M * 1000 * Avo;
Nc = max(n * Vpore_axon, 1e7);  % prevent flicker noise from vanishing

% Resistivity of Nanopore
R_total = rho * (((4 * L_pore) / (pi * d_pore ^ 2)) + (1 / d_pore));
R_atotal = rho * ((4 * L_pore) / pi * d_axon ^ 2);

% Conductance of Nanopore
G_total = 1 / R_total;
G_atotal = 1/ (R_total - R_atotal); 

%new current calc Idyn
I_dyn = G_atotal * (max(Vout) - min(Vout));  % approximate current amplitude

% Dielectric loss constant calc
w = 2 * pi * f3;            % angular freq
E = Er * E0;                % permitivity of soln
Dielect = sigma ./ (w * E);  % dielectric constant

% Non-ideal Capacitance, for SI_d
%plate_area = pi * (r_pore - r_axon)^2;  % <-- now varies with radius
C_m = (E0 * plate_area) / ((plate_L / Er) + (plate_L / E_SiN)); % capacitance between fluid and nanopore
C_s = (E0 * plate_area) / ((plate_L / E_Si) + (plate_L / E_SiN)); % capacitance between coating(if we have) and nanopore
C_n = C_m + C_s;

% % Transconductance calc
% g_m = sqrt(2 * mu_n * C_s * (d_pore / L_pore) * I);    % transconductance of amp
% 
% % Voltage Thermal Noise calc
% e_n = sqrt((8 * k * Temp * y) / g_m);      % input referred voltage noise of amp
e_n_input = 6e-7;   % assuming this is an input with our amplifier

% total capacitance for amplifier calc
C_cable = 2e-10 * .1;   %cable C, assuming 50pF/m and .1m cable length
C_amp = 5e-10;     %amp capacitance, assuming given, F
C_pore = (Er * E0 * pi * (r_pore - r_axon) ^ 2) / L_pore; % pore capacitance
C_atotal = C_cable + C_amp + C_pore; %total C for amp noise

% Spectral density calc for flicker noise
SI_f = (Hooge * I_dyn ^ 2) ./ (Nc * f1);

% Spectral density calc for white noise
SI_w = 4 * k * Temp * G_atotal + 2 * q * I_dyn;
SI_wv = SI_w * ones(size(f2));   % to make SI_w a vector

% Spectral density calc for dielectric noise
SI_d = 8 * pi * k * Temp * C_n .* f3;

% Spectral density calc for amplifier noise
SI_a = ((2 * pi * C_atotal * e_n_input) ^ 2) .* (f4 .^ 2);

% Bmax calc
Bmax = (I_dyn / (e_n_input * C_atotal)) ^ (2 / 3);  %integrating factor

% IRMS Calc
SI_total = SI_f + SI_w + SI_d + SI_a;   %total SI sum
%f = f1 + f2 + f3 + f4;    % all freq
Irms = sqrt(trapz(ftotal, SI_total));    % current of noise
Irms_pA = Irms * 1e12;
fprintf('Total Irms = %.3e A = %.2f pA\n', Irms, Irms_pA);


%% ===== IFFT-Based Noise Fluctuations — Now vs Frequency =====
Fs = 1/(.01e-3);    % Sampling frequency, Hz 
N = length(V_extra); % Length of the signal
df = Fs / N;        % Frequency resolution
f_full = (0:N/2) * df; % Frequency bins

% Interpolation domain (avoid zero frequency)
f_interp = linspace(1e-6, Fs/2, N/2 + 1);  % Small positive frequency value to avoid zero

% IFFT signal generation function
make_ifft_signal = @(Sxx_half) real(ifft( ...
    [sqrt(Sxx_half), fliplr(sqrt(Sxx_half(2:end-1)))] .* ...  % Scale the PSD and handle symmetry correctly
    [1, exp(1i * 2 * pi * rand(1, N/2-1)), 1, conj(exp(1i * 2 * pi * rand(1, N/2-1)))] ...
)) * sqrt(Fs);   % Scale by the sampling frequency

% ===== Define full-range PSDs for each noise type =====
SI_f_interp = (1 ./ f_interp); % Flicker noise ~ 1/f
SI_f_interp = SI_f_interp / max(SI_f_interp) * max(SI_f);  % Normalize to max of the flicker PSD

SI_w_interp = ones(size(f_interp)) * mean(SI_wv); % White noise (flat)

SI_d_interp = f_interp; % Dielectric noise ~ f
SI_d_interp = SI_d_interp / max(SI_d_interp) * max(SI_d);  % Normalize to max of the dielectric PSD

SI_a_interp = (f_interp.^2); % Amplifier noise ~ f^2
SI_a_interp = SI_a_interp / max(SI_a_interp) * max(SI_a);  % Normalize to max of the amplifier PSD

% ===== Generate synthetic fluctuations =====
noise_flicker = make_ifft_signal(SI_f_interp);
noise_white = make_ifft_signal(SI_w_interp);
noise_dielectric = make_ifft_signal(SI_d_interp);
noise_amplifier = make_ifft_signal(SI_a_interp);

% ===== FFT to go back to freq domain and get amplitudes =====
to_amp = @(x) abs(fft(x) / sqrt(Fs)); % Normalize the FFT to get the amplitude
amp_flicker = to_amp(noise_flicker);
amp_white = to_amp(noise_white);
amp_dielectric = to_amp(noise_dielectric);
amp_amplifier = to_amp(noise_amplifier);

% % Plot Noises
figure;
loglog(f1, SI_f, 'r', 'LineWidth', 2);
hold on
loglog(f2, SI_wv, 'b', 'LineWidth', 2);
loglog(f3, SI_d, 'g', 'LineWidth', 2);
loglog(f4, SI_a, 'm', 'LineWidth', 2);
loglog(ftotal, SI_total, 'k--', 'LineWidth', 2);
hold off
grid on;
xlabel('Frequency (Hz)');
ylabel('S_I (A^2/Hz)');
title('Current Noise Spectral Density vs Frequency');
legend('Flicker', 'White', 'Dielectric', 'Amplifier', 'Average');

% ===== Plot Noise Amplitudes vs Frequency =====
figure;
subplot(4,1,1);
title({'';'title of my plot'})
semilogx(f_full, 10^2 * amp_flicker(1:N/2+1), 'r');
xlabel('Frequency (Hz)'); ylabel('Noise_{f} (mV/Hz^{1/2})');
title('Flicker Noise Amplitude vs Frequency');

subplot(4,1,2);
title({'';'title of my plot'})
semilogx(f_full, amp_white(1:N/2+1), 'b');
xlabel('Frequency (Hz)'); ylabel('Noise_{w} (mV/Hz^{1/2})');
title('White Noise Amplitude vs Frequency');

subplot(4,1,3);
title({'';'title of my plot'})
semilogx(f_full, amp_dielectric(1:N/2+1), 'g');
xlabel('Frequency (Hz)'); ylabel('Noise_{d} (mV/Hz^{1/2})');
title('Dielectric Noise Amplitude vs Frequency');

subplot(4,1,4);
title({'';'title of my plot'})
semilogx(f_full, amp_amplifier(1:N/2+1), 'm');
xlabel('Frequency (Hz)'); ylabel('Noise_{a} (mV/Hz^{1/2})');
title('Amplifier Noise Amplitude vs Frequency');

sgtitle('Noise Fluctuations Across Frequency Spectrum');

% ===== Add all noise sources together =====
total_noise = noise_flicker + noise_white + noise_dielectric + noise_amplifier;

V_out_noise = ((10^10 *total_noise) .* Vout) + Vout;
V_out_amp = ((10^10 * noise_amplifier) .* Vout) + Vout;
V_out_d = ((10^10 * noise_dielectric) .* Vout) + Vout;
V_out_w = ((10^10 * noise_white) .* Vout) + Vout;
V_out_f = ((10^10 * noise_flicker) .* Vout) + Vout;

figure;
hold on;
plot(t, V1,           'k-',  'LineWidth',1.5);
plot(t, V_decay(1,:), 'b--', 'LineWidth',1.2);
plot(t, V_decay(2,:), 'g--', 'LineWidth',1.2);
plot(t, V_decay(3,:), 'm--', 'LineWidth',1.2);
plot(t, V2,           'r-',  'LineWidth',1.5);
plot(t, V_out_amp, 'm');
plot(t, V_out_d, 'g');
plot(t, V_out_w, 'b');
plot(t, V_out_f, 'r');
plot(t, Vout,         'c-',  'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
legend('Node 1','L/3','L/2','3L/4','Node 2', 'Amplifier Noise', 'Dielectric Noise', 'White Noise', 'Flicker Noise', 'Vout');
title('AP Propagation: with noise');
grid on; hold off;


%% === Clean SNR Optimization by Varying radius_mc_um ===
r_range = [radius_ax_um + 0.01, 50];  % μm range
r_values = linspace(r_range(1), r_range(2), 100);
snr_values = zeros(size(r_values));

% Constants (non-looped)
Fs = 1 / (0.01e-3);    % Sampling frequency (Hz)
f_interp = linspace(1, Fs/2, nSteps);  % Avoid DC
w_interp = 2 * pi * f_interp;

for i = 1:length(r_values)
    r_mc_um = r_values(i);
    radius_mc_cm = r_mc_um * 1e-4;
    d_pore = r_mc_um * 2e-6;
    r_pore = d_pore / 2;

    % --- Geometry Updates ---
    re_rad_per_cm = rho_ecf_rad_cm / (pi * (radius_mc_cm - radius_ax_cm)^2);
    C_pore = (Er * E0 * pi * (r_pore - r_axon)^2) / L_pore;
    C_atotal = C_cable + C_amp + C_pore;

    % --- Resistance & Conductance ---
    R_total = rho * (((4 * L_pore) / (pi * d_pore^2)) + (1 / d_pore));
    R_atotal = rho * ((4 * L_pore) / (pi * d_axon^2));  % constant
    %G_atotal = min(1 / max((R_total - R_atotal), eps), 1e-7);  % max 1 μS
    G_atotal = min(1 / max((R_total - R_atotal), eps), 1e-8);  % Clip to 1 μS

    % --- Nanopore Volume & Carrier Count ---
    Vpore = pi * r_pore^2 * L_pore;
    Vpore_axon = max(Vpore - V_axon, eps);
    Nc = max(n * Vpore_axon, 1e7);  % prevent flicker noise from vanishing

    % --- Run Vout ---
    [~, Vout_temp] = simulate_HH_node_RK4(I_stim1, t, dt, Vrest, Vecf, ri_ohm_per_cm, re_rad_per_cm, Cm, ...
         gNa, gK, gL, ENa, EK, EL, phi);

    % === SIGNAL (A): Estimate ΔI from Vout via conductance ===
    delta_V = max(Vout_temp) - min(Vout_temp);
    delta_I = G_atotal * delta_V * (radius_ax_um / r_mc_um);  % decay signal at large radii

    % === NOISE PSDs (A²/Hz) ===

    % Flicker Noise
    I_dyn = G_atotal * (max(Vout_temp) - min(Vout_temp));
    SI_f = (Hooge * I_dyn^2) ./ (Nc * f_interp);
    %SI_f = (Hooge * I^2) ./ (Nc * f_interp);

    % White Noise
    SI_w = 4 * k * Temp * G_atotal * ones(size(f_interp));  % flat PSD

    % Dielectric Noise
    C_m = (E0 * plate_area) / ((plate_L / Er) + (plate_L / E_SiN));
    C_s = (E0 * plate_area) / ((plate_L / E_Si) + (plate_L / E_SiN));
    C_n = C_m + C_s;
    SI_d = 8 * pi * k * Temp * C_n .* f_interp;  % A²/Hz

    % Amplifier Noise
    SI_a = ((2 * pi * C_atotal * e_n_input) .^ 2) .* (f_interp .^ 2) .* (r_mc_um / radius_ax_um)^2;


    % Total PSD
    SI_total_interp = SI_f + SI_w + SI_d + SI_a;

    % === Irms Calculation ===
    Irms_i = max(sqrt(trapz(f_interp, SI_total_interp)), 1e-12);  % 1 pA floor

    % === SNR Calculation ===
    snr_values(i) = delta_I / Irms_i;

    
    
end

% === Find and Plot Optimal Radius ===
[max_snr, idx_max] = max(snr_values);
optimal_radius_mc_um = r_values(idx_max);

fprintf('\n=== Optimized Pore Radius ===\n');
fprintf('radius_mc_um = %.3f μm\n', optimal_radius_mc_um);
fprintf('Max SNR = %.3f\n', max_snr);

figure;
plot(r_values, snr_values, 'b-', 'LineWidth', 2); hold on;
plot(optimal_radius_mc_um, max_snr, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Pore Radius (μm)');
ylabel('SNR');
title('SNR vs Pore Radius (radius\_mc\_um)');
grid on;
legend('SNR curve', sprintf('Optimum = %.2f μm', optimal_radius_mc_um));


% % Find optimal radius and plot
% [max_snr, idx_max] = max(snr_values);
% optimal_radius_mc_um = r_values(idx_max);
% 
% fprintf('\n=== Optimized Pore Radius ===\n');
% fprintf('radius_mc_um = %.3f μm\n', optimal_radius_mc_um);
% fprintf('Max SNR = %.3f\n', max_snr);
% 
% figure;
% plot(r_values, snr_values, 'b-', 'LineWidth', 2); hold on;
% plot(optimal_radius_mc_um, max_snr, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
% xlabel('Pore Radius (μm)');
% ylabel('SNR');
% title('SNR vs Pore Radius (radius\_mc\_um)');
% grid on;
% legend('SNR curve', sprintf('Optimum = %.2f μm', optimal_radius_mc_um));
