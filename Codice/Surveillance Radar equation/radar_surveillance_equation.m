clc;
clear;

% Radar surveillance parameters with explained dimensions:

R_max = [];                     % maximum achievable distance [m]
P_avg = [];                     % Mean transmitted power [W]
A_eq = [];                      % Equivalent area of antenna [m^2]
RCS = [];                       % Radar Cross Section [m^2]          
N_pos = [];                     % Antenna positions in one scan
teta_az_cov = [];               % Azimuth angle covered by radar
teta_el_cov = [];               % Elevation angle covered by radar
Omega_s = [];                   % Solid angle of covered area equal to azimuth coverage * elevation coverage [sterad]
teta_az = [];                   % Bw in azimuth of beam [deg]
teta_el = [];                   % Bw in elevation of beam [deg]
Omega_0 = teta_az * teta_el;    % Solid angle covered by beam [sterad] 
SNR_min = [];                   % Minimum SNR to  detect a target with a given p_d and p_fa
k_b = physconst('Boltzmann');
T_s = [];                       % System temperature (assume T_0 = 290 K ?)
instrument_range = [];          % assumed maximum unambiguous distance reachable by the radar [m]
rpm = [];                       % velocity of a rotating antenna measured in rpm
t_dwell = teta_az/(6*rpm);      % time that beam spent on target 
c = physconst('Lightspeed');
t_m = (2*instrument_range)/c;   % Duration of sawtooth signal for FMCW [s]
N = t_dwell/t_m;                % Number of periods during a complete scan
t_s = N_pos*t_dwell;            % Scanning time (required to scan all N positions of antenna) [s]

% Assuming T_s (System temperature) as usual condition T_s = T_0 = 290 K
% Calculate SNR_min for IRIS radar
R_max_iris = 1000;                            % [m] detection range of rcs 0.03 m^2 (DJI Inspire drone)
P_avg_iris = 12;                              % [W]
teta_az_iris_deg = 6;
teta_el_iris_deg = 60;
teta_az_iris = deg2rad(teta_az_iris_deg);     % rad
teta_el_iris = deg2rad(teta_el_iris_deg);     % rad
rpm_iris = 30;                                % [rpm]
t_dwell_iris = teta_az_iris_deg / (6 * rpm_iris);    % dwell time (Time spent on target during a scansion) [s]

% considering FMCW radar the number of pulses can be computed considering
% the modulation period t_m, (each modulation period is related to a
% window in which there is the echo/es of traget/s).
% The modulation period t_m can be deducted considering the maximum measurable
% distance that is Rmax = c*t_d_max/2 and t_d_max (delay) is set always minor
% than t_m (t_d_max = 0.1*t_m) and max distance is 5 km for IRIS

instrument_range_iris = 4000;                       % [m]
B_iris = 9.65e9 - 8.9e9;                            % Bandwidth [Hz]
f_0_iris = 9.65e9 - B_iris/2;                       % Central frequency [Hz]
lambda_iris = physconst('Lightspeed')/f_0_iris;     % [m]

% metodo per ricavare il periodo di un chirp
% supponendo di avere fissato ritardo massimo 
% misurabile 1/30 della durata del chirp
 t_d_max = (2*instrument_range_iris)/c;              % Delay max related to Rmax (Rmax = (0.1*t_m*c)/2) [s]
 t_m_iris= 10*t_d_max;
 %30*t_d_max;                               % t_d_max = 0.0333*t_m [s]

% altro metodo per ricavare il periodo di un chirp
% fissado velocità massima v_max = lambda/(4*t_m)
% assumo v max misurabile 140 km/h
 %v_max = 140/3.6;                                    % [m/s]
 %t_m_iris = lambda_iris/(4*v_max);                   % [s]

N_iris = t_dwell_iris/t_m_iris;                     % N di rampe in dwell time
teta_az_cov_iris = 2*pi;                            % [deg]
Omega_s_iris = teta_az_cov_iris * teta_el_iris;     % [sterad]
Omega_0_iris = teta_az_iris * teta_el_iris;         % [sterad]
N_pos_iris = Omega_s_iris/Omega_0_iris;
t_s_iris = (N_pos_iris * t_dwell_iris);             % [s]
G_iris = (4*pi)/Omega_0_iris;                       


A_eq_iris = ((lambda_iris^2)*G_iris)/(4*pi);        % [m^2]
RCS_iris = 0.01;                                    % [m^2]
T_s_iris = 290;                                     % assumed equal to 290 K

SNR_min_sw1 = 10^(21.144/10);                       % p_d = 0.9 e pf = 10e-6 sw1-2 tgt
SNR_min_sw0 = 10^(13.2/10);                         % p_d = 0.9 e pf = 10e-6 sw0 tgt


% SNR_min considerando Ae = ((lambda_iris^2)*G_iris)/(4*pi)
SNR_min_iris_dB = surveillance_SNR_min(R_max_iris, P_avg_iris, A_eq_iris, RCS_iris, t_s_iris, Omega_s_iris, T_s_iris, N_iris)
SNR_min_iris = 10^(SNR_min_iris_dB/10);

%pocihè ho calcolato anche N nel caso di iris posso calcolare precisamente
%il SNR dato dall'integrazione degli N impulsi in una scansione
SNR_iris = SNR_min_iris * N_iris;
SNR_iris_dB = 10*log10(SNR_iris)

S_min_iris_dB = S_min_radar_equation(R_max_iris, A_eq_iris, P_avg_iris, RCS_iris, lambda_iris)  

% fattore di potenza per apertura d'antenna assunzione SNR_min sw1-2 
P_avg_iris = radar_mean_power(SNR_min_iris, R_max_iris, RCS_iris, t_s_iris, Omega_s_iris, T_s_iris, A_eq_iris, N_iris)

%calcolo lo slope del segnale di iris considerando che B=750MHz e B = S*t_m
S_iris = B_iris/t_m_iris
f_IF_max_iris = S_iris*((2*instrument_range_iris)/c)


%Utilizzando la funzione di matlab si ottiene lo stesso risultato (modo per
%verificare correttezza dei risultati)
snr_iris = radareqsearchsnr(R_max_iris, P_avg_iris * A_eq_iris, Omega_s_iris, t_s_iris, 'RCS', RCS_iris, 'Ts', T_s_iris, 'Loss', 0)

%modo piu preciso per calcolare angolo solido 
omega_s = solidangle([-180;180],[0;60])
range = linspace(1,50e3,10000);
snr = radareqsearchsnr(range,P_avg_iris * A_eq_iris,omega_s,t_s_iris,'RCS',RCS_iris,'Ts',T_s_iris,'Loss',0);
snr = radareqsearchsnr(range,P_avg_iris * A_eq_iris,omega_s,t_s_iris,'RCS',RCS_iris,'Ts',T_s_iris,'Loss',0);
plot(range*0.001,snr)
grid on
xlabel('Range (km)')
ylabel('SNR (dB)')
title('SNR vs Range')


%%%%%%%%%%%%%%%%%%% Functions
% using as input parameters the previous variables with specified dimension
% in order: P_avg, A_eq, RCS, t_s, N, Omega_s, SNR_min, T_s. Provide as
% output the maximum achievable distance.
function [R_max] = surveillance_distance(P_avg, A_eq, RCS, t_s, Omega_s, SNR, T_s)

k_b = physconst('Boltzmann');

R_max = (P_avg * A_eq * RCS * t_s)/(4 * pi * Omega_s * SNR * k_b * T_s);

% R_max value is at fourth power, to obtain R_max do the fourth root
R_max = nthroot(R_max,4);

end



% Provide SNR_min by knowing all the others parameters as input in order:
% R_max, P_avg, A_eq, RCS, t_s, N, Omega_s, T_s.
function [SNR_min] = surveillance_SNR_min(R_max, P_avg, A_eq, RCS, t_s, Omega_s, T_s, N)

k_b = physconst('Boltzmann');

SNR_min =  (P_avg * A_eq * RCS * t_s) / (4 * pi * Omega_s * N * ((R_max)^4) * k_b * T_s);

SNR_min = 10*log10(SNR_min); % dB

end



% Compute Smin (minimum received power) Knowing R_max, P_avg, RCS, lambda
function [S_min] = S_min_radar_equation(R_max, A_eq, P_avg, RCS, lambda)


S_min = (P_avg * (A_eq^2) * RCS)/((4*pi) * (lambda^2) * ((R_max)^4));
S_min = 10*log10(S_min); % dB

end


% Compute mean power knowing SNRmin (before integration), R
% max, RCS, scansion time  t_s, coverage solid angle Omega_s, sys
% temperature T_s, equivalent area of radar A_e and number of pulses sent 
% on  target in one scansion N 
function [P_avg] = radar_mean_power(SNR_min, R_max, RCS, t_s,  Omega_s, T_s, A_e, N)

k_b = physconst('Boltzmann');

P_avg_aperture = (SNR_min * N * 4 * pi * Omega_s * ((R_max)^4) * k_b * T_s)/(RCS * t_s);
P_avg = P_avg_aperture/A_e;


end
