%%%%%%%%%%% dimensionamento radar
%   having as unkown parameter the mean power, assuming
%   theta azimut of 10 degrees, theta elevation 40 degrees, max distance 1
%   km, 60 rpm, 750MHz bandwidth and 9GHz operating frequency (siz stands
%   for sizing)
clear;
clc;
c = physconst('Lightspeed');
P_avg_vector_siz = zeros(10,1);
R_max_siz = (1000:100:1900);                            % [m] detection range


%%%%% Dimensioning radar with different r_max 
%
%
%
%
for n=1:1

    SNR_siz = 10^(20/10); % snr after integration (20dB)
    theta_siz_az = 36;    % [degrees]
    theta_siz_el = 60;   % [degrees]
    teta_az_siz = deg2rad(theta_siz_az);         % rad
    teta_el_siz = deg2rad(theta_siz_el);         % rad
    rpm_siz = 60;                                % [rpm]
    t_dwell_siz = theta_siz_az / (6 * rpm_siz);            % dwell time (Time spent on target during a scansion) [s]
    B_siz = 750e6; % bandwidth              
    lambda_siz = c/(9e9);     % [m] f0 9GHz

    t_d_max_siz = (2*R_max_siz(n))/c;              % Delay max related to Rmax (Rmax = (0.1*t_m*c)/2) [s]
    %t_m_siz  = 10*t_d_max_siz;                               % t_d_max = 0.1*t_m [s]
    t_m_siz = 2e-3;
    N_siz = t_dwell_siz/t_m_siz;
    teta_az_cov_siz = 2*pi;                            % [rad]
    Omega_s_siz = teta_az_cov_siz * teta_el_siz;       % [sterad]
    Omega_0_siz = teta_az_siz * teta_el_siz;           % [sterad]
    N_pos_siz = Omega_s_siz/Omega_0_siz;
    t_s_siz = (N_pos_siz * t_dwell_siz);               % [s]
    G_siz = (4*pi)/Omega_0_siz;                        %


    A_eq_siz = ((lambda_siz^2)*G_siz)/(4*pi);          % [m^2]
    RCS_siz = 0.01;                                    % [m^2]
    T_s_siz = 290;                                     % assumed equal to 290 K

    SNRmin_siz= SNR_siz/N_siz;                        %snr min before integration
    P_avg_siz = radar_mean_power(SNRmin_siz, R_max_siz(n), RCS_siz, t_s_siz, Omega_s_siz, T_s_siz, A_eq_siz, N_siz)
    P_avg_vector_siz(n) = P_avg_siz;
    S_siz = B_siz/t_m_siz;
    f_IF_max_siz = S_siz*((2*R_max_siz(n))/c);
end

stem(R_max_siz,P_avg_vector_siz)
xlabel('Max distance (m)')
ylabel('Mean power (W)')




%%%%%%%%%%%%%%%%%%% Functions
% Compute mean power knowing SNRmin (before integration), R
% max, RCS, scansion time  t_s, coverage solid angle Omega_s, sys
% temperature T_s, equivalent area of radar A_e and number of pulses sent 
% on  target in one scansion N 
function [P_avg] = radar_mean_power(SNR_min, R_max, RCS, t_s,  Omega_s, T_s, A_e, N)

k_b = physconst('Boltzmann');

P_avg_aperture = (SNR_min * N * 4 * pi * Omega_s * ((R_max)^4) * k_b * T_s)/(RCS * t_s);
P_avg = P_avg_aperture/A_e;


end
