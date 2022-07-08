clear;
clc;
close all;

%% System parameters
c = 3e8;       % Light Velocity [m/s]
f = 9e9;       % Frequency [Hz]
lambda = c/f;  % Wavelength [m] 
overs = 10;    % Zero-Padding (used later in the FFT)
SNR = 100;     % Signal-to-Noise Ratio [linear]

% Radar Waveform Parameters
N = 64;         % Number of Pulses in the FFT Window
M = 3000;      % Total Number of Pulses on Target
PRF = 1/(30e-6);    % Pulse Repetition Frequency [Hz]
PRT = 1/PRF;    % Pulse Repetition Time [s]
tmax = M/PRF;   % Time on Target [s]

% Resolutions
delta_t = N/(2*PRF);   % Time Resolution [s]
delta_f = PRF/N;       % Frequency Resolution [Hz]

% Time-Domain Discretization
slow_time_tot = -(tmax)/2:PRT:(tmax)/2-PRT;                               


%% Helicopter-Drone Features
N_R = 1;            % Number of Rotors 
N_B = 2;            % Number of Blades 

beta = deg2rad(15);      % Elevation Angle [rad] 
theta = pi/2 - beta;     % Aspect Angle [rad]

L = 0.6;                 % Blade Length [m] 
W = 0.1;                 % Blade Width [m]
H = 0.03;                % Blade Height [m]

Omega_rad = 2*pi*25;              % Rotation Rate of the Rotor [rad/s]
Omega_rev = Omega_rad/(2*pi);     % Rotation Rate of the Rotor [rev/s]
V_tip = 2*pi*Omega_rev*L;         % Blade Tip Velocity [m/s]

phi_0 = (2*pi).*rand(1);          % Initial Phase of the Rotor

rcs_drone = 0.12;                 % Average RCS [m^2]         
rcs_body = 0.25*rcs_drone;        % Body RCS Contribution [m^2]
rcs_blade = 0.75*rcs_drone;       % Blades RCS Contribution [m^2]

fD_body = 0;                      % Doppler Frequency of the Body [Hz]


% Variations of velocity of rotation of each blade 
% Velocità di rotazione Gaussiana con sigma 10% e valore atteso omega_rev
sigma = 0.1*Omega_rev;
%Omega_rev = Omega_rev + rand(1)*sigma; 

%% Received Signal
s_t = zeros(size(slow_time_tot));

% Received Signal from a Rotor with N_B blades
for k=0:N_B-1
    phi = Omega_rev*2*pi*slow_time_tot + phi_0 + k*2*pi/N_B;
    s_t = s_t + (lambda./(1i*4*pi*tan(phi))).*(1-exp(1i.*4*pi*L/lambda*cos(beta)*sin(phi)));
end


% RCS Amplitude Contribution
% s_t = ((sqrt(4*pi/((lambda)^2)))*(H*L)*sin(theta) + (sqrt(4*pi/((lambda)^2)))*(W*L)*cos(theta)).*s_t;
s_t = (sqrt(rcs_blade)).*s_t;

% Doppler Contribution (Translational Motion of the Target)
s_t = exp(1i*2*pi*fD_body.*slow_time_tot).*s_t;

% % Target Body Contribution
s_body = (sqrt(rcs_body)).*exp(1i*2*pi*fD_body.*slow_time_tot);
s_t = s_t + s_body;

% Signal Power Computation
n = length(s_t);
S_f = fftshift(fft(s_t,n,2));
Ps = (1/n)*sum((abs(S_f)).^2);  % Power in linear scale
Ps_dB = 10*log10(Ps);           % Power in dB scale

% Add the AWGN noise
s_t = awgn(s_t,100,'measured');

% Noise Power Computation
Pn_dB = Ps_dB - SNR;    % Power in dB scale
Pn = 10^(Pn_dB/10);     % Power in linear scale



%% Time-Domain Signal
figure
plot(slow_time_tot,10*log10((abs(s_t)).^2))
title('Target Time-Domain Signal');
xlabel('Time (s)');
ylabel('|s(t)|^2 [dB]');


%% Spectrogram 
% Hamming weights
hamm = hamming(N);

% Plot Spectrogram 
figure
spectrogram(s_t,hamm,N/2,N*overs,PRF,'yaxis','centered');
title('Target Spectrogram');
colormap("gray")


%% Spectrum 
n = length(s_t);
fd_ax = linspace(-PRF/2,PRF/2,n);
y = fftshift(fft(s_t,n,2));
y_dB = 10*log10(abs(y).^2);

figure
plot(fd_ax,y_dB)
title('Target Spectrum');
xlabel('Frequency (Hz)');
ylabel('|S(f)|^2 [dB]');

