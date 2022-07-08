clc;
clear;


c = physconst('lightspeed');            % Light Velocity [m/s]
f = 9e9;                                % Frequency [Hz]
lambda = c/f;                           % Wavelength [m] 
dr = 1;                                 % Range resolution [m]
B = c/(2*dr);                           % Bandwidth [Hz]
t_scan = 6;                             % Scan time [s]
SNR_db = 20;                            % Desidered SNR of received echo [dB]
SNR_linear = 10^(SNR_db/10);            % SNR conversion in linear
theta_el = 60;                          % Desidered coverage in elevation [deg]          
T_s = 290;                              % System temperature [Kelvin]
overs = 1;                              % Zero-Padding (used later in the FFT)
R_max = 1000;                           % Identification range [m]


% 1. Drone's parameters

Omega = 40*2*pi;                        % Angular velocity of rotor [rad/s]
rho = linspace(0.042,0.325,36);             % Blade lenght [m]
N_b = 2;                                % Number of blades per rotor
N_r = 1;                                % Number of rotors
v_tip = Omega*rho(length(rho));         % Blade tip velocity  [m/s]  
T_BF = (2*pi)/(N_b*Omega);              % Period between two peak [s]
R0 = 0;                                 % Distance of centre of rotation from radar [m]
f_max = 2*v_tip/lambda;                 % Maximum Doppler freqeuncy
rcs_drone = 0.01;                       % Average RCS [m^2]         
rcs_body = 0.25*rcs_drone;              % Body RCS Contribution [m^2]
rcs_blade = 0.75*rcs_drone;             % Blades RCS Contribution [m^2]
beta = deg2rad(15);                     % Elevation Angle [rad] 
v_body = 0;                             % Velocity of boddy [m/s]

% Minimum T_chirp duration to have cell migration
% T_chirp_max_migration = (lambda*B*dr)/(rho*Omega*c);
T_chirp = 81.16e-6;                                  % Chirp duration [s]
mu = B/T_chirp;                                     % Slope [Hz/s]
f_IF_max = mu*((2*R_max)/c); %+ 2*v_tip/lambda;     % Maximum received frequency [Hz]
f_ADC = 2*f_IF_max;                                 % Sampling freqeuncy [Hz]
T_ADC = 1/f_ADC;

t_dwell= 5*T_BF;                                    % Dwell time [s]

% Time-Domain Discretization (domain of signals)
ts = 0:T_ADC:t_dwell;
s_t = zeros(size(ts));      % Empty vector for signal

% Time-Domain signal
theta_0 = (2*pi).*rand(1,4);

% Distance from CM phase contribution
d_pos_rot= 250e-3;                                      % Max distance rotor from CM [m] (DJI matrice 300)
delta = ((-d_pos_rot/2) + (d_pos_rot)*rand(2,1));       % Absolute value of position of each rotor (Only elevation angle is considered)
delta_final = [delta(1);-delta(1);delta(2);-delta(2)];  % Position of each rotor

% Variations of velocity of rotation of each rotor 
% Velocità di rotazione Gaussiana con sigma 10% e valore atteso omega_rev
sigma = 0.1*Omega;
Omega_var = Omega + randn(1,4)*sigma;

% For each rotor
for r=1:N_r
    % For each blade
    for k=0:N_b-1
        % For each point
        for l=1:length(rho) 
            % Range point equation
            R_0 = R0 + (rho(l)*cos(Omega.*ts+theta_0(r)+k*2*pi/N_b)) + v_body.*ts; %delta_final(r);  
            % Beaten signal model
            s_t = s_t + (exp(1i*(4*pi/lambda).*R_0) .* exp(1i*(4*pi*mu/c).*mod(ts,T_chirp).*R_0).* exp(-1i*(4*pi*mu/c^2).*(R_0.^2)));
              ... .*sin(Omega.*ts+theta_0(r)+k*2*pi/N_b);
            ...* exp(+1i*(4*pi/lambda).*mod(ts,T_chirp).*v_r);
        
        end
    end
end

% RCS Amplitude Contribution of blades
s_t = (sqrt(rcs_blade)).*s_t;

%  Body contribution
R_body = R0 + v_body.*ts;
s_body = (sqrt(rcs_body))*(exp(1i*(4*pi/lambda).*R_body) .* exp(1i*(4*pi*mu/c).*mod(ts,T_chirp).*R_body).* exp(-1i*(4*pi*mu/c^2).*(R_body.^2)));
s_t = s_t + s_body;

% Signal Power Computation
n = length(s_t);
S_f = fftshift(fft(s_t,n,2));
Ps = (1/n)*sum((abs(S_f)).^2);  % Power in linear scale
Ps_dB = 10*log10(Ps);           % Power in dB scale

% Add the AWGN noise
s_t = awgn(s_t,SNR_linear,'measured');

%plot(ts,10*log10((abs(s_t)).^2))

% Start FMCW processing

N = int64(T_chirp/T_ADC);   % Samples per each chirp
M = length(s_t);            % Total samples

% Rearrangement vector s_t in a matrix Nx(M/N)
% Add zeros at the end of s_t in order to make divisble M per N
while mod(M,N)>0
    M = M+1;
end

s_t(M) = 0;

matrix = reshape(s_t,N,int64(M/N));
matrix = transpose(matrix);
Fs = int64(f_ADC);

% Do FFT per each row of matrix (per each chirp)
y = (fftshift(fft(matrix,N,2)));   % 2 along the row
freq = 0:(Fs/N):(Fs);
y_dB = 10*log10(abs(y).^2);

figure
title('Target Range Profile')
imagesc(y_dB)

% Set x-axis parameters of Time-Range plot
clear xticklabels
clear xticks
tick = 100;
ticks = [1:tick*dr:N];
xticks(ticks)
x_labels = [-R_max:tick*dr:R_max];
xticklabels(x_labels)
xlabel('Range [m]')

% Set y-axis parameters of Time-Range plot
clear yticklabels
clear yticks
y_tick = 50;
ticks = [0:y_tick:M/N];
yticks(ticks);
y_labels = double(ticks)*T_chirp;
y_labels = round(y_labels,4);
yticklabels(y_labels)
ylabel('Slow Time [s]')

c = colorbar;
c.Label.String = 'Power [dB]';

% To verify the max freqeuncy of each ramp
% prendo la frequenza di ogni massimo di ogni fft
max_f_ramps = zeros(1,int64(M/N));
for i=1:int64(M/N)
    [~,max_f_ramps(i)] = max(abs(y(i,:)));
    %plot(freq,y_dB(i,:))
    hold on
end

hold off

%plot(freq(max_f_ramps(:)))

% Do the second FFT
N_fft = 28;             % length  of STFT windows
overs = 1;              % zero padding
n_pages = (M/N)/(N_fft/4)-4;

%Build 3D matrix by dividing the initial matrix in sub-matrix of Nfft length
matrix2 = zeros(N_fft*overs,N,n_pages);

% Weigth of each submartix (necessary in case of overlap among window
% during the shifting)
hamm = hamming(N_fft);

% Do the first page by hand
matrix2(:,:,1) = (fftshift(fft(y(1:N_fft,:).*hamm,N_fft*overs,1)));
l = N_fft/4+1;
m = N_fft + N_fft/4;
for j=1:n_pages-1
    
    matrix2(:,:,j+1) = (fftshift(fft(y(l:m,:).*hamm,N_fft*overs,1))); % 1: FFT  along the columns
    l = l+ N_fft/4;
    m = m+ N_fft/4;
end

% Plot spectrogram
spect = matrix2(:,(R0*dr)+1,:); % take the column in range matrix where is the target
spect = reshape(spect,N_fft*overs,n_pages);
figure

imagesc(10*log10(abs(spect).^2))
df = 1/(N_fft*T_chirp);
dt = N_fft*T_chirp/4;
title('Target Spectrogram')

% Set x-axis parameters of Frequecy-Time plot
clear xticklabels
clear xticks
ticks = [1:10:n_pages];
xticks(ticks)
x_labels = double(ticks)*dt;
x_labels = round(x_labels,4);
xticklabels(x_labels)
xlabel('Slow Time [s]')

% Set y-axis parameters of Time-Range plot
clear yticklabels
clear yticks
ticks = [1:N_fft/4:overs*N_fft+1]; 
yticks(ticks);
y_labels = linspace(1/T_chirp,-1/T_chirp,length(ticks));
yticklabels(y_labels/2)
ylabel('Frequency [Hz]')

c = colorbar;
c.Label.String = 'Power [dB]';