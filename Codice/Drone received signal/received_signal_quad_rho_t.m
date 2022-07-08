clc;
clear;


c = physconst('lightspeed');            % Light Velocity [m/s]
f = 9e9;                                % Frequency [Hz]
lambda = c/f;                           % Wavelength [m] 
dr= 1;                                  % Range resolution [m]
B = c/(2*dr);                           % Bandwidth [Hz]
t_scan = 6;                             % Scan time [s]
SNR_db = 20;                            % Desidered SNR of received echo [dB]
SNR_linear = 10^(SNR_db/10);            % SNR conversion in linear
theta_el = 60;                          % Desidered coverage in elevation [deg]          
T_s = 290;                              % System temperature [Kelvin]
overs = 1;                             % Zero-Padding (used later in the FFT)
R_max = 100;


% 1. Drone's parameters

Omega = 116*2*pi;                        % Angular velocity of rotor [rad/s]
rho = linspace(0,0.05,9);                % Blade lenght [m]
N_b = 2;                                % Number of blades per rotor
N_r = 4;                                % Number of rotors
v_tip = Omega*rho(length(rho));         % Blade tip velocity  [m/s]  
T_BF = (2*pi)/(N_b*Omega);              % Period between two peak [s]
R0 = 0;                               % Distance of centre of rotation from radar [m]
f_max = 2*v_tip/lambda;                 % Maximum Doppler freqeuncy
rcs_drone = 0.01;                        % Average RCS [m^2]         
rcs_body = 0.25*rcs_drone;              % Body RCS Contribution [m^2]
rcs_blade = 0.75*rcs_drone;             % Blades RCS Contribution [m^2]
beta = deg2rad(15);                     % Elevation Angle [rad] 
fD_body = 0; 
N_bf = 10;

% Minimum T_chirp duration to have cell migration
T_chirp_max_migration = (lambda*B*dr)/(rho(length(rho))*Omega*c);

% Choose T_chirp in order to have cell migration
T_chirp = 1.3e-3;

% Radar parameters in function of T_chirp
mu = B/T_chirp;                                     % Slope [Hz/s]
f_IF_max = mu*((2*R_max)/c); %+ 2*v_tip/lambda;     % Maximum received frequency [Hz]
f_ADC = 2*f_IF_max;                                   % Sampling frequency  [Hz]         
T_ADC = 1/f_ADC;                                    % Sampling period   [s]



% Obesrvation time in order to visualize N_bf BF
t_dwell=0.1;                                  % Observation time [s]


% Time-Domain Discretization
ts = 0:T_ADC:t_dwell;
s_t = zeros(size(ts));



% Time-Domain signal
theta_0 = (2*pi).*rand(1,4);

% Taking into account different position of each rotors
d_pos_rot= 202e-3;                % Max distance rotor from CM [m] 

%delta = ((-d_pos_rot/2) + (d_pos_rot)*rand(2,1));         % Absolute value of position of each rotor (Only elevation angle is considered)
%delta_final = [delta(1);-delta(1);delta(2);-delta(2)];      % Position of each rotor

% Variations of velocity of rotation of each blade 
% Velocità di rotazione Gaussiana con sigma 10% e valore atteso omega_rev
%sigma = 0.1*Omega;
%Omega = Omega + randn(2)*sigma; 

% Equazione moto circolare uniforme
%a = rho*cos(Omega.*ts);
%b = -rho*Omega*sin(Omega.*ts);
%R_t = R0  +a+ b.*ts;
%v_r = -(rho*Omega)*sin(Omega.*ts+theta_0);

    
%Time domain signal
%s_t = exp( 1i.*((( -4*pi*mu/c).*mod(ts,T_chirp).*R_0) - ((4*pi/lambda)*R_0) + (4*pi*mu/c^2).*(R_0+v_r.*mod(ts,T_chirp)).^2)).* ...
       % exp( 1i*(-4*pi*mu/c*v_r.*(mod(ts,T_chirp).^2) + (4*pi/lambda*v_r).*mod(ts,T_chirp))) ;

for r=1:N_r
    for k=0:N_b-1
        for l=1:length(rho)
            R_0 = R0 + (rho(l)*cos(Omega.*ts + theta_0(r) + k*2*pi/N_b));
            s_t = s_t + (exp(+1i*(4*pi/lambda).*R_0) .* exp(+1i*(4*pi*mu/c).*mod(ts,T_chirp).*R_0).* exp(1i*(-4*pi*mu/c^2).*(R_0.^2)))...
                    .*sin(Omega.*ts + theta_0(r)+ k*2*pi/N_b);
        end
    end
end
%s_t =  exp(-1i*(4*pi/lambda).*R0) .* exp(-1i*(4*pi/lambda).*a) .* exp(-1i*(4*pi/lambda).*b.*mod(ts,T_chirp)) .* ...
%        exp(-1i*(4*pi*mu/c).*mod(ts,T_chirp).*R0) .* exp(-1i*(4*pi*mu/c).*mod(ts,T_chirp).*a) .* exp(-1i*(4*pi*mu/c).*b.*mod(ts,T_chirp).^2) .* ...
%        exp(1i*(4*pi*mu/c^2).*R0^2) .* exp(1i*(4*pi*mu/c^2).*a.^2);... .* exp(1i*(4*pi*mu/c^2).*b.^2.*mod(ts,T_chirp).^2);


% Body contribution
s_body = (sqrt(rcs_body)).*(exp(+1i*(4*pi/lambda).*R0) .* exp(+1i*(4*pi*mu/c).*mod(ts,T_chirp).*R0).* exp(1i*(-4*pi*mu/c^2).*(R0.^2)));
s_t = s_t + s_body;

% RCS Amplitude Contribution (considering one single rotor in helicopter
% case)
s_t = (sqrt(rcs_blade/N_bf)).*s_t;

%instf = instfreq(s_t,f_ADC);

% Maximum value of migration [m]
%r_err_max_chen = ((rho/(lambda*mu)).*sqrt(((c^2)*(Omega^2))+(lambda^2)*(mu^2)));


% Signal Power Computation
n = length(s_t);
S_f = fftshift(fft(s_t,n,2));
Ps = (1/n)*sum((abs(S_f)).^2);  % Power in linear scale
Ps_dB = 10*log10(Ps);           % Power in dB scale

% Add the AWGN noise
s_t = awgn(s_t,SNR_linear,'measured');


%plot(ts,10*log10((abs(s_t)).^2))

% Build the matrix in which put the signal

N = int64(T_chirp/T_ADC);   %campioni per ogni rampa
M = length(s_t);            %campioni totali

% riarrangio il vettore s_t in una matrice Nx(M/N)
% aggiungo M punti nulli al segnale per renderlo divisibile senza resto
% per N
while mod(M,N)>0
    M = M+1;
end

s_t(M) = 0;

matrix = reshape(s_t,N,int64(M/N));
matrix = transpose(matrix);
Fs = int64(f_ADC);

% faccio la fft per ogni riga della matrice (per ogni rampa)

y = (fftshift(fft(matrix,N,2)));   % 2 lungo le righe
freq = -Fs/2:(Fs/N):Fs/2;
y_dB = 10*log10(abs(y).^2);

figure
imagesc(y_dB)
title('Target Range Profile')

% Set x-axis parameters of Time-Range plot
clear xticklabels
clear xticks

tick = 5;
ticks = [1:tick*dr:N];
xticks(ticks)
x_labels = [-R_max:tick*dr:R_max];
xticklabels(x_labels)
xlabel('Range Bin')

% Set y-axis parameters of Time-Range plot
clear yticklabels
clear yticks

y_tick = 2;
ticks = [0:y_tick:M/N];
yticks(ticks);
y_labels = double(ticks)*T_chirp;
y_labels = round(y_labels,4);
yticklabels(y_labels)
ylabel('Slow Time [s]')

c = colorbar;
c.Label.String = 'Power [dB]';

% prendo la frequenza di ogni massimo di ogni fft
max_f_ramps = zeros(1,int64(M/N));
for i=1:int64(M/N)
    [~,max_f_ramps(i)] = max(abs(y(i,:)));
    %plot(freq,y_dB(i,:))
    %hold on
end

%hold off
%plot(freq(max_f_ramps(:)))



% prendo la frequenza di ogni massimo di ogni fft
max_f_ramps = zeros(1,int64(M/N));
for i=1:int64(M/N)
    [~,max_f_ramps(i)] = max(abs(y(i,:)));
    %plot(freq,y_dB(i,:))
    hold on
end

hold off

%plot(freq(max_f_ramps(:)))

% faccio ulteriore fft per misurare velocità ( analisi f-t caso senza
% migrazione)



N_fft = 12; % lunghezza fast time sotto matrici

n_pages = (M/N)/N_fft-1 ;

%creo la matrice 3D suddividendo in sottomatrici di lunghezza N_fft quella
%iniziale
matrix2 = zeros(N_fft*overs,N,n_pages);

% faccio la seconda fft della matrice iniziale, sotto matrice per sotto
% matrice

% faccio la prima pagina a mano
matrix2(:,:,1) = (fftshift(fft(y(1:N_fft,:),N_fft*overs,1)));

for j=1:n_pages-1
    matrix2(:,:,j+1) = (fftshift(fft(y(j*N_fft:((j+1)*N_fft)-1,:),N_fft*overs,1))); % 1 lungo le colonne
end


prova = matrix2(:,R0*dr+1,:);
prova = reshape(prova,N_fft*overs,n_pages);
figure
%imagesc((abs(prova).^2))

imagesc(10*log10(abs(prova).^2))
df = 1/(N_fft*T_chirp);
dt = N_fft*T_chirp;
title('Target Spectrogram')


% Set x-axis parameters of Frequecy-Time plot
clear xticklabels
clear xticks

ticks = [1:1:n_pages];
xticks(ticks)
x_labels = double(ticks)*N_fft*T_chirp;
x_labels = round(x_labels,4);
xticklabels(x_labels)
xlabel('Slow Time [s]')

% Set y-axis parameters of Time-Range plot
clear yticklabels
clear yticks


ticks = [0:N_fft:overs*N_fft];
yticks(ticks);
y_labels = linspace(-1/T_chirp,1/T_chirp,length(ticks));
yticklabels(-y_labels/2)
ylabel('Frequency [Hz]')

c = colorbar;
c.Label.String = 'Power [dB]';