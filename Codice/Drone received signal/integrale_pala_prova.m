clc;
clear;


c = physconst('lightspeed');            % Light Velocity [m/s]
f = 9e9;                                % Frequency [Hz]
lambda = c/f;                           % Wavelength [m] 
dr= 1;                                  % Range resolution [m]
B = c/(2*dr);                           % Bandwidth [Hz]
t_scan = 2;                             % Scan time [s]
SNR_db = 20;                            % Desidered SNR of received echo [dB]
SNR_linear = 10^(SNR_db/10);            % SNR conversion in linear
theta_el = 60;                          % Desidered coverage in elevation [deg]          
T_s = 290;                              % System temperature [Kelvin]
overs = 10;                             % Zero-Padding (used later in the FFT)
R_max = 1000;


% 1. Drone's parameters

Omega = 25*2*pi;                        % Angular velocity of rotor [rad/s]
rho = 0.6;                              % Blade lenght [m]
N_b = 1;                                % Number of blades per rotor
N_r = 1;                                % Number of rotors
v_tip = Omega*0.6;                      % Blade tip velocity  [m/s]  
T_BF = (2*pi)/(N_b*Omega);              % Period between two peak [s]
R0 = 500;                              % Distance of centre of rotation from radar [m]
f_max = 2*v_tip/lambda;                 % Maximum Doppler freqeuncy
rcs_drone = 0.3;                        % Average RCS [m^2]         
rcs_body = rcs_drone;                   % Body RCS Contribution [m^2]
rcs_blade = 0.75*rcs_drone;             % Blades RCS Contribution [m^2]
beta = deg2rad(15);                     % Elevation Angle [rad] 
fD_body = 0;  

% Minimum T_chirp duration to have cell migration
% T_chirp_max_migration = (lambda*B*dr)/(rho*Omega*c);
T_chirp = 68.02e-6;

mu = B/T_chirp;                                     % Slope [Hz/s]
f_IF_max = mu*((2*R_max)/c); %+ 2*v_tip/lambda;       % Maximum received frequency [Hz]
f_ADC = 2*f_IF_max;     
T_ADC = 1/f_ADC;

t_dwell=5*T_BF;
% Time-Domain Discretization
ts = 0:T_ADC:t_dwell;
s_t = zeros(size(ts));


for k=0:N_b
    a_var = cos(Omega*ts + k*2*pi/N_b)/lambda;
    b_var = (mu/c)*cos(Omega*ts + k*2*pi/N_b).*mod(ts,T_chirp);
    c_var = (mu*2*R0*cos(Omega*ts + k*2*pi/N_b)/c^2);
    d_var = (mu*cos(Omega*ts + k*2*pi/N_b).^2)/c^2;


    s_t = s_t + (exp(1i*(4*pi/lambda).*R0) .* exp(1i*(4*pi*mu/c).*mod(ts,T_chirp).*R0).* exp(-1i*(4*pi*mu/c^2).*(R0.^2))).*...
        (((-1)^(3/4)).*exp(pi*1i*((a_var+b_var-c_var).^2/d_var)) .* ((erfz( ((-1^(1/4))*sqrt(1i)*(a_var+b_var-c_var-2*d_var*rho.*cos(Omega*ts+k*2*pi/N_b)))/sqrt(d_var))) - ...
        erfz( (sqrt(pi)*(-1^(1/4))*(a_var+b_var-c_var))/sqrt(d_var))  ) ./ sqrt(d_var));
end


s_body = (exp(1i*(4*pi/lambda).*R0) .* exp(1i*(4*pi*mu/c).*mod(ts,T_chirp).*R0).* exp(-1i*(4*pi*mu/c^2).*(R0.^2))); 
s_t = s_t + s_body;



N = int64(T_chirp/T_ADC);   %campioni per ogni rampa
M = length(s_t);            %campioni totali

% riarrangio i il vettore s_t in una matrice Nx(M/N)

while mod(M,N)>0
    M = M+1;
end

s_t(M) = 0;

matrix = reshape(s_t,N,int64(M/N));
matrix = transpose(matrix);
Fs = int64(f_ADC);


% faccio la fft per ogni riga della matrice (per ogni rampa)

y = (fftshift(fft(matrix,N,2)));   % 2 lungo le righe
freq = 0:(Fs/(2*N)):(Fs/2);
y_dB = 10*log10(abs(y).^2);
figure
imagesc(y_dB)

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



N_fft = 42; % lunghezza fast time sotto matrici

n_pages = (M/N)/N_fft;

%creo la matrice 3D suddividendo in sottomatrici di lunghezza N_fft quella
%iniziale
matrix2 = zeros(N_fft*2,N,n_pages);

% faccio la seconda fft della matrice iniziale, sotto matrice per sotto
% matrice

% faccio la prima pagina a mano
matrix2(:,:,1) = (fftshift(fft(y(1:N_fft,:),N_fft*2,1)));

for j=1:n_pages-1
    matrix2(:,:,j+1) = (fftshift(fft(y(j*N_fft:((j+1)*N_fft)-1,:),N_fft*2,1))); % 1 lungo le colonne
end


prova = matrix2(:,501,:);
prova = reshape(prova,N_fft*2,(M/N)/N_fft);
figure
%imagesc((abs(prova).^2))

imagesc(10*log10(abs(prova).^2))


