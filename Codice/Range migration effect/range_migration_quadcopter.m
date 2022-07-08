clc;
clear;

%quadcopter drone
c = physconst('Lightspeed');
rho = 0.1;                  % ampiezza del moto di rotazione (raggio della rotazione) [m]
omega = 100*2*pi;            % velocità di rotazione [rad/s]     
f0 = 9e9;                % freq centrale [Hz]
lambda = c/f0;           % lunghezza d'onda [m]    
B = 750e6;               % banda del segnale utilizzato [Hz]
tc = 2e-3;               % durata del chirp           
mu = (B/tc);             % slope (verificare perchè è moltiplicato per 2PI)
tm = [0:0.0001:0.05];   % intervallo temporale
R0 = 0;                  % distanza centro di rotazione
N_r = 4;
N_b = 2;

%implementazione formula elaborata del libro
phi = atan((lambda*mu)/(c*omega)) -pi/2 + j*2*pi/N_b;    % elemento di fase

% formula come la presenta il chen (passaggi extra per riscriverla)
r_var_chen = ( R0 + rho*cos(omega*tm) + ((rho/(lambda*mu))*sqrt(((c^2)*(omega^2))+((lambda^2)*(mu^2))).*cos(omega.*tm + phi) )); % contributo sul range che varia nel tempo dovuto alla rotazione
%plot(tm,r_var_chen);
%ylabel('Formula chen')
% max quando arg coseno multiplo di k2pi
r_err_max_chen = ((rho/(lambda*mu))*sqrt(((c^2)*(omega^2))+(lambda^2)*(mu^2)))

% formula senza elaborazioni
r_err = ((c*rho*omega*sin(omega.*tm))/(lambda*mu));
r_err_max = max(r_err);

%   dimensione di una cella in range
cell = c/(2*B);

% max durata chirp per non incorrere in cell migration
Tchirp = ((lambda*B)*sqrt(3*cell^2 - rho^2))/(rho*c*omega)

%plot(tm,r_var_chen)

phi_0 = 2*pi.*rand(1,4);

% Distance from CM phase contribution

d_pos_rot= 202e-3;                % Max distance rotor from CM [m] 

delta = ((-d_pos_rot/2) + (d_pos_rot)*rand(2,1));         % Absolute value of position of each rotor (Only elevation angle is considered)
delta_final = [delta(1);-delta(1);delta(2);-delta(2)];      % Position of each rotor

% Received Signal from a Rotor with N_B blades
for r=1:N_r
    for j=0:N_b-1
        phi = atan((lambda*mu)/(c*omega)) -pi/2 + (j*2*pi/N_b) + phi_0(r);
        r_var_chen = ( R0 + delta_final(r) + rho*cos(omega*tm+phi_0(r)) + ((rho/(lambda*mu))*sqrt(((c^2)*(omega^2))+...
                        ((lambda^2)*(mu^2))).*cos(omega.*tm + phi) ));
        plot(tm,r_var_chen)
        hold on
        yline(R0+delta_final(r));
    end
end

axis([0 0.05 -2 2])
ylabel('Range variations [m]')
xlabel('Slow time [s]')

