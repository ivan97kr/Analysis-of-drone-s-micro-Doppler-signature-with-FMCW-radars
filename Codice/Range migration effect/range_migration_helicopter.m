clc;
clear;

%quadcopter drone
c = physconst('Lightspeed');
rho = 0.6;                  % ampiezza del moto di rotazione (raggio della rotazione) [m]
omega = 25*2*pi;            % velocità di rotazione [rad/s]     
f0 = 9e9;                % freq centrale [Hz]
lambda = c/f0;           % lunghezza d'onda [m]    
B = 250e6;               % banda del segnale utilizzato [Hz]
tc = 1e-3;               % durata del chirp           
mu = (B/tc);             % slope (verificare perchè è moltiplicato per 2PI)
tm = [0:tc:0.1];   % intervallo temporale
R0 = 0;                  % distanza centro di rotazione

%implementazione formula elaborata del libro
phi = atan((lambda*mu)/(c*omega)) -pi/2;    % elemento di fase

% formula come la presenta il chen (passaggi extra per riscriverla)
r_var_chen = ( R0 + ((rho/(lambda*mu))*sqrt(((c^2)*(omega^2))+((lambda^2)*(mu^2))).*cos(omega.*tm + phi) )); % contributo sul range che varia nel tempo dovuto alla rotazione
%plot(tm,r_var_chen);
%ylabel('Formula chen')
% max quando arg coseno multiplo di k2pi
r_err_max_chen = ((rho/(lambda*mu))*sqrt(((c^2)*(omega^2))+(lambda^2)*(mu^2)));

% formula senza elaborazioni
r_err = ((c*rho*omega*sin(omega.*tm))/(lambda*mu));
r_err_max = max(r_err);

%   dimensione di una cella in range
cell = c/(2*B);

r = R0 + rho*cos(omega.*tm) + (((c*rho*omega)/(lambda*mu))*sin(omega.*tm));
r2 = R0 + + rho*cos((omega.*tm)+pi) + (((c*rho*omega)/(lambda*mu))*sin((omega.*tm)+pi));
r3 = R0;

% max durata chirp per non incorrere in cell migration
Tchirp = (lambda*B*cell)/(rho*c*omega);

plot(tm,r)
hold on
plot(tm,r2)
hold on
yline(r3)
%axis([0 100e-3 -2 2])
%ylim([996 1004])
%yticks(994:0.6:1005)
%xticks(0:tc:0.1)
ylabel('Range variations [m]')
xlabel('Seconds [s]')
%grid on
