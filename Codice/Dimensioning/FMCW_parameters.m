clc;
clear;

% prima verifica sulla freqeunza di campionamento
% Data la risoluzione in range desiderata DeltaR si calcola la banda

c = physconst('Lightspeed');
f = 9e9;                                       % central frequency [Hz]
lambda = c/f;                                  % wavelenght [m]
DeltaR = 0.6;                                    % range resolution [m]
B = c/(2*DeltaR);                              % Bandwidth
R_max = 1e3;                                   % max range [m]

T_chirp = [1e-6:0.1e-6:1e-3];                  % Sweep time interval [s]
T_chirp_final = [];
T_chirp_err = [];


for p=1:length(T_chirp)


                              
    mu = B/T_chirp(p);                                % slope [Hz/s]

    % sampling frequency to achieve max range
    fs_range = ( (R_max*4*B)/(c*T_chirp(p)) );


    %verifica che non si ha range cell migration nel caso helicopter drone
    %piu stringente (velocita punta della pala maggiore)

    rho = 0.6;                             % blade length [m]
    Omega = 25*2*pi;                       % rotation rate [rad/s]
    v_r_max = Omega*rho;                   % radial velocity of blade tip [m/s]


    %valore della migrazione
    RM = ((c*v_r_max)/(lambda*mu));

    if RM>DeltaR
        r = length(T_chirp_err);
        T_chirp_err(r+1) = T_chirp(p);
        text = ['Range cell migration presente per T_chirp =',num2str(T_chirp(p))];
        disp(text);
        continue;  
    else
        q = length(T_chirp_final);
        T_chirp_final(q+1) = T_chirp(p);
        continue;
    end
end

% sampling frequency to achieve max range and max blade tip velocity
%fs_range_vel = ( (R_max*4*B)/(c*T_chirp) ) + ((2*v_r_max)/lambda)

T_chirp_max = max(T_chirp_final)
T_chirp_min = min(T_chirp_final)

T_chirp_err_max = max(T_chirp_err)
T_chirp_err_min = min(T_chirp_err)
