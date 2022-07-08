clc;
clear;


% Fixed parameters

c = physconst('Lightspeed');
R_max = 1e3;                            % Maximum range [m]


% 1. Drone's parameters

Omega = 2*pi*40;                        % Angular velocity of rotor [rad/s]
rho = 0.600;                              % Blade lenght [m]
N_b = 2;                                % Number of blades per rotor
N_r = 1;                                % Number of rotors
v_tip = Omega*rho;                      % Blade tip velocity  [m/s]  
T_migr = (2*pi)/(N_b*Omega);        % Period between two peak [s]
R0 = 0;                                 % Distance of centre of rotation from radar [m]
RCS = 0.03;                             % RCS of drone [m^2]

% 2-3. Input radar parameters

f = 9e9;                                % Central frequency [Hz]
lambda = c/f;                           % Wavelenght [m]
dr= 1;                                % Range resolution [m]
B = c/(2*dr);                           % Bandwidth [Hz]
t_scan = 2;                             % Scan time [s]
SNR_db = 20;                            % Desidered SNR of received echo [dB]
SNR_linear = 10^(SNR_db/10);            % SNR conversion in linear
theta_el = 60;                          % Desidered coverage in elevation [deg]          
T_s = 290;                              % System temperature [Kelvin]
L_db = 10;                              % Noise figure [dB]
L_linear = 10^(L_db/10);                % Noise figure [linear]
f_max = 2*v_tip/lambda;


% 4. Input rho-t plane desidered parameters

N_min = 1;                             % Number of range-cell desidered on peak value of migration
M = 1;                                 % Number of time-cell desidered between two peak
N_mig_peak = 5;                        % Number of peak migration desidered to see


% 5. Compute dr,dt e t_dwell. (dr is already given as input parameter,
% dt is the chirp duration)


% Compute dwell time in order to satisfy N_mig_peak

t_dwell_min = N_mig_peak*T_migr;       % Minimum dwell time [s]


% 6. Checks that there is a T_chirp that satisfies N_min and M 

T_chirp_max = (T_migr)/(M);                                                                 % Upper limit for T_chirp to satisfy M [s]
T_chirp_min = (lambda*N_min)/(2*rho*Omega);
%T_chirp_min = sqrt( ((N_min^2-(rho^2/dr^2))*(lambda^2*B^2*dr^2))/(rho^2*c^2*Omega^2) );    % Lower limit for T_chirp to satisfy N_min [s]


% If min and max value are incompatible then stop

if any(T_chirp_min > T_chirp_max)
    disp('No one value of T_chirp satisfy the conditions.') 
    return
else 
    disp('T_chirp interval that satisfy required resolution conditions exist')
 end

% 7. Fix T_chirp value by maximizing M and N_min

i = 0;

while i == 0
    
    % increase N_min and M by 1 at time
    N_min = N_min + 1;
    M = M + 1;

    % recompute limits with new N_min and M
    T_chirp_max = (T_migr)/(M);
    T_chirp_min = (lambda*N_min)/(2*rho*Omega);
    %T_chirp_min = sqrt( ((N_min^2-(rho^2/dr^2))*(lambda^2*B^2*dr^2))/(rho^2*c^2*Omega^2) );

    % check if the maximum has been reached
    if any(T_chirp_min > T_chirp_max)
        disp('Maximum value of T_chirp found.') 
        N_min = N_min - 1;
        M = M - 1;
        i = 1; % stop while cycle
    end

end

% Final limits on T_chirp values

T_chirp_max = (T_migr)/(M);          % [s]   
T_chirp_min = (lambda*N_min)/(2*rho*Omega);
%T_chirp_min = sqrt( ((N_min^2-(rho^2/dr^2))*(lambda^2*B^2*dr^2))/(rho^2*c^2*Omega^2) );     %[s]


% Check if T_chirp max staisfy the unanmbiguous radar range condition

k_range = 10;                   % Input parameter k, choosen in order to achieve T_chirp>k*t_delay_max
t_delay_max = (2*R_max/c);      % Maximum 2-way delay of echoes received from R_max [s]

if any(T_chirp_max < k_range*t_delay_max)
    disp('T_chirp_max = ', num2str(T_chirp_max),' is lower than k = ', num2str(k_range), ' times t_delay_max ')
    return;
else 
    disp('Choosen T_chirp value satisfy unanbiguous range condition')
end



% Fix T_chirp to its maximum value

T_chirp = T_chirp_max;             % Chirp duration [s]
mu = B/T_chirp;                    % Slope [Hz/s]

rho_migr = (c*rho*Omega)/(lambda*mu);
%rho_migr = ((rho/(lambda*mu))*sqrt(((c^2)*(Omega^2))+(lambda^2)*(mu^2)));     % Values of migration peaks [m]

N_migr_cell = rho_migr/dr;         % Number of range cells on peak of migration with choosed T_chirp
N_time_cell = T_migr/T_chirp;      % Number of time cells between two peak of migration


% Fix time resolution of range-slow time grid

dt = T_chirp;  % [s]


% 8. Plot of range variations over time on grid with computed
%    resolutions dr and dt

ts = [0:T_chirp:t_dwell_min];                      % Slow time interval [s]
r_var = zeros(length(ts),N_b);                     % Initialzied array of vector for each rotation points around R0
phi_0 = (2*pi).*rand(1,1);


% Suppose N_b = 2 points each opposite to the other
for k=0:N_b-1

     % Phase of range equation
     phi = atan((lambda*mu)/(c*Omega)) -pi/2 + phi_0;

     % Range variation formula
     r_var(:,k+1) = ( R0 + ((rho/(lambda*mu))*sqrt(((c^2)*(Omega^2))+((lambda^2)*(mu^2))).*cos(Omega.*ts + phi + k*2*pi/N_b) ));
    
     plot(ts,r_var(:,k+1))
     yticks(-R0-100:dr:R0+100)
     xticks(0:T_chirp:t_dwell_min)
     ylabel('Range variations [m]')
     xlabel('Seconds [s]')
     hold on
     grid on
end

% Plot fixed rotor contribution at R0 distance from radar

yline(R0)


% 9. Compute radar parameters obtained with T_chirp and T_dwell computed

rpm = 60/t_scan;                             % [Rip per minute]
theta_az = 6*rpm*t_dwell_min;                % Angle resolution [rad]
text = ['Theta azimut = ', num2str(theta_az),' Degrees'];
disp(text)
d = 65*(lambda/theta_az);                    % Antenna dimension [m]
text = ['Antenna dimension= ', num2str(d),' m'];
disp(text)


% Compute required tx power 

Omega_0 = deg2rad(theta_az) * deg2rad(theta_el);  % Solid angle of coverage [sterad]
G = (4*pi)/Omega_0;                               % Gain of antenna
A_eq = ((lambda^2)*G)/(4*pi);                     % Equivalent area of antenna [m^2]  
N = t_dwell_min/T_chirp;


P_tx = p_tx_req(R_max, SNR_linear, G, RCS, T_chirp, T_s, lambda,L_linear);
text = ['P_tx = ', num2str(P_tx),' W'];
disp(text)


% Compute required sampling frequency

f_IF_max = mu*((2*R_max)/c) + 2*v_tip/lambda;       % Maximum received frequency [Hz]
f_ADC = 2*f_IF_max;                                 % Computed sampling frequency [Hz]

text = ['Required slope is: ', num2str(mu*1e-9), ' GHz/s'];
disp(text)
text = ['Required sampling freqeuncy is: ', num2str(f_ADC*1e-6), ' MHz'];
disp(text)




% Functions:

% Provide P_tx necessary knowing all the other parameters from radar
% surveillance equation for a CW radar

function [P_tx] = p_tx_req(R_max, SNR, G, RCS, T_chirp, T_s, lambda,L)

k_b = physconst('Boltzmann');

P_tx =   ( SNR * (4*pi)^3 * ((R_max)^4) * k_b * T_s*L) / ( (G^2) * lambda^2 * RCS * T_chirp );


end











