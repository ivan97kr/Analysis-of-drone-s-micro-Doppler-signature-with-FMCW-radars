clc;
clear;


% Fixed parameters

c = physconst('Lightspeed');
k_b = physconst('Boltzmann');
R_max = 1e3;                            % Maximum range [m]

% 1. Input radar parameters

f = 9e9;                                % Central frequency [Hz]
lambda = c/f;                           % Wavelenght [m]
dr= 1;                                    % Range resolution [m]
B = c/(2*dr);                           % Bandwidth [Hz]
t_scan = 2;                             % Scan time [s]
SNR_db = 20;                            % Desidered SNR of received echo [dB]
SNR_linear = 10^(SNR_db/10);            % SNR conversion in linear
theta_el = 60;                          % Desidered coverage in elevation [deg]          
T_s = 290;                              % System temperature [Kelvin]
overs = 10;                             % Zero-Padding (used later in the FFT)


% 1. Drone's parameters

Omega = 2*pi*40;                        % Angular velocity of rotor [rad/s]
rho = 0.325;                              % Blade lenght [m]
N_b = 2;                                % Number of blades per rotor
N_r = 1;                                % Number of rotors
v_tip = Omega*rho;                      % Blade tip velocity  [m/s]  
T_BF = (2*pi)/(N_b*Omega);              % Period between two peak [s]
R0 = 0;                                 % Distance of centre of rotation from radar [m]
f_max = 2*v_tip/lambda;                 % Maximum Doppler freqeuncy
rcs_drone = 0.3;                       % Average RCS [m^2]         
rcs_body = 0.25*rcs_drone;              % Body RCS Contribution [m^2]
rcs_blade = 0.75*rcs_drone;             % Blades RCS Contribution [m^2]
beta = deg2rad(15);                     % Elevation Angle [rad] 
fD_body = 0;  

% 2. Radar constraint in order to see a drone at R_max and to correct
% measure v_tip

k_r = 10; 
T_chirp_min_range = k_r*(2*R_max/c);                        % Unambiguous range measurment 
T_chirp_max_migration = (lambda*B*dr)/(rho*Omega*c);        % To avoid cell migration effect
T_chirp_max_velocity = 1/(2*f_max);                         % Unambiguous velocity measurment

% Take minimum among the 2 max limits
T_chirp_max_radar = min(T_chirp_max_velocity,T_chirp_max_migration);

% Check that interval is consistent 
if any(T_chirp_min_range > T_chirp_max_radar)
    disp('No one value of T_chirp satisfy the conditions of radar constraints.') 
    return
end

% We have now an interval of T_chirp values that are compatible with a
% radar able to measure v_tip and detect a target at R_max


% 4. Input f-t plane desidered parameters

N_BF = 5;                               % Number of desired BF to see
n_t = 5;                                % Number of time-pixel desidered between two BF
n_f = 5;                                % Number of freq-pizel desidered between two BF


% 5. Compute df,dt e tdwell

t_dwell_min = N_BF*T_BF;       % Minimum dwell time [s]
dt_req = T_BF/n_t;             % Maximum time resolution [s]
df_req = f_max/n_f;            % Maximum frequency resolution [Hz]


% Compute the interval of existence of N_fft values (window of STFT) and
% integer cast

N_fft_min = 1/(T_chirp_max_radar*df_req);
N_fft_max = dt_req/T_chirp_min_range;
N_fft_min_index = int8(N_fft_min);
N_fft_max_index = int8(N_fft_max);

% Check that interval of N_fft is consistent 
if any(N_fft_max_index<N_fft_min_index)
    disp('No one value of N_fft is compatible with T_chirp first constraint and desired resolution.')
    return;
end

% I need only the even N values
if (not(mod(N_fft_min_index,2) == 0))
    N_fft_min_index = N_fft_min_index+1;
end

% Data matrix in which each row correspond to N_fft multiple of 2 values (e.g. 2,4,6...) and contain:
% [T_chirp, df, dt, P_tx, Slope, f_adc]
data_matrix = zeros(N_fft_max_index-N_fft_min_index,6);

for N_fft=N_fft_min:2:N_fft_max


    % check if N_fft value chosen is in the existence interval
    if any(N_fft< N_fft_min | N_fft>N_fft_max)
        disp('N_fft value chosen does not belongs to the existence interval')
        return;
    end

    % 6. Compute T_chirp limits coming from resolution requirements

    T_chirp_max_res = (T_BF)/(n_t*N_fft);
    T_chirp_min_res = n_f/(N_fft*f_max);

    % If min and max value are incompatible then stop

    if any(T_chirp_min_res > T_chirp_max_res)
        disp('No one value of T_chirp satisfy the conditions.') 
        return
    end

    % 7. Fix T_chirp value by maximizing n_t and n_f

    i = 0;
    n_t = 5;
    n_f = 5;

    while i == 0
    
        % increase N_min and M by 1 at time
        n_f = n_f + 1;
        n_t = n_t + 1;

        % recompute df e dt with new nf e nt values
        dt_req = T_BF/n_t;             % Maximum time resolution [s]
        df_req = f_max/n_f;            % Maximum frequency resolution [Hz]

        % recompute N_fft existence interval
        N_fft_min_opt = 1/(T_chirp_max_radar*df_req);
        N_fft_max_opt = dt_req/T_chirp_min_range;
    
        % recompute limits with new n_t and n_f
        T_chirp_max_res = (T_BF)/(n_t*N_fft);
        T_chirp_min_res = n_f/(N_fft*f_max);

        % check if the maximum has been reached
        if any(T_chirp_min_res > T_chirp_max_res | T_chirp_max_res<T_chirp_min_range)%N_fft< N_fft_min_opt | N_fft>N_fft_max_opt)
            disp('Maximum value of T_chirp found.') 
            n_f = n_f - 1;
            n_t = n_t - 1;
            i = 1; % stop while cycle
        end

    end

    % Final T_chirp interval by resolution constraint

    T_chirp_max_res = (T_BF)/(n_t*N_fft);          % [s]   
    T_chirp_min_res = n_f/(N_fft*f_max);


    % Choose T_chirp max value compatible with T_chirp radar constraint
    while (T_chirp_max_res>T_chirp_max_radar)
        T_chirp_max_res = T_chirp_max_res-1e-6;  % decrease by 1 microsec at time
    end



    % Check if T_chirp max staisfy the unanmbiguous radar range condition

    if any(T_chirp_max_res < T_chirp_min_range)
        disp('Max T_chirp does not satisfy unambiguous range condition') 
        return
    end



    T_chirp = T_chirp_max_res;             % Chirp duration [s]
    mu = B/T_chirp;                        % Slope [Hz/s]

    % Fix time and frequency resolution 

    df_final = 1/(N_fft*T_chirp);         % Frequency resolution with chirp duration chosen
    dt_final = N_fft*T_chirp;           % Time resolution with chirp duration chosen

    n_t_final = T_BF/dt_final;            % Number of time pixel between 2 BF
    n_f_final = f_max/df_final;           % Number of frequency pixel between 2 BF

    % Recompute N_fft final existence interval
    N_fft_min_final = 1/(T_chirp_max_radar*df_final);
    N_fft_max_final = dt_final/T_chirp_min_range;


    % 8. Compute radar parameters obtained with T_chirp and T_dwell computed

    rpm = 60*t_scan;                             % [Rip per minute]
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
    SNR_fft = SNR_linear/N_fft;

    P_tx = p_tx_req(R_max, SNR_fft, G, rcs_drone, T_chirp, T_s, lambda);
    text = ['P_tx = ', num2str(P_tx),' W'];
    disp(text)



    % Compute required sampling frequency

    f_IF_max = mu*((2*R_max)/c) + 2*v_tip/lambda;       % Maximum received frequency [Hz]
    f_ADC = 2*f_IF_max;                                 % Computed sampling frequency [Hz]

    text = ['Required slope is: ', num2str(mu*1e-9), ' GHz/s'];
    disp(text)
    text = ['Required sampling freqeuncy is: ', num2str(f_ADC*1e-6), ' MHz'];
    disp(text)

    % Save tuple of values in data_matrix for row N_fft-N_fft_min+1
    index = int8(N_fft-N_fft_min_index+1);
    data_matrix(index,:)=[T_chirp,df_final,dt_final,P_tx,mu,f_ADC];
end

% Remove zero rows
data_matrix( all(~data_matrix,2), : ) = [];

% Ask user to select a N_fft value multiple of 2, between the interval of
% existence of N_fft
N_fft_input = 0;
while (N_fft_input<N_fft_min_index || N_fft_input>N_fft_max_index || not(mod(N_fft_input,2)) == 0)
    text = ['Choose a N_fft value multiple of 2 between ',num2str(N_fft_min_index),' and ',num2str(N_fft_max_index),' \n'];
    N_fft_input = input(text);
    if (N_fft_input<N_fft_min_index || N_fft_input>N_fft_max_index)
        disp('Please respect the interval!')
    end
        % I need only the even N values
    if (not(mod(N_fft_input,2) == 0))
        disp('Please select N_fft multiple of 2!')
    end
end
N_fft = N_fft_input;
index_final = N_fft-N_fft_min_index+1-((N_fft-N_fft_min_index)/2);
values = data_matrix(index_final,:);
T_chirp = values(1);
df_final = values(2);
dt_final = values(3);
P_tx = values(4);
mu = values(5);
f_ADC = values(6);
text = ['With N_fft = ', num2str(N_fft),' the max value of T_chirp achievable is: ',num2str(T_chirp*1e6), ' microsec',...
            newline, 'The max resolutions are:',newline,'df = ',num2str(df_final),' Hz, dt = ',num2str(dt_final*1e3), ' ms', ...
            newline, 'The necessary tx power is:  ',num2str(P_tx),' W',newline, 'The slope is: ',num2str(mu*1e-9), ' GHz/s', ...
            newline, 'The sampling freq. is: ', num2str(f_ADC*1e-6), ' MHz'];
disp(text);

% 9. Plot of freqeuncy variations over time selecting a desired value of
% N_fft from data_matrix

ts = -(t_dwell_min)/2:T_chirp:(t_dwell_min)/2-T_chirp;                      % Slow time interval [s]

%% Received Signal
s_t = zeros(size(ts));                     % Initialzied array of vector for each rotation points around R0
phi_0 = (2*pi).*rand(1);

% Received Signal from a Rotor with N_B blades
for k=0:N_b-1
    phi = Omega*ts + phi_0 + k*2*pi/N_b;
    s_t = s_t + (lambda./(1i*4*pi*tan(phi))).*(1-exp(1i.*4*pi*rho/lambda*cos(beta)*sin(phi)));
end

% RCS amplitude contribution
s_t = (sqrt(rcs_blade)).*s_t;

% Doppler Contribution (Translational Motion of the Target)
s_t = exp(1i*2*pi*fD_body.*ts).*s_t;

% % Target Body Contribution
s_body = (sqrt(rcs_body)).*exp(1i*2*pi*fD_body.*ts);
s_t = s_t + s_body;

% Signal Power Computation
n = length(s_t);
S_f = fftshift(fft(s_t,n,2));
Ps = (1/n)*sum((abs(S_f)).^2);  % Power in linear scale
Ps_dB = 10*log10(Ps);           % Power in dB scale

% Add the AWGN noise
s_t = awgn(s_t,SNR_linear,'measured');

% Noise Power Computation
Pn_dB = Ps_dB - SNR_db;    % Power in dB scale
Pn = 10^(Pn_dB/10);     % Power in linear scale

%% Time-Domain Signal
figure
plot(ts,10*log10((abs(s_t)).^2))
title('Target Time-Domain Signal');
xlabel('Time (s)');
ylabel('|s(t)|^2 [dB]');

%% Spectrogram 
% Hamming weights
N_fft = int8(N_fft)
hamm = hamming(N_fft);

% Plot Spectrogram 
figure
spectrogram(s_t,hamm,N_fft/2,N_fft*overs,1/T_chirp,'yaxis','centered');
title('Target Spectrogram');
colormap("gray")


%% Spectrum 
n = length(s_t);
fd_ax = linspace(-1/(T_chirp*2),1/(T_chirp*2),n);
y = fftshift(fft(s_t,n,2));
y_dB = 10*log10(abs(y).^2);

figure
plot(fd_ax,y_dB)
title('Target Spectrum');
xlabel('Frequency (Hz)');
ylabel('|S(f)|^2 [dB]');


% Functions:

% Provide P_tx necessary knowing all the other parameters from radar
% surveillance equation for a CW radar

function [P_tx] = p_tx_req(R_max, SNR, G, RCS, T_chirp, T_s, lambda)

k_b = physconst('Boltzmann');

P_tx =   ( SNR * (4*pi)^3 * ((R_max)^4) * k_b * T_s) / ( (G^2) * lambda^2 * RCS * T_chirp );


end
