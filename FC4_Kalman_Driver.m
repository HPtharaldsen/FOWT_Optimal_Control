% FC4 - SS_Model_Driver
close all; clear all; clc

%% Set Top-Level Linear directory
linear_dir = 'C:\Umaine Google Sync\GitHub\FOWT_Optimal_Control\Models\FOCAL_C4\Linear_Files';

%% Load in Simulation or Test Data
load('Test_Results.mat','test_results');
load('OpenFAST_Results.mat','sim_results');

%% Prepare Time Vector
% Extract time vectors
test_time = test_results.Time;
sim_time = sim_results.Time;

% Smooth time vector
time = linspace(min(test_time),max(test_time),length(test_time))';

% Prepare to shift wind input
dt = max(test_time)/length(test_time);
causality_shift_index = floor(29.95/dt);

%% Prepare Wave Input
eta = test_results.Wave1Elev;

%% Prepare wind input
% Wind Case
wind_case = 3;

% Load wind file
wind_file_path = 'C:\Umaine Google Sync\GitHub\FOWT_Optimal_Control\Wind_Files';
wind = getFOCALWindVector(wind_file_path,wind_case,test_time);

% Time-shift wind to account for hydro causalization time
wind = [zeros(causality_shift_index,1);
        wind(1:end-causality_shift_index)];

%% Prepare Control Input Values
% Blade pitch command (collective)
c_pitch = test_results.pitch1Position*(pi/180);

% Ind. Pitch
idv_pitch = zeros(size(c_pitch,1),3);

% Generator torque command
gen_torque = test_results.genTorqueSetpointActual - 1.867*10^7;

%% Load in Platform Model
% Define Path
platform_folder = '1 - Platform';
platform_dir = sprintf('%s\\%s',linear_dir,platform_folder);

% Load in raw files
load(sprintf('%s\\FOCAL_C4_A.mat',platform_dir),'A');
load(sprintf('%s\\FOCAL_C4_B.mat',platform_dir),'B');
load(sprintf('%s\\FOCAL_C4_C.mat',platform_dir),'C');
load(sprintf('%s\\FOCAL_C4_D.mat',platform_dir),'D');
load(sprintf('%s\\FOCAL_C4_X_OP.mat',platform_dir));
load(sprintf('%s\\FOCAL_C4_Y_OP.mat',platform_dir));

% Remove rotor azimuth state & select inputs
state_range = [1:10,12:22];
position_range = [1:6];
velocity_range = [11:16];
control_range = [301,2107:2112,2191:2193,2195,2196];

x_OP = x_OP(state_range);

A = A(state_range,state_range);
B = B(state_range,control_range);
C = C(:,state_range);
D = 0*D(:,control_range);

% Scale outputs
% C(44:46,:) = C(44:46,:)*10^-5; % convert moorings to dN
% C(12:17,:) = C(12:17,:)*10^-3; % Tower Base Forces/Moments to MN

% Discretize Platform
platform_sys_c = ss(A,B,C,D);
platform_sys_d = c2d(platform_sys_c,dt,'zoh');
[A_platform,B_platform,C_platform,D_platform] = ssdata(platform_sys_d);

% Clear out A,B,C,D matrices
clear A B C D platform_sys_d platform_sys_c

%% Load in Hydrodynamics Model (FC4)
% Define Path
hydro_folder = '2 - Hydrodynamics';
hydro_dir = sprintf('%s\\%s',linear_dir,hydro_folder);

% Load in raw files
load(sprintf('%s\\FOCAL_C4_HD_A.mat',hydro_dir),'A');
load(sprintf('%s\\FOCAL_C4_HD_B.mat',hydro_dir),'B');
load(sprintf('%s\\FOCAL_C4_HD_C.mat',hydro_dir),'C');
load(sprintf('%s\\FOCAL_C4_HD_D.mat',hydro_dir),'D');
load(sprintf('%s\\FOCAL_C4_Hydro_OP.mat',hydro_dir),'Hydro_OP');

% Convert Hydro_OP type
Hydro_OP = cell2mat(Hydro_OP);

% Trim inputs
B = B(:,[37,7,8,9,10,11,12]);
D = D(:,[37,7,8,9,10,11,12]);

% Discretize hydrodynamics model
hydro_sys_c = ss(A,B,C,D);
hydro_sys_d = c2d(hydro_sys_c,dt,'zoh');
[A_hydro,B_hydro,C_hydro,D_hydro] = ssdata(hydro_sys_d);

% Clear extra variables
clear A B C D hydro_sys_c hydro_sys_d

%% Define System Measurements for Correction
% Combine measurements to single matrix
pitch = test_results.PtfmPitch;
roll = test_results.PtfmRoll;

FA_bending = test_results.towerBotMy*10^-6; % MN
SS_bending = test_results.towerBotMx*10^-6; % MN

rotor_speed = test_results.genSpeed*(30/pi);

FA_nacelle_acceleration = test_results.accelNacelleAx;

system_measurements = [pitch,roll,rotor_speed,FA_nacelle_acceleration];


% Form measurement function (H) from SS output
H = C_platform([22,21,8,32],:); % angular displacement
h_OP = y_OP([22,21,8,32]);
% H = C_platform([9,8,4,39,40,41],:); % tower bending

clear pitch roll pitch_acceleration roll_acceleration rotor_speed mooring_tension_1 mooring_tension_2 mooring_tension_3

%% Compute Measurement Covariance Matrix
% Low-pass filter @ 6Hz
f_sample = length(test_results.Time)/max(test_results.Time);
filtered_measurements = lowpass(system_measurements,1,f_sample);
% filtered_measurements = system_measurements;
filtered_measurements = [getRamp(zeros(1,size(system_measurements,2)),filtered_measurements(1,:),714);...
                         filtered_measurements];

measurement_noise = highpass(system_measurements,1,f_sample);

% Covariance of measurements
measurement_covariance = var(measurement_noise(23906:end-1000,:));
R = diag(measurement_covariance);

% R(1,1) = 0.0031;
% R(2,2) = 0.0031;
% R(4,4) = 0.4;
% R(5,5) = 0.4;
% R(6,6) = 0.65;
% R(7,7) = 0.04;

R = 0.001*eye(size(R));

%% Load in P & Q Matrices
kalman_dir = 'C:\Umaine Google Sync\GitHub\FOWT_Optimal_Control\Models\FOCAL_C4\Linear_Files\5 - Kalman Files';

load(sprintf('%s\\FC4_Q.mat',kalman_dir));
load(sprintf('%s\\FC4_P.mat',kalman_dir));
% P = 0*P;
P = zeros(size(A_platform));
% qi = [1,2,3,4,5,6,11,12,13,14,15,16,21];

% Q = Q(qi,qi);

Q = eye(size(Q));

%% Adjust Q Values
% Q(7:10,7:10) = 10*Q(7:10,7:10);
% Q(17:20,17:20) = 10*Q(17:20,17:20);

%% Prepare Platform Operating Point
OP = [6.82;
      -0.827;
      -0.1;
      3.11; % roll (4)
      5.45; % pitch (5)
      -0.532;
      7.5578589; % gen speed (7)
      -0.01703];  % Nacelle FA Ax (8)

%% Simulate System (Kalman Filter)
disp('Beginning Kalman filter simulation...')
% Initialization (zero IC)
if exist('test_time','var')
    ss_time = test_time;
else
    ss_time = sim_time;
end
x_HD = zeros(size(A_hydro,1),1);
x = zeros(size(A_platform,1),1);
Y = zeros(size(C_platform,1),length(ss_time)-1);

% Loop over simulation time
for i = 1:length(ss_time)-1

    % Separate platform position/velocity
    platform_positions = x(position_range);
    platform_velocities = x(velocity_range);

    % Define HydroDyn Input
    u_hydro = [eta(i);
               platform_velocities];

    % Update HydroDyn States
    x_HD = A_hydro*x_HD + B_hydro*u_hydro;

    % Extract resultant forces for platform input
    hydro_out = C_hydro*x_HD + D_hydro*u_hydro;
    platform_forces = hydro_out(2:end);
    % platform_forces = hydro_out;

    % Form platform input vector
    u_platform = [wind(i);
                  platform_forces;
                  idv_pitch(i,:)';
                  gen_torque(i);
                  c_pitch(i)];

    % Do prediction step
    [x,P] = predict(x,P,A_platform,dGain(1,Q),B_platform,u_platform);

    % Get "measurements"
    z = filtered_measurements(i,:);

    % Do update step
    [x,P,K] = update(H,P,R,z',x,h_OP);

    % Store platform outputs
    Y(:,i) = C_platform*x + y_OP;
end

Y = Y';

% Save P for inspection
Pk = P;

%% Simulate System (State-Space)
disp('Beginning state-space simulation...')
% Initialization (zero IC)
if exist('test_time','var')
    ss_time = test_time;
else
    ss_time = sim_time;
end
x_HD = zeros(size(A_hydro,1),1);
x = zeros(size(A_platform,1),1);
Y_raw = zeros(size(C_platform,1),length(ss_time)-1);
platform_positions = zeros(6,1);
platform_velocities = zeros(6,1);


% Loop over simulation time
for i = 1:length(ss_time)-1

    % Separate platform position/velocity
    platform_positions = x(position_range);
    platform_velocities = x(velocity_range);

    % Define HydroDyn Input
    u_hydro = [eta(i);
               platform_velocities];

    % Update HydroDyn States
    x_HD = A_hydro*x_HD + B_hydro*u_hydro;

    % Extract resultant forces for platform input
    hydro_out = C_hydro*x_HD;
    platform_forces = hydro_out(2:end);
    % platform_forces = hydro_out;

    % Form platform input vector
    u_platform = [wind(i);
                  platform_forces;
                  idv_pitch(i,:)';
                  gen_torque(i);
                  c_pitch(i)];

    % Do prediction step
    [x,P] = predict(x,P,A_platform,dGain(1,Q),B_platform,u_platform);

    % Store platform outputs
    Y_raw(:,i) = C_platform*x + y_OP;
end

Y_raw = Y_raw';

%% Plot Results
close all; 

% Plot parameters
tmax = 7500;

% Plot Platform Surge
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Surge')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_raw(:,18),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,18),'DisplayName','Kalman')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmSurge,'DisplayName','Experiment')
end
legend

% Plot Platform Sway
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Sway')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_raw(:,19),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,19),'DisplayName','Kalman')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmSway,'DisplayName','Experiment')
end
legend

% Plot Platform Heave
figure
% subplot(4,1,1)
gca; hold on; box on;
xlim([0,tmax])
title('Platform Heave [m]')
plot(ss_time(1:end-1)-29.95,Y_raw(:,20),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,20),'DisplayName','Kalman');
% plot(sim_time,sim_results.PtfmHeave,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmHeave,'DisplayName','Experiment')
end
legend

% Plot Platform Roll
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Roll')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_raw(:,21),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,21),'DisplayName','Kalman')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmRoll,'DisplayName','Experiment')
end
legend

% Plot Platform Pitch
figure
% subplot(4,1,1)
gca; hold on; box on;
xlim([0,tmax])
title('Platform Pitch [deg]')
plot(ss_time(1:end-1)-29.95,Y_raw(:,22),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,22),'DisplayName','Kalman');
% plot(sim_time,sim_results.PtfmPitch,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmPitch,'DisplayName','Experiment')
end
legend

% Plot Platform Yaw
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Yaw')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_raw(:,23),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,23),'DisplayName','Kalman')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
try
    plot(test_time,test_results.PtfmYaw,'DisplayName','Experiment')
end
legend

% Plot tower fore-aft bending moment
figure
gca; hold on; box on;
xlim([0,tmax])
title('Tower FA Bending Moment [kN-m]')
plot(ss_time(1:end-1)-29.95,Y_raw(:,13),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,13),'DisplayName','Kalman');
% plot(sim_time,sim_results.TwrBsMyt,'DisplayName','OpenFAST');
try
    plot(test_time,test_results.towerBotMy*10^-3,'DisplayName','Experiment')
end
legend

% % Plot tower bending spectrum
% tpsd = myPSD(Y(:,13),f_sample,25);
% epsd = myPSD(test_results.towerBotMy*10^-6,f_sample,25);
% spsd = myPSD(Y_raw(:,13),f_sample,25);
% rat = tpsd(:,2)./epsd(1:end-1,2);
% diffpsd = tpsd(:,2)/mean(rat(1405:2071));
% figure; gca; hold on;
% title('Tower Bending PSD')
% plot(spsd(:,1),spsd(:,2),'DisplayName','State-Space');
% plot(tpsd(:,1),tpsd(:,2),'DisplayName','Kalman Filter'); 
% plot(epsd(:,1),epsd(:,2),'DisplayName','Experiment');
% plot(tpsd(:,1),diffpsd,'DisplayName','Kalman Scaled');
% xlim([0,0.2]); 
% ylim([0,6.6*10^5]);
% legend

% Plot rotor speed
figure
% subplot(4,1,3)
gca; hold on; box on;
xlim([0,tmax])
title('Rotor Speed [RPM]')
% xlim([0 500])
plot(ss_time(1:end-1)-29.95,Y_raw(:,8),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,8),'DisplayName','Kalman');
% plot(sim_time,sim_results.RotSpeed,'DisplayName','OpenFAST')
try
    plot(test_time,(test_results.genSpeed*(30/pi)),'DisplayName','Experiment');
end
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (1)')
plot(ss_time(1:end-1)-29.95,Y_raw(:,44),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,44),'DisplayName','Kalman Filter');
plot(test_time,test_results.leg1MooringForce-3.28*10^5,'DisplayName','Experiment');
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (2)')
plot(ss_time(1:end-1)-29.95,Y_raw(:,45),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,45),'DisplayName','Kalman Filter');
plot(test_time,test_results.leg2MooringForce-3.17*10^5,'DisplayName','Experiment');
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (3)')
plot(ss_time(1:end-1)-29.95,Y_raw(:,46),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,46),'DisplayName','Kalman Filter');
plot(test_time,test_results.leg3MooringForce-3.4*10^5,'DisplayName','Experiment');
legend

% Plot nacelle acceleration
figure
gca; hold on; box on;
title('Nacelle FA Acceleration')
plot(ss_time(1:end-1)-29.95,Y_raw(:,32),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,32),'DisplayName','Kalman Filter');
plot(test_time,test_results.accelNacelleAx,'DisplayName','Experiment');
legend

% Plot wave elevation
figure
% subplot(4,1,3)
gca; hold on; box on;
xlim([0,tmax])
title('Wave Elevation [m]')
plot(ss_time,eta)

% Plot Wind Speed
figure
% subplot(4,1,4)
gca; hold on; box on;
plot(ss_time,wind,'DisplayName','Experiment')
title('Wind Speed [m/s]')
xlabel('Time [s]')
legend

%% Clear window
clc

%% Functions --------------------------------------------------------- %%
% Prediction (Labbe, 2020, pg 212)
function [x,P] = predict(x,P,F,Q,B,u)
    x = F*x + B*u; %predict states
    P = F*P*F' + Q; %predict process covariance
end

% Update
function [x,P,K] = update(H,P,R,z,x,OP)
    S = H*P*H' + R; % Project system uncertainty into measurement space & add measurement uncertainty
    K = P*H'*inv(S);
    y = z-(H*x + OP); % Error term
    x = x+K*y;
    KH = K*H;
    P = (eye(size(KH))-KH)*P;
end

% Mean relative error
function e = MRE(model,experiment)
    
    % Remove mean from each
    model = rMean(model);
    experiment = rMean(experiment);

    % Compute mean relative error
    e = mean(abs(model-experiment)/(max(experiment)-min(experiment)));
end