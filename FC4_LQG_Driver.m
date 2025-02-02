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
wind_case = 2;

% Load wind file
wind_file_path = 'C:\Umaine Google Sync\GitHub\FOWT_Optimal_Control\Wind_Files';
wind = getFOCALWindVector(wind_file_path,wind_case,test_time);

% Time-shift wind to account for hydro causalization time
wind = [zeros(causality_shift_index,1);
        wind(1:end-causality_shift_index)];

%% Prepare Control Input Values
% Blade pitch command (collective)
c_pitch = 1.1*test_results.pitch1Position*(pi/180);

% Ind. Pitch
idv_pitch = zeros(size(c_pitch,1),3);

% Generator torque command
gen_torque = test_results.genTorqueSetpointActual;

%% Load in Platform Model
% Define Path
platform_folder = '1 - Platform';
platform_dir = sprintf('%s\\%s',linear_dir,platform_folder);

% Load in raw files
load(sprintf('%s\\FOCAL_C4_A.mat',platform_dir),'A');
load(sprintf('%s\\FOCAL_C4_B.mat',platform_dir),'B');
load(sprintf('%s\\FOCAL_C4_C.mat',platform_dir),'C');
load(sprintf('%s\\FOCAL_C4_D.mat',platform_dir),'D');
% load(sprintf('%s\\FOCAL_C4_Output_OP.mat',platform_dir),'y_op');

% Remove rotor azimuth state & select inputs
A = A([1:10,12:end],[1:10,12:end]);
B = B([1:10,12:end],[301,2107:2112,2191,2192,2193,2195,2196,3943:3954]);
C = C(:,[1:10,12:end]);
D = 0*D(:,[301,2107:2112,2191,2192,2193,2195,2196,3943:3954]);

% Scale outputs
C(34:36,:) = C(34:36,:)*10^-5; % convert moorings to dN
C(8:13,:) = 0.6241*C(8:13,:)*10^-3; % 1/1.2 gain to imrove ss freq response

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

%% Prepare Controller Gains
% Separate control matrices
Bcontrol = B_platform(:,[8:10]);
Dcontrol = D_platform(:,[8:10]);

% Bintegral = Bplatform(:,10);

% Form discrete time state space for integral control
% Isys = ss(A_platform,B);

% Weighting matrices
Qcontrol = zeros(size(A_platform));
    Qcontrol(4,4) = 1; % reduce roll
    Qcontrol(5,5) = 1; % reduce pitch

Rcontrol = eye(3); % 3 blade pitch values

% Compute control gains (discrete time)
[Kcontrol,Scontrol,econtrol] = dlqr(A_platform,Bcontrol,Qcontrol,Rcontrol);

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
Y_opt = zeros(size(C_platform,1),length(ss_time)-1);
platform_positions = zeros(6,1);
platform_velocities = zeros(6,1);
speed_error = 0;

% Loop over simulation time
for i = 1:length(ss_time)-1

    % Separate platform position/velocity
    platform_positions = x(1:6);
    platform_velocities = x(11:16);

    % Define HydroDyn Input
    u_hydro = [eta(i);
               platform_velocities];

    % Update HydroDyn States
    x_HD = A_hydro*x_HD + B_hydro*u_hydro;

    % Extract resultant forces for platform input
    hydro_out = C_hydro*x_HD;
    platform_forces = hydro_out(2:end);
    % platform_forces = hydro_out;

    % Compute optimal individual pitch commands
    opt_pitch = -Kcontrol*x;

    % Compute collective pitch command
    % col_pitch = col_pitch + Igain*speed_error;

    % Form platform input vector
    u_platform = [wind(i);
                  platform_forces;
                  opt_pitch;
                  gen_torque(i);
                  c_pitch(i);
                  platform_positions;
                  platform_velocities];

    % Do prediction step
    [x] = predict(x,A_platform,B_platform,u_platform);

    % Store platform outputs
    Y_opt(:,i) = C_platform*x;

    % Add to integral rotor speed error
    rotor_speed_target = (15*10^6)/gen_torque(i)
    speed_error = speed_error + (rotor_speed_target - Y_opt(4,i))
end

Y_opt = Y_opt';

%% Simulate System - Baseline State Space
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
    platform_positions = x(1:6);
    platform_velocities = x(11:16);

    % Define HydroDyn Input
    u_hydro = [eta(i);
               platform_velocities];

    % Update HydroDyn States
    x_HD = A_hydro*x_HD + B_hydro*u_hydro;

    % Extract resultant forces for platform input
    hydro_out = C_hydro*x_HD;
    platform_forces = hydro_out(2:end);
    % platform_forces = hydro_out;

    % Compute optimal individual pitch commands
    opt_pitch = -Kcontrol*x;

    % Form platform input vector
    u_platform = [wind(i);
                  platform_forces;
                  opt_pitch*0;
                  gen_torque(i);
                  c_pitch(i);
                  platform_positions;
                  platform_velocities;];

    % Do prediction step
    [x] = predict(x,A_platform,B_platform,u_platform);

    % Store platform outputs
    Y_raw(:,i) = C_platform*x;
end
Y = Y_raw';

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
plot(ss_time(1:end-1)-29.95,Y_opt(:,14),'DisplayName','Optimal')
plot(ss_time(1:end-1)-29.95,Y(:,14),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmSurge,'DisplayName','Experiment')
% end
legend

% Plot Platform Sway
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Sway')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_opt(:,15),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,15),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmSway,'DisplayName','Experiment')
% end
legend

% Plot Platform Heave
figure
% subplot(4,1,1)
gca; hold on; box on;
xlim([0,tmax])
title('Platform Heave [m]')
plot(ss_time(1:end-1)-29.95,Y_opt(:,16),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,16),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmHeave,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmHeave,'DisplayName','Experiment')
% end
legend

% Plot Platform Roll
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Roll')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_opt(:,17),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,17),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmRoll,'DisplayName','Experiment')
% end
legend

% Plot Platform Pitch
figure
% subplot(4,1,1)
gca; hold on; box on;
xlim([0,tmax])
title('Platform Pitch [deg]')
plot(ss_time(1:end-1)-29.95,Y_opt(:,18),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,18),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmPitch,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmPitch,'DisplayName','Experiment')
% end
legend

% Plot Platform Yaw
figure
% subplot(4,1,2)
gca; hold on; box on;
title('Platform Yaw')
xlim([0,tmax])
plot(ss_time(1:end-1)-29.95,Y_opt(:,19),'DisplayName','State-Space')
plot(ss_time(1:end-1)-29.95,Y(:,19),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.PtfmSurge,'DisplayName','OpenFAST')
% try
%     plot(test_time,test_results.PtfmYaw,'DisplayName','Experiment')
% end
legend

% Plot tower fore-aft bending moment
figure
gca; hold on; box on;
xlim([0,tmax])
title('Tower FA Bending Moment [kN-m]')
plot(ss_time(1:end-1)-29.95,Y_opt(:,9),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,9),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.TwrBsMyt,'DisplayName','OpenFAST');
% try
%     plot(test_time,test_results.towerBotMy*10^-6,'DisplayName','Experiment')
% end
legend

% % Plot tower bending spectrum
% tpsd = myPSD(Y(:,9),f_sample,25);
% epsd = myPSD(test_results.towerBotMy*10^-6,f_sample,25);
% spsd = myPSD(Y_raw(:,9),f_sample,25);
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

% Plot generator power
figure; 
gca; hold on; box on;
xlim([0,tmax]);
title('Generator Power [kW]');
plot(ss_time(1:end-1)-29.95,-Y_opt(:,2),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,-Y(:,2),'DisplayName','Baseline SS')
% plot(test_results.Time,test_results.genSpeed.*test_results.genTorque*10^-6)
legend

% Plot rotor speed
figure
% subplot(4,1,3)
gca; hold on; box on;
xlim([0,tmax])
title('Rotor Speed [RPM]')
% xlim([0 500])
plot(ss_time(1:end-1)-29.95,Y_opt(:,4),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,4),'DisplayName','Baseline SS')
% plot(sim_time,sim_results.RotSpeed,'DisplayName','OpenFAST')
% try
%     plot(test_time,(test_results.genSpeed*(30/pi)),'DisplayName','Experiment');
% end
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (1)')
plot(ss_time(1:end-1)-29.95,Y_opt(:,34),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,34),'DisplayName','Baseline SS')
% plot(test_time,10^-5*test_results.leg1MooringForce-3.28*10^1,'DisplayName','Experiment');
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (2)')
plot(ss_time(1:end-1)-29.95,Y_opt(:,35),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,35),'DisplayName','Baseline SS')
% plot(test_time,10^-5*test_results.leg2MooringForce-3.17*10^1,'DisplayName','Experiment');
legend

% Plot lead mooring tension
figure
gca; hold on; box on;
title('Mooring Tension (3)')
plot(ss_time(1:end-1)-29.95,Y_opt(:,36),'DisplayName','State-Space');
plot(ss_time(1:end-1)-29.95,Y(:,36),'DisplayName','Baseline SS')
% plot(test_time,10^-5*test_results.leg3MooringForce-3.4*10^1,'DisplayName','Experiment');
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
function [x] = predict(x,F,B,u)
    x = F*x + B*u; %predict states
end

% Mean relative error
function e = MRE(model,experiment)
    
    % Remove mean from each
    model = rMean(model);
    experiment = rMean(experiment);

    % Compute mean relative error
    e = mean(abs(model-experiment)/(max(experiment)-min(experiment)));
end