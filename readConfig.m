function [settings] = readConfig(configFilePath)

    try
        configData = jsondecode(fileread(configFilePath));
    catch
        error('Error reading or parsing the JSON configuration file.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           LOAD DATA & Calibration       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scenarioName = configData.ScenarioName;
    settings.scenarioName = scenarioName;
    
    settings.dataPath = fullfile(configData.DataFolderPath,[scenarioName, '.mat']);
    
    calibrationFileName = configData.CalibrationFileName;
    settings.calibrationFilePath = fullfile(configData.CalibrationFolderPath,[calibrationFileName, '.mat']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             MAG & ZUPT & POS Aiding           %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings.magAiding   = configData.MagAiding.IsTurnedOn;
    settings.magAidingStartTime = configData.MagAiding.StartTime;
    settings.magAidingEndTime = configData.MagAiding.EndTime;
    
    settings.zuptAiding  = configData.ZuptAiding.IsTurnedOn;
    settings.zuptAidingStartTime = configData.ZuptAiding.StartTime;
    settings.zuptAidingEndTime = configData.ZuptAiding.EndTime;
    
    settings.posAiding = configData.PosAiding.IsTurnedOn;
    settings.posAidingStartTime = configData.PosAiding.StartTime;
    settings.posAidingEndTime = configData.PosAiding.EndTime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             SENSOR PARAMETERS           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % sensor locations
    settings.sensor_locs = configData.SensorConfig.Locations';
    % active sensor number
    settings.availableSensorIdx = configData.SensorConfig.ActiveSensors';
    % number of sensors 
    settings.numSensors = length(settings.availableSensorIdx);
    % sampling frequnecy
    settings.fs = configData.SensorConfig.Frequency;
    % sampling interval
    settings.dT = 1 / settings.fs;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Observability Constrained              %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings.ObservabilityConstrained = configData.ObservabilityConstrained;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     Outlier Rejection (NOT USED NOW!)   %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %settings.threshold =  chi2inv(.95, 3 * settings.numSensors);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             POLYNOMIAL ORDER            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings.polyOrder = configData.PolynomialModel.Order;
    settings.dimTheta  = settings.polyOrder^2 + 4 * settings.polyOrder + 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             STATE MASKS                 %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %settings.magSensorBiasInclude = false;
    settings.est_acc_bias = configData.BiasEstimation.est_acc_bias;
    settings.est_gyro_bias = configData.BiasEstimation.est_gyro_bias;
    settings.est_mag_bias = configData.BiasEstimation.est_mag_bias;
    [settings.numErrorStates, settings.errorStateMask] = makeErrorStateMask(settings);
    [settings.numStates, settings.stateMask] = makeStateMask(settings);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             GRAVITY VECTOR              %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    settings.g = configData.Gravity;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    INIT  uncertainties （standard deviation)     %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    settings.init_pos_std = configData.InitUncertainty.Init_Pos_Std;                                        % Position [m]
    settings.init_vel_std = configData.InitUncertainty.Init_Vel_Std;                                        % Velocity [m/s]
    settings.init_q_std = (pi/180)*configData.InitUncertainty.Init_Ori_Std*ones(3,1);                                 % Attitude (roll,pitch,yaw) [rad]
    settings.init_accBias_std = configData.InitUncertainty.Init_Acc_Bias_Std;                               % Accelerometer biases [m/s^2]
    settings.init_gyroBias_std = (pi/180)*configData.InitUncertainty.Init_Gyro_Bias_Std;                    % Gyro biases [rad/s]       
    settings.init_coeff_std = configData.InitUncertainty.Init_Theta_Std;                                    % Coefficients    
    settings.init_magBias_std = configData.InitUncertainty.Init_Mag_Bias_Std;                               % magnetometer biases [mu T]


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             FILTER PARAMETERS           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %-------------------------------------------- -------------------%
    %--------------- Process noise covariance (Q) -------------------%
    %-------------------------------------------- -------------------%
    % IMU measurement noise
    settings.sigma_acc_w = configData.FilterParameters.STD_Acc_Noise;                                            % unit: [m/s^2]
    settings.sigma_gyro_w = pi/180*configData.FilterParameters.STD_Gyro_Noise;                                          % unit: [rad/s]

    % IMU bias random walk noise 
    settings.sigma_acc_bias_rw = configData.FilterParameters.STD_Acc_Bias_Random_Walk;                           % unit: [m/s^(5/2)]
    settings.sigma_gyro_bias_rw = pi/180*configData.FilterParameters.STD_Gyro_Bias_Random_Walk;                         % unit: [rad/s^(3/2)]
    % MAG bias random walk noise
    settings.sigma_mag_bias_rw =  configData.FilterParameters.STD_Mag_Bias_Random_Walk;                          % unit: [mu T/s^(3/2)]
    
    % Theta random walk noise
    settings.sigma_coeff_w = configData.FilterParameters.STD_Theta_Random_Walk';
    if length(settings.sigma_coeff_w) ~= settings.dimTheta
        error('Dimension of theta random walk noise does not match the dimension of theta')
    end

    %-------------------------------------------- -------------------%
    % ---------- Measurement noise covariance (R) -------------------%
    %-------------------------------------------- -------------------%
    
    % Position aiding 
    settings.RPos = (configData.FilterParameters.STD_Pos)^2 * eye(3);                     % unit: [m^2]
    % Zupt aiding
    settings.RZupt = (configData.FilterParameters.STD_Zupt)^2 * eye(3);                    % unit:[(m/s)^2]

    displayInfo(settings);

end
% 
function [numErrorStates, masks] = makeErrorStateMask(settings)
    % create error state mask based on settings
    % error state mask is a logical vector, which indicates which states
    % are included in error state vector
    % error state vector is defined as:
    %   delta_x = [delta_p; delta_v; delta_epsilon; delta_acc_bias; delta_gyro_bias; delta_theta; delta_mag_bias]
    %   where delta_p, delta_v, delta_epsilon, delta_acc_bias, delta_gyro_bias are 3x1 vectors
    %   delta_mag_bias is a (numSensors - 1) * 3 x 1 vector
    %   delta_theta is a dimTheta x 1 vector

    % Define the dimensions of each state
    dims = struct('pos', 3, 'vel', 3, 'epsilon', 3, ...
                'acc_bias', 3, 'gyro_bias', 3, 'theta', settings.dimTheta,...
                'mag_bias', 3 * (settings.numSensors - 1));

    numErrorStates = 9 + settings.dimTheta;
    % Update numErrorStates based on toggles in the settings 
    if(settings.est_acc_bias == true)
        numErrorStates = numErrorStates + dims.acc_bias;
    end
    if(settings.est_gyro_bias == true)
        numErrorStates = numErrorStates + dims.gyro_bias;
    end
    if(settings.est_mag_bias == true)
        numErrorStates = numErrorStates + dims.mag_bias;
    end
    % For each state
    idx = 0;
    for state = fieldnames(dims)'
        if(strcmp(state{1},'pos')) || strcmp(state{1},'vel') || strcmp(state{1}, 'epsilon') || strcmp(state{1},'theta')
            masks.(state{1}) = false(numErrorStates, 1);
            masks.(state{1})(idx + 1 : idx + dims.(state{1})) = true;
            idx = idx + dims.(state{1});
        else
            % If the state is to be estimated
            if settings.(['est_' state{1}])
                % Set the corresponding elements in masks to true
                masks.(state{1}) = false(numErrorStates, 1);
                masks.(state{1})(idx + 1 : idx + dims.(state{1})) = true;
                idx = idx + dims.(state{1});
            end
        end
    end
end



function [numStates, masks] = makeStateMask(settings)
    % create state mask based on settings
    % state mask is a logical vector, which indicates which states
    % are included in state vector
    % state vector is defined as:
    %   x = [p; v; q_nb; acc_bias; gyro_bias; theta; mag_bias]
    %   where p, v, acc_bias, gyro_bias are 3x1 vectors
    %   q_nb is a 4x1 quaternion
    %   mag_bias is a (numSensors - 1) * 3 x 1 vector
    %   theta is a dimTheta x 1 vector

    % Define the dimensions of each state
    dims = struct('pos', 3, 'vel', 3, 'q_nb', 4, ...
                'acc_bias', 3, 'gyro_bias', 3, 'theta', settings.dimTheta,...
                'mag_bias', 3 * (settings.numSensors - 1));

    numStates = 10 + settings.dimTheta;
    % Update numStates based on toggles in the settings 
    if(settings.est_acc_bias == true)
        numStates = numStates + dims.acc_bias;
    end
    if(settings.est_gyro_bias == true)
        numStates = numStates + dims.gyro_bias;
    end
    if(settings.est_mag_bias == true)
        numStates = numStates + dims.mag_bias;
    end
    % For each state
    idx = 0;
    for state = fieldnames(dims)'
        if(strcmp(state{1},'pos')) || strcmp(state{1},'vel') || strcmp(state{1}, 'q_nb') || strcmp(state{1},'theta')
            masks.(state{1}) = false(numStates, 1);
            masks.(state{1})(idx + 1 : idx + dims.(state{1})) = true;
            idx = idx + dims.(state{1});
        else
            % If the state is to be estimated
            if settings.(['est_' state{1}])
                % Set the corresponding elements in masks to true
                masks.(state{1}) = false(numStates, 1);
                masks.(state{1})(idx + 1 : idx + dims.(state{1})) = true;
                idx = idx + dims.(state{1});
            end
        end
    end
end


function [] = displayInfo(settings)

    fprintf('******************************************\n'),
    fprintf('         DATA & CALIB Configuration       \n'),
    fprintf('******************************************\n'),
    fprintf('Processing: %s.\n', settings.dataPath);
    fprintf('Use calibration file: %s.\n', settings.calibrationFilePath);

    fprintf('******************************************************\n'),
    fprintf('         MAG & ZUPT & POS AidingConfiguration         \n'),
    fprintf('******************************************************\n'),
    if (settings.magAiding)
        fprintf('MAG aiding starts at %d s, ends at %d s.\n', settings.magAidingStartTime, settings.magAidingEndTime);
    end
    if (settings.zuptAiding)
        fprintf('ZUPT aiding starts at %d s, ends at %d s.\n', settings.zuptAidingStartTime, settings.zuptAidingEndTime);
    end
    if (settings.posAiding)
        fprintf('POS aiding starts at %d s, ends at %d s.\n', settings.posAidingStartTime, settings.posAidingEndTime);
    end

    fprintf('******************************************\n'),
    fprintf('           SENSOR Configuration           \n'),
    fprintf('******************************************\n'),
    displaySensorConfig(settings);

    fprintf('******************************************\n'),
    fprintf('          POLYNOMIAL Configuration        \n'),
    fprintf('******************************************\n'),
    fprintf('Use %d order model, theta dimension: %d.\n', settings.polyOrder, settings.dimTheta);

    fprintf('******************************************\n'),
    fprintf('             State Configuration          \n'),
    fprintf('******************************************\n'),
    fprintf('Estimate magnetometer bias: ');
    if settings.est_mag_bias
        fprintf('Yes!\n');
    else
        fprintf('No!\n');
    end 
end


function [] = displaySensorConfig(settings)
    r=settings.sensor_locs;
    plot(r(1,:),r(2,:),'bo','MarkerSize',4);
    r=r(:,settings.availableSensorIdx);
    hold on;
    plot(r(1,:),r(2,:),'ro','MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','r')
    grid on;
    axis equal;
    set(gca, 'YDir' , 'reverse' )
    ylabel('y-axis [m]','FontName','Times New Roman');
    xlabel('x-axis [m]', 'FontName','Times New Roman');
    legend("", "Active sensor");
    title('Current configuration','FontName','Times New Roman');
end