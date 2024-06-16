function out_data=magaidedINS(in_data,settings)

    % Copy data to variables with shorter name
    u=in_data.u;
    mag=in_data.mag_array.field;
    t=in_data.t;
    gt=in_data.gt;

    magEnable = settings.magAiding;
    magStartTime = settings.magAidingStartTime;
    magEndTime = settings.magAidingEndTime;

    zuptEnable  = settings.zuptAiding;
    zuptStartTime = settings.zuptAidingStartTime;
    zuptEndTime = settings.zuptAidingEndTime;

    posEnable = settings.posAiding;
    posStartTime = settings.posAidingStartTime;
    posEndTime = settings.posAidingEndTime;

    stateMask = settings.stateMask;

    % Get measurement equations
    [H, HPos, HZupt, Phi]=getHs(settings);

    %% Initialization
    % Initialize the navigation state
    z = reshape(mag(1, :), 3, []);
    z = z(:, settings.availableSensorIdx);
    % initiate polynomial magnetic model
    m = polyMagModel(settings.polyOrder);
    r = settings.sensor_locs;
    m = m.set_phi(r(:, settings.availableSensorIdx));
    m = m.init_theta(z(:));
    [x_h, P]=init_navigation_state(m,  u(1, 1:3), gt, settings);

    % Get process noise covariance and measurement noise covariance
    Q = getProcessNoiseCov(settings);
    settings.Q = Q;
    RPos = settings.RPos;
    RZupt = settings.RZupt;

    % Allocate memory for the output data
    N=size(u,1);
    out_data.x_h=zeros(settings.numStates,N);
    out_data.x_h(:,1)=x_h;
    out_data.diag_P=zeros(settings.numErrorStates,N);
    out_data.diag_P(:,1)=diag(P);
    out_data.cov_p=zeros(3,3,N);
    out_data.cov_p(:,:,1)=P(settings.errorStateMask.pos, settings.errorStateMask.pos);
    out_data.cov_q=zeros(3,3,N);
    out_data.cov_q(:,:,1)=P(settings.errorStateMask.epsilon, settings.errorStateMask.epsilon);

    out_data.res=zeros(3*settings.numSensors, N);
    out_data.nis=zeros(1, N);
    out_data.debug.theta = zeros(settings.dimTheta, N);
    out_data.debug.theta_pred = zeros(settings.dimTheta, N);
    out_data.debug.theta(:, 1) = x_h(end-settings.dimTheta+1 : end);
    out_data.debug.theta_pred(:, 1) = x_h(end-settings.dimTheta+1 : end);
    out_data.debug.noisevar = zeros(N, 1);

    % ob = zeros(15*settings.numErrorStates, settings.numErrorStates);
    % ob_idx = 16;
    % ob(1:15, :) = H;
    xk_km1 = x_h;
    null_k= [zeros(3,1); ...
         -vect2skew(xk_km1(stateMask.vel))*settings.g;...
         q2r(xk_km1(stateMask.q_nb))'*settings.g;...
         zeros(6, 1); ...
         zeros(settings.dimTheta,1)];
    %% Information fusion
    for k=2:N
        
        % Sampling period
        Ts=t(k)-t(k-1);

        % calculate psi
        psi = calc_psi(x_h, u(k - 1, :), Ts, settings);

        % update theta
        m = m.update_theta(psi);
        out_data.debug.theta_pred(:, k) = m.get_theta();
        % Get state space model matrices
        if(settings.ObservabilityConstrained)  && (t(k) >= posStartTime)      
            %[F,G]=state_space_model_null2(x_h, xk_km1, u(k - 1, :), Ts, m, settings);
            %[F, G] = state_space_model_null_opt_q(x_h, xk_km1, u(k - 1, :), Ts, m, settings);
            [F, G] = state_space_model_null_final(x_h, xk_km1, u(k - 1, :), Ts, m, settings);
        else
            [F,G]=state_space_model(x_h, u(k - 1, :), Ts, m, settings);
        end
        %[F,G]=state_space_model(x_h, u(k - 1, :), Ts, m, settings);
        %[F,G]=state_space_model_null(x_h, xk_km1, u(k - 1, :), Ts, m, settings);
        % Check observality
        % null_k=F*null_k;
        % % check H*null_k is zero matrix
        % if all(abs(H*null_k) < 1e-6, 'all')
        %     disp('H*null_k is a zero matrix')
        % else
        %     disp('H*null_k is not a zero matrix')
        % end
        % disp(max(abs(H*null_k)))
        % ob(ob_idx:ob_idx+14, :) = ob(ob_idx - 15:ob_idx-1, :) * F;
        % ob_idx = ob_idx + 15;
        % rank_ob = rank(ob);
        % if rank_ob > 19 && (k>settings.numErrorStates)
        %     fprintf('Gain rank! Current timestep: %d, rank: %d!\n', k, rank_ob);
        % end
        % % Shift elements in ob if it's full
        % if ob_idx > size(ob, 1)
        %     ob(1:end-15, :) = ob(16:end, :);
        %     ob_idx = ob_idx - 15;
        % end
        % Update the nominal state state
        [x_h]=nominalStateProp(x_h, u(k - 1, :), Ts, m, settings);
        % save xk_km1 for later use
        xk_km1 = x_h; 
        %disp(null_k(4:6) + vect2skew(xk_km1(stateMask.vel))*settings.g)
        %null_k(4:6) = -vect2skew(xk_km1(stateMask.vel))*settings.g;
        %disp(null_k(7:9) - q2r(xk_km1(stateMask.q_nb))'*settings.g)
        %null_k(7:9) = q2r(xk_km1(stateMask.q_nb))'*settings.g;
        % Time update of the Kalman filter state covariance.
        P=F * P * F'+ G * Q * G';
        
        % Position aiding update
        if posEnable && (t(k) >= posStartTime && t(k) <= posEndTime) && (~any(isnan(gt.pos(:, k)), 'all')) 
            K = P * HPos' / (HPos * P * HPos' + RPos); 
            delta_z = gt.pos(:, k) - x_h(stateMask.pos);
            P = (eye(size(P)) - K * HPos) * P * (eye(size(P)) - K * HPos)' + K * RPos * K';
            delta_x = K * delta_z;
            [x_h, P]=correctNominalState(x_h, P, delta_x, settings);
        end    

        % Zupt aiding update
        if zuptEnable && (t(k) >= zuptStartTime && t(k) <= zuptEndTime) && (zupt(k)) 
            K = P * HZupt' / (HZupt * P * HZupt' + RZupt); 
            delta_z = [0; 0; 0] - x_h(stateMask.vel);
            P = (eye(size(P)) - K * HZupt) * P * (eye(size(P)) - K * HZupt)' + K * RZupt * K';
            delta_x = K * delta_z;
            [x_h, P]=correctNominalState(x_h, P, delta_x, settings);
        end    

        % Mag aiding update
        if magEnable && (t(k) >= magStartTime && t(k) <= magEndTime)
            z = reshape(mag(k, :), 3, []);
            z = z(:, settings.availableSensorIdx);
            delta_z = z(:) - Phi * x_h(settings.stateMask.theta);
            if(settings.est_mag_bias)
                delta_z = delta_z - [x_h(settings.stateMask.mag_bias); zeros(3, 1)];
            end
            out_data.res(:, k) = delta_z;
            [trueCoeff, resVar] = m.LS_coeff(z(:));
            out_data.debug.noisevar(k) = resVar;
            out_data.debug.theta(:, k) = trueCoeff;
            R =  resVar * eye(3*settings.numSensors);
            K = P * H' / (H * P * H' + R);   
            P = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K * R * K';
            delta_x = K * delta_z;
            [x_h, P]=correctNominalState(x_h, P, delta_x, settings);
            % update theta in model 
            m = m.set_theta(x_h(end-settings.dimTheta + 1 : end));
        end

        % Save the data to the output data structure
        out_data.x_h(:,k)=x_h;
        out_data.diag_P(:,k)=diag(P);
        out_data.cov_q(:,:,k) = P(settings.errorStateMask.epsilon, settings.errorStateMask.epsilon);
        out_data.cov_p(:,:,k) = P(settings.errorStateMask.pos, settings.errorStateMask.pos);
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          SUB-FUNCTIONS                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Process Noise  Covariance     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q] = getProcessNoiseCov(settings)
    Q = blkdiag(settings.sigma_acc_w^2 * eye(3), ...
        settings.sigma_gyro_w^2 * eye(3));
    
    if settings.est_acc_bias
        Q = blkdiag(Q, settings.sigma_acc_bias_rw^2 * eye(3));
    end
    
    if settings.est_gyro_bias
        Q = blkdiag(Q, settings.sigma_gyro_bias_rw^2 * eye(3));
    end
    
    Q = blkdiag(Q, diag(settings.sigma_coeff_w.^2));
    
    if(settings.est_mag_bias)
        Q = blkdiag(Q, settings.sigma_mag_bias_rw^2 * eye((settings.numSensors - 1) * 3));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Init navigation state     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_h, P]=init_navigation_state(m, u, gt, settings)
    % extract initial pose from gt
    if isfield(gt,'init_pose')
        init_pos = gt.init_pose(1:3) ;   
        %init_pos = gt.pos(:, 1); 
    elseif isfield(gt,'pos')
        init_pos = gt.pos(:, 1) ;   
    else
        init_pos = [0;0;0];
    end

    if isfield(gt,'init_pose')                                      
        init_vel = gt.init_pose(4:6);                                           
    else
        init_vel = zeros(3, 1);
    end



    if isfield(gt,'init_pose')                                      
        init_q = gt.init_pose(7:10);                                           
    elseif isfield(gt,'ori')
        init_q = rotm2quat(gt.ori(:, :, 1))';  
    else
        init_q = initialize_pose(u);
    end

    % Initial nominal state vector (position, velocity, orientation, inertial sensor biases, theta, magnetometer bias) 
    % assume zero initial velocity
    x_h=[init_pos; init_vel; init_q];
    if(settings.est_acc_bias)
        x_h = [x_h; zeros(3, 1)];
        %x_h = [x_h; -0.0451; 0.0391; 0.0665];
    end
    if(settings.est_gyro_bias)
        x_h = [x_h; zeros(3, 1)];
        %x_h = [x_h; -1.71e-4; -5.29e-4; 1.4e-3];
    end
    if settings.magAiding
        x_h=[x_h; m.get_theta()];
    else
        x_h=[x_h; zeros(settings.dimTheta, 1)];
    end
    
    if(settings.est_mag_bias)
        x_h = [x_h; zeros(3 * (settings.numSensors - 1), 1)];
    end

    % Initial error state uncertainties (position, velocity, orientation, inertial sensor biases, theta, magnetometer bias) 
    P = diag([settings.init_pos_std^2 * ones(3, 1); 
            settings.init_vel_std^2 * ones(3, 1);  
            settings.init_q_std.^2]);
    if(settings.est_acc_bias)
        P = blkdiag(P, settings.init_accBias_std^2 * eye(3));
    end
    if(settings.est_gyro_bias)
        P = blkdiag(P, settings.init_gyroBias_std^2 * eye(3));
    end
    %        settings.init_accBias_std^2 * ones(3, 1); 
    %        settings.init_gyroBias_std^2 * ones(3, 1);
    P = blkdiag(P, settings.init_coeff_std^2 * eye(settings.dimTheta));
    if(settings.est_mag_bias)
        P = blkdiag(P, settings.init_magBias_std^2 * eye(3 * (settings.numSensors - 1)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Calculate Psi       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = calc_psi(x_h, u, dt, settings)
    mask = settings.stateMask;
    
    vk = x_h(mask.vel);
    R_nb = q2r(x_h(mask.q_nb));
    if(settings.est_acc_bias)
        acc_bias = x_h(mask.acc_bias);
    else
        acc_bias = zeros(3, 1);
    end
    %acc_bias = x_h(mask.acc_bias);
    if(settings.est_gyro_bias)
        gyro_bias = x_h(mask.gyro_bias);
    else
        gyro_bias = zeros(3, 1);
    end
    %gyro_bias = x_h(mask.gyro_bias);
    
    %  parse u
    acc_m = u(1:3)';
    omega_m = u(4:end)';

    acc_nav = R_nb * (acc_m - acc_bias) + settings.g;
    dp = vk * dt + 1/2 * acc_nav * dt^2;
    dp_body = R_nb.' * dp;

    omega_h = omega_m - gyro_bias;
    rotvec = omega_h * dt;

    psi = [dp_body; rotvec];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Nominal state propogation        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = nominalStateProp(xk, u, dt, m, settings)
    %   INPUT:
    %                 xk:  current state estimate at time k
    %                  u:  accelerometer reading and gyroscope reading
    %                 dt:  sampling interval
    %           settings:  struct setting 
    %   
    %   OUTPUT:
    %                  x:  predict nominal state at time k + 1            

        masks = settings.stateMask;        
        
        % parse xk
        pk = xk(masks.pos);
        vk = xk(masks.vel);
        q_nb = xk(masks.q_nb);
        R_nb = q2r(q_nb);
        if(settings.est_acc_bias)
            acc_bias = xk(masks.acc_bias);
        else
            acc_bias = zeros(3, 1);
        end
        %acc_bias = xk(masks.acc_bias);
        if(settings.est_gyro_bias)
            gyro_bias = xk(masks.gyro_bias);
        else
            gyro_bias = zeros(3, 1);
        end
        %gyro_bias = xk(masks.gyro_bias);
        %  parse u
        acc_m = u(1:3)';
        omega_m = u(4:end)';
    
        acc_nav = R_nb * (acc_m - acc_bias) + settings.g;
        dp = vk * dt + 1/2 * acc_nav * dt^2;
        
        % nominal state \hat{x}_k as in eq. 2b and eq. 14
        x = zeros(size(xk));
        x(masks.pos) = pk + dp;
        x(masks.vel) = vk + acc_nav * dt;
        if(settings.est_acc_bias)
            x(masks.acc_bias) = acc_bias;
        end
        %x(masks.acc_bias) = acc_bias;
        if(settings.est_gyro_bias)
            x(masks.gyro_bias) = gyro_bias;
        end
        %x(masks.gyro_bias) = gyro_bias;
        
        omega_h = omega_m - gyro_bias;
        rotvec = omega_h * dt;
        dq = rotvec2quat(rotvec');
        qw = dq(1);
        qv = dq(2:4)';
        x(masks.q_nb) = (qw*eye(4) + [0 -qv'; qv -vect2skew(qv)]) * q_nb;
        x(masks.theta) = m.get_theta();

        if(settings.est_mag_bias)
            x(masks.mag_bias) = xk(masks.mag_bias);
        end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Measurement Equation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, HPos, HZupt, Phi]=getHs(settings)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             Measurement Matrix          %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Phi = [];
    for i = 1 : length(settings.availableSensorIdx)
        idx = settings.availableSensorIdx(i);
        Phi =  [Phi; get_ABnull(settings.sensor_locs(:, idx), settings.polyOrder)];
    end
    
    H = zeros(settings.numSensors * 3, settings.numErrorStates);
    H(:, settings.errorStateMask.theta) = Phi;
    if(settings.est_mag_bias)
        H(1:end-3, settings.errorStateMask.mag_bias) = eye(3 * (settings.numSensors - 1));
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             Auxiliary aiding            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HZupt = zeros(3, settings.numErrorStates);
    HZupt(:, settings.errorStateMask.vel) = eye(3);

    HPos = zeros(3, settings.numErrorStates);
    HPos(:, settings.errorStateMask.pos) = eye(3);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Nominal State Correction         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xh, P]=correctNominalState(xh, P, delta_x, settings)
    
    stateMask = settings.stateMask;
    errorStateMask =  settings.errorStateMask;

    xh(stateMask.pos) = xh(stateMask.pos) +  delta_x(errorStateMask.pos);
    xh(stateMask.vel) = xh(stateMask.vel) +  delta_x(errorStateMask.vel);
    qv = 1/2 * delta_x(errorStateMask.epsilon);
    xh(stateMask.q_nb) = (eye(4) + [0 -qv'; qv -vect2skew(qv)]) * xh(stateMask.q_nb);
    xh(stateMask.q_nb) = xh(stateMask.q_nb).' / norm(xh(stateMask.q_nb)); 
    if(settings.est_acc_bias)
        xh(stateMask.acc_bias) = xh(stateMask.acc_bias) +  delta_x(errorStateMask.acc_bias);
    end
    %xh(stateMask.acc_bias) = xh(stateMask.acc_bias) +  delta_x(errorStateMask.acc_bias);
    if(settings.est_gyro_bias)
        xh(stateMask.gyro_bias) = xh(stateMask.gyro_bias) +  delta_x(errorStateMask.gyro_bias);
    end
    %xh(stateMask.gyro_bias) = xh(stateMask.gyro_bias) +  delta_x(errorStateMask.gyro_bias);
    xh(stateMask.theta)     = xh(stateMask.theta) +  delta_x(errorStateMask.theta);

    if(settings.est_mag_bias)
        xh(stateMask.mag_bias)     = xh(stateMask.mag_bias) +  delta_x(errorStateMask.mag_bias);
    end

    % G = blkdiag(eye(6), eye(3) - 1 / 2 * vect2skew(delta_x(7:9)) ,eye(6 + settings.dimTheta));
    
    % if(settings.magSensorBiasInclude)
    %     G = blkdiag(G, eye(3 * (settings.numSensors - 1)));
    % end

    %P = G * P * G.';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Outlier rejection (NOT IN USE!)    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isOutlier, nis] = outlier(mag, x_h, Phi, H, P, R, settings)
    S = (H * P * H' + R);   
    z = reshape(mag, 3, []);
    z = z(:, settings.availableSensorIdx);
    delta_z = z(:) - Phi * x_h(settings.stateMask.theta);
    if(settings.magSensorBiasInclude)
        delta_z = delta_z - [x_h(settings.stateMask.mag_bias); zeros(3, 1)];
    end
    
    nis = delta_z' / S * delta_z;
    if nis > settings.threshold 
        isOutlier = true;
    else
        isOutlier = false;
    end
end

function [q0] = initialize_pose(u)
    % initialize orientation based accelerometer
    yaw   = 0;      
    pitch = atan2(-u(1), sqrt(u(2)^2 + u(3)^2));
    roll  = atan2(u(2), u(3));
    q0 = eul2quat([yaw pitch roll])';
end