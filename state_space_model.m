function [F,G]=state_space_model(x_h,u,dt,m,settings)
    masks = settings.stateMask;
    numErrorStates = settings.numErrorStates;
    errorMasks = settings.errorStateMask;           

    % parse x_h
    vk = x_h(masks.vel);
    R_nb = q2r(x_h(masks.q_nb));
    if(settings.est_acc_bias)
        acc_bias = x_h(masks.acc_bias);
    else
        acc_bias = zeros(3,1);
    end
    if(settings.est_gyro_bias)
        gyro_bias = x_h(masks.gyro_bias);
    else
        gyro_bias = zeros(3,1);
    end
    %  parse u
    acc_m = u(1:3)';
    omega_m = u(4:end)';

    % remove bias from gyro readings
    omega_h = omega_m - gyro_bias;
    rotvec = omega_h * dt;

    J1J2 = m.J1J2;

    % M as in eq. 23d
    M = zeros(settings.dimTheta + 6, numErrorStates);                        
    M(1:settings.dimTheta, errorMasks.theta) = eye(settings.dimTheta);
    M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.vel)   = R_nb.' * dt;
    M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.epsilon)  = vect2skew(R_nb.' *  dt * (vk  + settings.g * dt / 2));
    if(settings.est_acc_bias)
        M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.acc_bias)   = -dt^2 / 2 * eye(3);
    end
    if(settings.est_gyro_bias)
        M(settings.dimTheta + 4: settings.dimTheta + 6, errorMasks.gyro_bias)  = -eye(3)*dt;
    end

    % construct F
    F = zeros(numErrorStates);
    F(errorMasks.pos, errorMasks.pos) = eye(3);
    F(errorMasks.pos, errorMasks.vel) = eye(3) * dt;
    F(errorMasks.vel, errorMasks.vel) = eye(3);
    F(errorMasks.vel, errorMasks.epsilon) = -R_nb * vect2skew(acc_m - acc_bias) * dt;
    if(settings.est_acc_bias)
        F(errorMasks.vel, errorMasks.acc_bias) = -R_nb * dt;
        F(errorMasks.acc_bias, errorMasks.acc_bias) = eye(3);
    end
    
    F(errorMasks.epsilon, errorMasks.epsilon) = eye(3) - vect2skew(rotvec);
    if(settings.est_gyro_bias)
        F(errorMasks.epsilon, errorMasks.gyro_bias) = -eye(3)*dt;
        F(errorMasks.gyro_bias, errorMasks.gyro_bias) = eye(3);
    end
    F(errorMasks.theta, :) = m.pinvA * ([m.B J1J2] * M) ;

    % construct G
    G = zeros(numErrorStates, size(settings.Q, 1));
    G(errorMasks.vel, 1:3) = eye(3) * dt;
    G(errorMasks.epsilon, 4:6) = eye(3) * dt;
    if(settings.est_acc_bias)
        G(errorMasks.acc_bias, 7:9) = eye(3) * sqrt(dt);
    end
    if(settings.est_gyro_bias)
        G(errorMasks.gyro_bias, 10:12) = eye(3) * sqrt(dt);
    end
    G(errorMasks.theta, 1 : 3) = -m.pinvA * (J1J2(:, 1:3) * (dt^2/2));
    G(errorMasks.theta, 4 : 6) = -m.pinvA * (J1J2(:, end-2:end))*dt;
    idx = 6;
    if(settings.est_acc_bias)
        idx = idx + 3;
    end
    if(settings.est_gyro_bias)
        idx = idx + 3;
    end
    G(errorMasks.theta, idx + 1 : idx + settings.dimTheta) = eye(settings.dimTheta);

    if(settings.est_mag_bias)
        G(errorMasks.mag_bias, end - 3 * (settings.numSensors - 1) + 1: end) = eye(3 * (settings.numSensors - 1)) * sqrt(dt);
    end

    % if(settings.magSensorBiasInclude)
    %     F(errorMasks.mag_bias, errorMasks.mag_bias) = eye(3 * (settings.numSensors - 1));
    %     G(errorMasks.mag_bias, end - 3 * (settings.numSensors - 1) + 1: end) = eye(3 * (settings.numSensors - 1)) * sqrt(dt);
    % end
end