function [mse, pos_cov, yaw_cov, x_cov, y_cov, z_cov] = compute_MSE_NEES(in_data, out_data)
% Compute the RMSE and NEES for the given data
%   in_data: struct with the following fields
    gt_pos = in_data.gt.pos;
    gt_ori = in_data.gt.ori;



    pos_est = out_data.x_h(1:3, :);
    q_est = out_data.x_h(7:10, :);


    delta_p = pos_est - gt_pos;
    %delta_q = zeros(3, size(pos_est, 2));
    % for i = 1 : size(pos_est, 2)
    %     delta_quat = quatmultiply(quatconj(q_est(:,i)'), rotm2quat(gt_ori(:,:,i)));
    %     if(delta_quat(1) < 0)
    %         delta_quat = -delta_quat;
    %     end
    %     delta_q(:, i) = quatlog(delta_quat);
    % end

    oriEst = out_data.x_h(7:10, :);
    oriRef = in_data.gt.ori;
    delta_q = angdiff(quat2eul(oriEst','ZYX'), rotm2eul(oriRef,'ZYX'));
    delta_q = delta_q(:, 1);                               % only consider yaw angle


    mse.pos = diag(delta_p' * delta_p);
    mse.x = delta_p(1,:)'.^2;
    mse.y = delta_p(2,:)'.^2;
    mse.z = delta_p(3,:)'.^2;
    mse.ori = delta_q .* delta_q;

    pos_cov = zeros(size(pos_est, 2), 1);
    x_cov = zeros(size(pos_est, 2), 1);
    y_cov = zeros(size(pos_est, 2), 1);
    z_cov = zeros(size(pos_est, 2), 1);
    for i = 1 : size(pos_est, 2)
        pos_cov(i, :) = trace(squeeze(out_data.cov_p(1:3, 1:3, i)));
        x_cov(i) = out_data.cov_p(1, 1, i);
        y_cov(i) = out_data.cov_p(2, 2, i);
        z_cov(i) = out_data.cov_p(3, 3, i);
    end

    yaw_cov = zeros(size(q_est, 2), 1);
    eulCovOCMAINS = zeros(in_data.numFrames, 3);
    for i = 1 : in_data.numFrames
        J = ypr_jacobian_quat(out_data.x_h(7:10, i)');
        eulCovOCMAINS(i, :) = diag(J * out_data.cov_q(:, :, i) * J');
        yaw_cov(i) = eulCovOCMAINS(i, 1);                % unit in rad
    end

end

