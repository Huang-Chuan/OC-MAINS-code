function J = ypr_jacobian_quat(quatRef)
    
    eulAngRef = quat2eul(quatRef, 'ZYX');
    J = zeros(3);

    dqx = [1, 1/2 * 1e-7, 0, 0];
    dqy = [1, 0, 1/2 * 1e-7, 0];
    dqz = [1, 0,   0,  1/2 * 1e-7];
    
    eulAng1 = quat2eul(quatmultiply(quatRef, dqx), 'ZYX');
    eulAng2 = quat2eul(quatmultiply(quatRef, dqy), 'ZYX');
    eulAng3 = quat2eul(quatmultiply(quatRef, dqz), 'ZYX');
    
    
    J(:, 1) = angdiff(eulAng1, eulAngRef).' / 1e-7;
    J(:, 2) = angdiff(eulAng2, eulAngRef).' / 1e-7;
    J(:, 3) = angdiff(eulAng3, eulAngRef).' / 1e-7;

end