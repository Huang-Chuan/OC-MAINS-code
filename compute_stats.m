function [stat] = compute_stats(in_data, out_data, settings)
    % 
    stat.duration = in_data.t(end) - in_data.t(1);    
    stat.dist = sum(vecnorm(diff(in_data.gt.pos')'));
    
    

    speed_t = vecnorm(calc_vel(in_data));

    t_start = 60;
    idx = in_data.t > t_start;

    stat.duration_test = stat.duration - t_start;
    stat.dist_test = sum(vecnorm(diff(in_data.gt.pos(:, idx)')'));
    stat.height = mean(in_data.gt.pos(3, idx));
    
    posEst = out_data.x_h(1:3, idx);
    posTrue = in_data.gt.pos(:, idx);
    posCov = out_data.diag_P(1:3, idx);

    magdata = in_data.mag_array.field(idx, :);

    % horizontal error 
    stat.h_err = mean(vecnorm(posEst(1:2, :) - posTrue(1:2, :)));
    fprintf('Horizontal error: %.2f m.\n', stat.h_err);
    stat.h_err_end = mean(vecnorm(posEst(1:2, end) - posTrue(1:2, end)));
    fprintf('Horizontal error (end): %.2f m.\n',stat.h_err_end);
    % vertical error
    stat.v_err = mean(abs(posEst(3, :) - posTrue(3, :)));
    fprintf('Vertical error: %.2f m.\n', stat.v_err);
    stat.v_err_end = mean(vecnorm(posEst(3, end) - posTrue(3, end)));
    fprintf('Vertical error (end): %.2f m.\n',stat.v_err_end);
    % NEES
    pos_err = posEst - posTrue;
    nees = 0;
    for i = 1 : size(pos_err, 2)
        nees = nees + pos_err(:, i)' * (1./posCov(:, i) .* pos_err(:, i)); 
    end
    nees = nees / size(pos_err, 2);
    fprintf('NEES: %d.\n', nees);
    
    % Magnetic field info
    sigstd = 0;
    for i = 1 : size(magdata, 1)
        y = reshape(magdata(i, :)', 3, []);
        sigstd = sigstd + sqrt(mean((y-mean(y,2)).^2,'all'));
    end
    sigstd = sigstd / size(magdata, 1);
    fprintf('Magnetic field standard deviation: %f.\n', sigstd);
    
    % speed rmse
    speed_magOdometry = sqrt(sum(out_data.x_h(4:6, :).^2,1));
    speed_rmse = mean(abs(speed_t(idx) - speed_magOdometry(idx)));
    fprintf('The speed rmse is: %f.\n', speed_rmse);

    stat.speed_rmse = speed_rmse;
end