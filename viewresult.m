function [] = viewresult(out_data, in_data, settings)
    timeVector = in_data.t(1:end);
    X = out_data.x_h(:, 1:end)';
    Ps = out_data.diag_P(:, 1:end).';

    stateMask = settings.stateMask;
    errorStateMask = settings.errorStateMask;
    fieldRes = out_data.res(:,(1:end));

    % extract ground truth
    if isfield(in_data.gt,'pos')
        pos_true = in_data.gt.pos(:, 1:end);
        has_pos_ref = true;
    else
        pos_true = [];
        has_pos_ref = false;
    end
    if isfield(in_data.gt,'ori')
        ori_true = in_data.gt.ori(:, :, 1:end);
        has_ori_ref = true;
    else
        ori_true = [];
        has_ori_ref = false;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%               position plot             %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    tiledlayout(3,2);
    % plot horizontal plane trajectory
    nexttile([2, 2]);
    h1=plot(X(:,1),X(:,2),'k','DisplayName','Trajectory');
    hold on;
    h2=plot(X(1,1),X(1,2),'rs','DisplayName','Start point');
    h3=plot(X(end,1),X(end,2),'r*','DisplayName','End point');
    h4=error_ellipse(diag(Ps(end,1:2)),X(end,[1 2]),'conf',0.95,'style','cyan');
    if has_pos_ref
        h5=plot(pos_true(1,:),pos_true(2,:),'b-.','DisplayName','Ground truth');
    end


    rotm = quat2rotm(X(end, 7:10));
    quiver(X(end,1), X(end,2), rotm(1, 1), rotm(2, 1), 'r--');
    quiver(X(end,1), X(end,2), rotm(1, 2), rotm(2, 2), 'b--');
    % plot final heading
    if has_ori_ref
        quiver(X(end,1), X(end,2), ori_true(1, 1, end), ori_true(2, 1, end), 'r');
        quiver(X(end,1), X(end,2), ori_true(1, 2, end), ori_true(2, 2, end), 'b');
    end
    title('2D Trajectory')
    
    if has_pos_ref
        legend([h1,h2,h3,h4,h5],{'Trajectory', 'Start point', 'End point', '95% conf.','Groud truth'})
    else
        legend([h1,h2,h3,h4],{'Trajectory', 'Start point', 'End point', '95% conf.'})
    end
    xlim([-8 8]);
    ylim([-8 8]);
    
    xlabel('x [m]')
    ylabel('y [m]')
    
    grid on
    box on
    
    % plot vertical height
    nexttile([1, 2]);
    h1=plot(timeVector, X(:, 3),'k');
    hold on
    if has_pos_ref
        h2=plot(timeVector, pos_true(3,:),'b-.');
    end
    title('Heigth')
    if has_pos_ref
        legend([h1,h2],{'Estimated height', 'True height'})
    else
        legend(h1,{'Estimated height'})
    end
    xlabel('time [s]')
    ylabel('z [m]')
    grid on
    box on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           height and speed plot         %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,1,1)
    plot(timeVector, X(:, 3),'r')
    hold on;
    curve1 =  X(:, 3)+sqrt(Ps(:, 3));
    curve2 =  X(:, 3)-sqrt(Ps(:, 3));
    fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5],'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
    %ylim([-2 2])
    title('Heigth')
    xlabel('time [s]')
    ylabel('z [m]')
    grid on
    box on

    subplot(2,1,2)
    plot(timeVector,sqrt(sum(X(:, 4:6).^2, 2)))
    title('Speed')
    xlabel('time [s]')
    ylabel('|v| [m/s]')
    grid on
    box on

    figure
    for ii=1:3
        subplot(3, 1, ii)
        plot(timeVector,X(:, 3 + ii),'r')
        if ii==1
            title('Velocity')
        end
        hold on;
        curve1 =  X(:, 3 + ii)+sqrt(Ps(:, 3 + ii));
        curve2 =  X(:, 3 + ii)-sqrt(Ps(:, 3 + ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel('Speed [m/s]')
        xlabel('time [s]')
        grid minor
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%              Orientation plot           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    labels = ['Yaw', 'Pitch', "Roll"];
    eul_deg = rad2deg(quat2eul(X(:, 7:10),'ZYX'));
    eulCov = zeros(length(timeVector), 3);
    for i = 1 : length(timeVector)
        J = ypr_jacobian_quat(X(i, 7:10));
        eulCov(i, :) = diag(J * out_data.cov_q(:, : , i) * J');
    end
    figure
    for ii=1:3
        subplot(3,1,ii)
        plot(timeVector,eul_deg(:, ii),'r')
        hold on
        if has_ori_ref
            plot(timeVector,in_data.gt.rpys(ii, :), 'k')
        end
        if ii==1
            title('Attitude')
        end
        hold on;
        % compute correct jacobian for each estimated orientation
        curve1 =  eul_deg(:, ii) + rad2deg(sqrt(eulCov(:, ii)));
        curve2 =  eul_deg(:, ii) - rad2deg(sqrt(eulCov(:, ii)));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel(strcat(labels(ii), '[deg]'))
        xlabel('time [s]')
        legend('est', 'gt')
        grid minor
    end

    % plot orientation uncertainty
    % figure;
    % labels = ["$\epsilon_x$", "$\epsilon_y$", "$\epsilon_z$"];
    % for ii=1:3
    %     subplot(3,1,ii)
    %     plot(timeVector, rad2deg(sqrt(out_data.diag_P(6+ii, :))), 'r');
    %     ylabel(strcat(labels(ii), '[deg]'),'interpreter','latex');
    %     xlabel('time [s]')
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%          Acceleration bias plot         %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    acc_bias = X(:, 11:13);
    acc_bias_var = Ps(:, 10:12);
    figure;
    for ii = 1 : 3
        subplot(3, 1, ii)
        plot(timeVector,acc_bias(:, ii),'r')
        if ii==1
            title('Accelerometer bias')
        end
        hold on;
        curve1 =  acc_bias(:, ii) + sqrt(acc_bias_var(:, ii));
        curve2 =  acc_bias(:, ii) - sqrt(acc_bias_var(:, ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        hold on
        ylabel('Bias [m/s^2]')
        xlabel('time [s]')
        grid minor
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           Gyroscope bias plot           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gyro_bias = rad2deg(X(:, 14:16));
    gyro_bias_var = Ps(:, 13:15);
    figure;
    for ii = 1 : 3
        subplot(3, 1, ii)
        plot(timeVector, gyro_bias(:, ii),'r')
        if ii==1
            title('Gyro bias')
        end
        hold on;
        curve1 =  gyro_bias(:, ii) + rad2deg(sqrt(gyro_bias_var(:, ii)));
        curve2 =  gyro_bias(:, ii) - rad2deg(sqrt(gyro_bias_var(:, ii)));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel('Bias [deg/s]')
        xlabel('time [s]')
        grid minor
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%        coefficient      plot            %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coeff = X(:, stateMask.theta);
    coeff_var = Ps(:, errorStateMask.theta);
    f=figure;
    t = tiledlayout(3, 5);
    for i = 1 : min(15, size(coeff, 2))
        nexttile,
        hold on;
        plot(timeVector, coeff(:, i), 'r', 'LineWidth', 2);
        
        curve1 =  coeff(:, i) + sqrt(coeff_var(:, ii));
        curve2 =  coeff(:, i) - sqrt(coeff_var(:, ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
    
        title(strcat('$\theta_{', num2str(i), '}$'),'FontSize',12,'FontName','Times New Roman','interpreter','latex');
        grid minor;
        
        box on
    end
    axes(f, 'visible', 'off');
    title(t, 'Coefficient');
    xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel(t, ' ','FontSize',12,'FontName','Times New Roman')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%          residual       plot            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    f=figure;
    tiledlayout(5, 6);
    for i = 1 : settings.numSensors
        nexttile,
        hold on;
        plot(timeVector, fieldRes(3*(i-1)+1, :), 'r', 'LineWidth', 2);
        plot(timeVector, fieldRes(3*(i-1)+2, :), 'g', 'LineWidth', 2);
        plot(timeVector, fieldRes(3*(i-1)+3, :), 'b', 'LineWidth', 2);
        grid minor;
        
        box on
    end
    axes(f, 'visible', 'off');
    title(t, 'Residual');
    xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel(t, ' ','FontSize',12,'FontName','Times New Roman')
end
