function [h1, h2] = ploterror(in_data, out_data)
    
    error_struct = computerror(in_data, out_data);
    timeVector = error_struct.t;

    h1=figure;
    subplot(3,1,1)
    hold on;
    a = area(timeVector, error_struct.stdPosErr(:, 1));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    %set(gca, 'YScale', 'log')
    plot(timeVector, error_struct.rmsePosError(:, 1), 'r','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('position error (X)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
    subplot(3,1,2)
    hold on;
    a = area(timeVector, error_struct.stdPosErr(:, 2));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    %set(gca, 'YScale', 'log')
    plot(timeVector, error_struct.rmsePosError(:, 2), 'r','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('position error (Y)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
    subplot(3,1,3)
    hold on;
    a = area(timeVector, error_struct.stdPosErr(:, 3));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    %set(gca, 'YScale', 'log')
    plot(timeVector, error_struct.rmsePosError(:, 3), 'r','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('position error (Z)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             orientation plot            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eulCov = zeros(in_data.numFrames, 3);
    for i = 1 : in_data.numFrames
        J = ypr_jacobian_quat(out_data.x_h(7:10, i)');
        eulCov(i, :) = diag(J * out_data.cov_q(:, :, i) * J');
    end


    h2=figure;
    subplot(3,1,1)
    hold on;
    curve1 =  rad2deg(sqrt(eulCov(:, 1)));
    curve2 =  zeros(size(error_struct.stdQError(:, 1)));
    inBetween = [curve1; flipud(curve2)];
    tt = [timeVector fliplr(timeVector)];
    h = fill(tt, inBetween, [.5 .5 .5]);
    h.FaceAlpha = 0.8;
    plot(timeVector,error_struct.rmseQError(:, 1)  * 180 / pi, 'r','LineWidth', 2, 'DisplayName','Upd');
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [$^{\circ}$]','FontSize',12,'FontName','Times New Roman', 'interpreter','latex');
    title('Orientation error (Yaw)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
    %

    subplot(3,1,2)
    hold on;
    curve1 =  rad2deg(sqrt(eulCov(:, 2)));
    curve2 =  zeros(size(error_struct.stdQError(:, 2)));
    inBetween = [curve1; flipud(curve2)];
    tt = [timeVector fliplr(timeVector)];
    h = fill(tt, inBetween, [.5 .5 .5]);
    h.FaceAlpha = 0.8;
    plot(timeVector, error_struct.rmseQError(:, 2)  * 180 / pi, 'r','LineWidth', 2, 'DisplayName','Upd');
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [$^{\circ}$]','FontSize',12,'FontName','Times New Roman', 'interpreter','latex');
    title('Orientation error (Pitch)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on



    subplot(3,1,3)
    hold on;
    curve1 =  rad2deg(sqrt(eulCov(:, 3)));
    curve2 =  zeros(size(error_struct.stdQError(:, 3)));
    inBetween = [curve1; flipud(curve2)];
    tt = [timeVector fliplr(timeVector)];
    h = fill(tt, inBetween, [.5 .5 .5]);
    h.FaceAlpha = 0.8;
    plot(timeVector, error_struct.rmseQError(:, 3)* 180 / pi, 'r','LineWidth', 2, 'DisplayName','Upd');
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [$^{\circ}$]','FontSize',12,'FontName','Times New Roman', 'interpreter','latex');
    title('Orientation error (Roll)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
end


function error_struct = computerror(in_data, out_data)
    
    posEst = out_data.x_h(1:3, :);
    oriEst = out_data.x_h(7:10, :);

    posRef = in_data.gt.pos;
    oriRef = in_data.gt.ori;

    %ref_Pos = ref.positions';
    delta_q = angdiff(quat2eul(oriEst','ZYX'), rotm2eul(oriRef,'ZYX'));
    mseError = [(posEst - posRef)'.^2  ...
                 delta_q.^2 ];
    
    error_struct.rmsePosError = sqrt(mseError(:, 1:3));
    error_struct.rmseQError  =  sqrt(mseError(:, 4:6));

    error_struct.stdPosErr =   sqrt(out_data.diag_P(1:3, :))';
    error_struct.stdQError  =  sqrt(out_data.diag_P(7:9, :))';

    error_struct.t = in_data.t;
end