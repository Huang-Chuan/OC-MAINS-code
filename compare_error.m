function [h1] = compare_error(in_data, out_data1, out_data2)
    % out_data1: magains output data
    % out_data2: INS output data
    
    h0=figure;
    plot(out_data1.x_h(1,:),out_data1.x_h(2,:), 'k');
    hold on
    plot(out_data2.x_h(1,:),out_data2.x_h(2,:), 'b');
    plot(in_data.gt.pos(1,:),in_data.gt.pos(2,:), 'm');
    xlim([-8 10])
    ylim([-6 6])
    grid minor;
    box on


    error_struct1 = computerror(in_data, out_data1);
    error_struct2 = computerror(in_data, out_data2);

    timeVector = error_struct1.t;






    h1=figure("Position",[320 0 930 900]);
    subplot(3,1,1)
    hold on;
    a = area(timeVector, error_struct1.stdPosErr(:, 1));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    set(gca, 'YScale', 'log')
    plot(timeVector, error_struct1.rmsePosError(:, 1), 'r','LineWidth', 2);
    plot(timeVector, error_struct2.rmsePosError(:, 1), 'k--','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('Position error (X)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
    subplot(3,1,2)
    hold on;
    a = area(timeVector, error_struct1.stdPosErr(:, 2));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    set(gca, 'YScale', 'log')
    plot(timeVector, error_struct1.rmsePosError(:, 2), 'r','LineWidth', 2);
    plot(timeVector, error_struct2.rmsePosError(:, 2), 'k--','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('Position error (Y)','FontSize',12,'FontName','Times New Roman')
    grid minor;
    box on
    subplot(3,1,3)
    hold on;
    a = area(timeVector, error_struct1.stdPosErr(:, 3));
    a.FaceColor=[.5 .5 .5];
    a.FaceAlpha = 0.3;
    set(gca, 'YScale', 'log')
    plot(timeVector, error_struct1.rmsePosError(:, 3), 'r','LineWidth', 2);
    plot(timeVector, error_struct2.rmsePosError(:, 3), 'k--','LineWidth', 2);
    xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel('error [m]','FontSize',12,'FontName','Times New Roman')
    title('Position error (Z)','FontSize',12,'FontName','Times New Roman')
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
    delta_q = delta_q(:, [3 2 1]);
    mseError = [(posEst - posRef)'.^2  ...
                 delta_q.^2 ];
    
    error_struct.rmsePosError = sqrt(mseError(:, 1:3));
    error_struct.rmseQError  =  sqrt(mseError(:, 4:6));

    error_struct.stdPosErr =   sqrt(out_data.diag_P(1:3, :))';
    error_struct.stdQError  =  sqrt(out_data.diag_P(7:9, :))';

    error_struct.t = in_data.t;
end