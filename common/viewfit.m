timeVector = (1 : size(imu, 1))*settings.dT;
% plot residuals from filter and polynomial 

Phi = settings.Phi;
stateMask = settings.stateMask;
coeffEst = X(:, stateMask.theta);
resFilter = mag' - Phi * coeffEst.';

coeffPolyFit = (Phi.'*Phi) \ (Phi.'*mag.');
resFit = mag' - Phi * coeffPolyFit;
figure;
t = tiledlayout(5, 6);
for i = 1 : 30
    nexttile,
    resF = vecnorm(resFilter(3 * (i - 1) + 1 : 3 * i, :));
    plot(timeVector, resF, 'LineWidth', 2, 'color','r');
    hold on;
    resP = vecnorm(resFit(3 * (i - 1) + 1 : 3 * i, :));
    plot(timeVector, resP, 'LineWidth', 2, 'color','g');
    grid minor;
    box on
end
title(t,'residuals')

figure;
t = tiledlayout(5, 6);
for i = 1 : 30
    nexttile,
    resF = vecnorm(resFilter(3 * (i - 1) + 1 : 3 * i, :));
    magNorm = vecnorm(mag(:, 3 * (i - 1) + 1 : 3 * i).');
    plot(timeVector, resF./magNorm * 100, 'LineWidth', 2, 'color','r');
end
title(t,'residuals (percentage)')