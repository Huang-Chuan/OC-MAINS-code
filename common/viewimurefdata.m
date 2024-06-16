function [h1, h2, h3] = viewimurefdata(timeVector, imu, ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    numberOfSamples = size(imu, 1);
    Rn_b = ref.rotations;
    imu = imu.';

    specificForce = zeros(numberOfSamples, 3);
 
    for i = 1 : numberOfSamples
        specificForce(i, :) = (squeeze(Rn_b(:, :, i)) * imu(1:3, i))';
    end
    
    eul = rad2deg(unwrap(rotm2eul(Rn_b, 'XYZ')));
    

    h1=figure;
    plot(timeVector, specificForce);
    legend('Forward', 'Right', 'Down');

    h2=figure;
    plot(timeVector, eul);
    legend('Roll', 'Pitch', 'Yaw');    

    h3=figure;
    t = tiledlayout(3, 1);
    nexttile;
    plot(timeVector, ref.positions(1, :), 'r', 'LineWidth', 2);
    nexttile;
    plot(timeVector, ref.positions(2, :), 'r', 'LineWidth', 2);
    nexttile;
    plot(timeVector, ref.positions(3, :), 'r', 'LineWidth', 2);
    title(t, 'trajectory')
end