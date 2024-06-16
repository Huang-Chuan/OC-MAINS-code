function [h1, h2, h3, h4] = viewimudata(timeVector, imu, mag)
%VIEWIMU plots IMU & magnetometer array data 
% timevector: Number of timestamps x 1 vector
%        imu: Number of timestamps x 6 matrix   
%             imu[:, 1:3]  accelerometer data
%             imu[:, 4:6]  gyroscope     data
%        mag: Number of timestamps x (3 x number of available magnetometers) matrix
  
  h1 = figure; clf; hold on;
  plot(timeVector, imu(:, 1), 'r-')
  plot(timeVector, imu(:, 2), 'g-')
  plot(timeVector, imu(:, 3), 'b-')
  plot(timeVector, sqrt(sum(imu(:, 1 : 3).^2, 2)), 'k-')
  title('Acceleration')

  h2 = figure; clf; hold on;
  plot(timeVector, imu(:, 4), 'r-')
  plot(timeVector, imu(:, 5), 'g-')
  plot(timeVector, imu(:, 6), 'b-')
  plot(timeVector, sqrt(sum(imu(:, 4 : 6).^2, 2)), 'k-')
  title('Angular velocity')


  h3 = figure; clf, hold on;
  numMag = size(mag, 2) / 3;
  for i = 0 : numMag - 1
    plot(timeVector, mag(:, 3*i+1), 'r-')
    plot(timeVector, mag(:, 3*i+2), 'g-')
    plot(timeVector, mag(:, 3*i+3), 'b-')
    plot(timeVector, sqrt(sum(mag(:, 3*i+(1:3)).^2, 2)), 'k-')
  end
  grid on
  title('Magnetometer readings');
  xlabel('time [s]')
  ylabel('field [\mu T]');

  h4 = figure; clf, hold on;
  t = tiledlayout(5, 6);
  for i = 0 : numMag - 1
    nexttile,
    scatter3(mag(:,3*i+1),mag(:,3*i+2),mag(:,3*i+3));
  end
  title('Magnetometer reading in 3D');


  h5 = figure;
  t = tiledlayout(5, 6);
  for i = 0 : numMag - 1
      nexttile,
      plot(timeVector, sqrt(sum(mag(:, 3*i+(1:3)).^2, 2)), 'k-')
      xlabel('time [s]')
      ylabel('field [\mu T]');
      ylim([20 80])
  end
  grid on
  

  title(t,'Magnetometic filed norm')


end



