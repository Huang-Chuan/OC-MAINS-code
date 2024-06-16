function imu=load_imu_data(imuarray_data)
    % number of init samples
    nr_init_samples = 200;
    nr_IMU_nr = 32;

    imu=zeros(length(imuarray_data), 6);
    for ii=nr_IMU_nr
        imu=imu+imuarray_data(:, (ii-1)*6+(1:6));
    end
    imu=imu./length(nr_IMU_nr);

    imu(:, 1:3)=imu(:, 1:3);
    imu(:, 4:6)=imu(:, 4:6);
    imu(:, 4:6)=imu(:, 4:6)-mean(imu(1:nr_init_samples, 4:6));
end
