%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           
% Main script for comparing the MAINS algorithm and OC-MAINS algorithm on real-world dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clearvars;
rng(1);
addpath('common/');
N = 12;
%% Load filter settings
disp('Loads settings')
settings = readConfig('config/config_ocmains_realexp.json');
settings.ObservabilityConstrained = false;
%% Load data
disp('Loads data')
load(settings.dataPath);
%%
mains_stat = repmat(struct('pos_err_sq', [], 'pos_cov', [], ...
                            'yaw_err_sq', [] , 'yaw_cov', []) , 1, N);
mains_data = repmat(struct('in_data', [], 'out_data', []), 1, N);


in_data_arr = repmat(in_data, N, 1);
%%
for i = 1 : N
    in_data_arr(i).gt.init_pose = [in_data.gt.pos(:,1) + settings.init_pos_std * randn(3, 1);...
                        0.5;-0.03;0.01;
                        rotm2quat(in_data.gt.ori(:,:,1) *  SO3_exp(settings.init_q_std(1) * [0;0;randn]))'];
end

for i = 1 : N
    disp('Runs the MAG-aided INS')
    % add random noise to initial pose
    in_data_i = in_data_arr(i);
    out_data=magaidedINS(in_data_i,settings);
    [MSE, pos_uncertainty_i, yaw_uncertainty_i, x_cov, y_cov, z_cov] = compute_MSE_NEES(in_data_i, out_data);
    mains_stat(i).pos_err_sq = MSE.pos;
    mains_stat(i).x_err_sq = MSE.x;
    mains_stat(i).y_err_sq = MSE.y;
    mains_stat(i).z_err_sq = MSE.z;

    mains_stat(i).pos_cov = pos_uncertainty_i;
    mains_stat(i).x_cov = x_cov;
    mains_stat(i).y_cov = y_cov;
    mains_stat(i).z_cov = z_cov;

    mains_stat(i).yaw_err_sq = MSE.ori;
    mains_stat(i).yaw_cov = yaw_uncertainty_i;

    mains_data(i).in_data = in_data_i;
    mains_data(i).out_data = out_data;
end

stat.mains.pos_err = sqrt(mean([mains_stat.pos_err_sq],2));
stat.mains.pos_std = sqrt(mean([mains_stat.pos_cov],2));
stat.mains.yaw_err = sqrt(mean([mains_stat.yaw_err_sq],2));
stat.mains.yaw_std = sqrt(mean([mains_stat.yaw_cov],2));
%%
ocmains_stat = repmat(struct('pos_err_sq', [], 'pos_cov', [], ...
                            'yaw_err_sq', [] , 'yaw_cov', []) , 1, N);
ocmains_data = repmat(struct('in_data', [], 'out_data', []), 1, N);
settings.ObservabilityConstrained = true;

for i = 1 : N
    disp('Runs OC-MAINS')
    in_data_i = in_data_arr(i);
    out_data=magaidedINS(in_data_i,settings);
    %ploterror(in_data,OCMAINS_out_data);
    [MSE, pos_uncertainty_i, yaw_uncertainty_i, x_cov, y_cov, z_cov] = compute_MSE_NEES(in_data_i, out_data);

    ocmains_stat(i).pos_err_sq = MSE.pos;
    ocmains_stat(i).x_err_sq = MSE.x;
    ocmains_stat(i).y_err_sq = MSE.y;
    ocmains_stat(i).z_err_sq = MSE.z;

    ocmains_stat(i).pos_cov = pos_uncertainty_i;
    ocmains_stat(i).x_cov = x_cov;
    ocmains_stat(i).y_cov = y_cov;
    ocmains_stat(i).z_cov = z_cov;
    ocmains_stat(i).yaw_err_sq = MSE.ori;
    ocmains_stat(i).yaw_cov = yaw_uncertainty_i;

    ocmains_data(i).in_data = in_data;
    ocmains_data(i).out_data = out_data;
end


%% Compute the error
stat.mains.pos_err = sqrt(mean([mains_stat.pos_err_sq],2));
stat.mains.x_err = sqrt(mean([mains_stat.x_err_sq], 2));
stat.mains.y_err = sqrt(mean([mains_stat.y_err_sq], 2));
stat.mains.z_err = sqrt(mean([mains_stat.z_err_sq], 2));

stat.mains.pos_std = sqrt(mean([mains_stat.pos_cov],2));
stat.mains.x_std = sqrt(mean([mains_stat.x_cov], 2));
stat.mains.y_std = sqrt(mean([mains_stat.y_cov], 2));
stat.mains.z_std = sqrt(mean([mains_stat.z_cov], 2));

stat.mains.yaw_err = sqrt(mean([mains_stat.yaw_err_sq],2));
stat.mains.yaw_std = sqrt(mean([mains_stat.yaw_cov],2));

stat.ocmains.pos_err = sqrt(mean([ocmains_stat.pos_err_sq],2));
stat.ocmains.x_err = sqrt(mean([ocmains_stat.x_err_sq], 2));
stat.ocmains.y_err = sqrt(mean([ocmains_stat.y_err_sq], 2));
stat.ocmains.z_err = sqrt(mean([ocmains_stat.z_err_sq], 2));
stat.ocmains.pos_std = sqrt(mean([ocmains_stat.pos_cov],2));
stat.ocmains.x_std = sqrt(mean([ocmains_stat.x_cov], 2));
stat.ocmains.y_std = sqrt(mean([ocmains_stat.y_cov], 2));
stat.ocmains.z_std = sqrt(mean([ocmains_stat.z_cov], 2));
stat.ocmains.yaw_err = sqrt(mean([ocmains_stat.yaw_err_sq],2));
stat.ocmains.yaw_std = sqrt(mean([ocmains_stat.yaw_cov],2));


figure;
plot(0.01*(1:length(stat.mains.pos_err)), stat.mains.pos_err,'k');
hold on
plot(0.01*(1:length(stat.mains.pos_err)), stat.mains.pos_std,'k--','LineWidth',2);
plot(0.01*(1:length(stat.ocmains.pos_err)), stat.ocmains.pos_err,'color','#EDB120');
plot(0.01*(1:length(stat.ocmains.pos_err)), stat.ocmains.pos_std,'color','#EDB120','LineStyle','--','LineWidth',2);
legend('MAINS RMSE', 'MAINS perceived uncertainty', 'OC-MAINS RMSE', 'OC-MAINS perceived uncertainty');
ylabel('RMSE [m]','FontSize',12,'FontName','Times New Roman')
xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
title('Position RMSE and Estimator Uncertainty','FontSize',12,'FontName','Times New Roman')
grid on;
box on

figure;
plot(0.01*(1:length(stat.mains.yaw_err)), rad2deg(stat.mains.yaw_err),'k');
hold on
plot(0.01*(1:length(stat.mains.yaw_err)),  rad2deg(stat.mains.yaw_std), 'k--','LineWidth',2);
plot(0.01*(1:length(stat.ocmains.yaw_err)), rad2deg(stat.ocmains.yaw_err),'color','#EDB120');
plot(0.01*(1:length(stat.mains.yaw_err)),  rad2deg(stat.ocmains.yaw_std),'color','#EDB120','LineStyle','--','LineWidth',2);
legend('MAINS RMSE', 'MAINS perceived uncertainty', 'OC-MAINS RMSE', 'OC-MAINS perceived uncertainty');
ylabel('RMSE [$^{\circ}$]','FontSize',12,'FontName','Times New Roman', 'interpreter','latex');
xlabel('time [s]','FontSize',12,'FontName','Times New Roman')
title('Yaw RMSE and estimator perceived uncertainty','FontSize',12,'FontName','Times New Roman')
grid on;
box on

f=figure;
f.Position = [100,100,1147,889];
subplot(3,1,1)
plot(0.01*(1:length(stat.mains.x_err)), stat.mains.x_err,'k');
hold on
plot(0.01*(1:length(stat.mains.x_std)), stat.mains.x_std,'k--','LineWidth',2);
plot(0.01*(1:length(stat.ocmains.x_err)), stat.ocmains.x_err,'color','#EDB120');
plot(0.01*(1:length(stat.ocmains.x_std)), stat.ocmains.x_std,'color','#EDB120','LineStyle','--','LineWidth',2);
legend('MAINS RMSE', 'MAINS perceived uncertainty', 'OC-MAINS RMSE', 'OC-MAINS perceived uncertainty', 'FontSize', 12);
ylabel('RMSE (X) [m]','FontSize',20,'FontName','Times New Roman')
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 20); % X-axis ticks
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 20); % Y-axis ticks

title('Position RMSE and estimator perceived uncertainty','FontSize',24,'FontName','Times New Roman')
grid on;
box on
subplot(3,1,2)
plot(0.01*(1:length(stat.mains.y_err)), stat.mains.y_err,'k');
hold on
plot(0.01*(1:length(stat.mains.y_std)), stat.mains.y_std,'k--','LineWidth',2);
plot(0.01*(1:length(stat.ocmains.y_err)), stat.ocmains.y_err,'color','#EDB120');
plot(0.01*(1:length(stat.ocmains.y_std)), stat.ocmains.y_std,'color','#EDB120','LineStyle','--','LineWidth',2);
ylabel('RMSE (Y) [m]','FontSize',20,'FontName','Times New Roman')
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 20); % X-axis ticks
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 20); % Y-axis ticks
grid on;
box on
subplot(3,1,3)
plot(0.01*(1:length(stat.mains.z_err)), stat.mains.z_err,'k');
hold on
plot(0.01*(1:length(stat.mains.z_std)), stat.mains.z_std,'k--','LineWidth',2);
plot(0.01*(1:length(stat.ocmains.z_err)), stat.ocmains.z_err,'color','#EDB120');
plot(0.01*(1:length(stat.ocmains.z_std)), stat.ocmains.z_std,'color','#EDB120','LineStyle','--','LineWidth',2);
ylabel('RMSE (Z) [m]','FontSize',20,'FontName','Times New Roman')
xlabel('time [s]','FontSize',24,'FontName','Times New Roman')
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 20); % X-axis ticks
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 20); % Y-axis ticks
grid on;
box on

magdata=load('magfield\true_mag.mat');
magdata.mag_norm = magdata.mag_pred;
addpath('plotfunctions')
plot_trajectory_realexp(mains_data(1).out_data, ocmains_data(1).out_data, mains_data(1).in_data, magdata)






