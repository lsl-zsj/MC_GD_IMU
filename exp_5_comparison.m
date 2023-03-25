

clear all
addpath(genpath(pwd));


load('gait_0.5_magd.mat');
IMU=gait;

fs=IMU.fs;
sample_freq=fs;

Accelerometer=-IMU.Acceleration;
Gyroscope=IMU.Gyroscope;
Magnetic=IMU.Magnetic*100;
len=length(Accelerometer);

% plot the raw data
for i=1:len
    MagNorm(i)=norm(Magnetic(i,:));
end
t=0:1/fs:1/fs*(length(Accelerometer)-1);
time=[t;t;t];
time=time';

figure
x1=subplot(3,1,1);
plot(time,Accelerometer)
legend('x','y','z','interpreter','latex')
ylabel('$y_{acc}$','interpreter','latex')
set(gca,'FontSize',16)
x2=subplot(3,1,2);
plot(time,Gyroscope)
legend('x','y','z','interpreter','latex')
ylabel('$y_{gyr}$','interpreter','latex')
set(gca,'FontSize',16)
x3=subplot(3,1,3);
plot(time,Magnetic,time,MagNorm,'black')
legend('x','y','z','norm','interpreter','latex')
ylabel('$y_{mag}$','interpreter','latex')
xlabel('$time/s$','interpreter','latex')
set(gca,'FontSize',16)
linkaxes([x1,x2,x3],'x')




%% ESKF
MagSth=80;
ahrs=orientation_estimation_ahrs_fun_xsens(Accelerometer,Gyroscope,Magnetic,fs,MagSth);
euler_eskf=eulerd(ahrs.Quat,'ZXY','frame');

%% MKMC
% sigma_1=2.01;
% sigma_2=0.1351;
sigma_1=1.6188;
sigma_2=0.4234;
sigma1=2*sigma_1*sigma_1;
sigma2=2*sigma_2*sigma_2;
xigma_x=[10^8 10^8 10^8 10^8 10^8 10^8 sigma1 sigma1 sigma1 sigma2 sigma2 sigma2]; 
xigma_y=[10^8 10^8 10^8 10^8 10^8 10^8];
mkmc_ahrs=orientation_estimation_ahrs_mkmc_fun_xsens(Accelerometer,Gyroscope,Magnetic,fs,xigma_x,xigma_y,MagSth);
euler_mkmc=eulerd(mkmc_ahrs.Quat,'ZXY','frame');


%% MadgwickAHRS
AHRS = MadgwickAHRS('SamplePeriod', 1/fs, 'Beta', 0.1);

time=0:1/fs:1/fs*(len-1);
quat = zeros(length(time), 4);
for t = 1:length(time)
    AHRS.Update(Gyroscope(t,:), Accelerometer(t,:), Magnetic(t,:));	% gyroscope units must be radians
    quat(t, :) = AHRS.Quaternion;
end

% Plot algorithm output as Euler angles
for i=1:length(quat)
Quat(i)=quaternion(quat(i,1),quat(i,2),quat(i,3),quat(i,4));
end

euler_gd=eulerd(Quat,'ZXY','frame');


%% ShileiLiAHRS
Sigma_acc=0.1;
Sigma_mag=0.05;
Sigma_acc=0.0760;
Sigma_mag=0.015;
AHRS = ShileiLiAHRS('SamplePeriod', 1/fs, 'Beta', 0.1,'Sigma_g',Sigma_acc,'Sigma_m',Sigma_mag);

time=0:1/fs:1/fs*(len-1);
quat = zeros(length(time), 4);
Err = zeros(length(time), 6);
for t = 1:length(time)
    AHRS.Update(Gyroscope(t,:), Accelerometer(t,:), Magnetic(t,:));	% gyroscope units must be radians
    quat(t, :) = AHRS.Quaternion;
    Err(t,:)=AHRS.Err;
end

% Plot algorithm output as Euler angles
for i=1:length(quat)
Quat(i)=quaternion(quat(i,1),quat(i,2),quat(i,3),quat(i,4));
end

euler_cgd=eulerd(Quat,'ZXY','frame');

% figure
% hold on
% plot(euler_gd)
% plot(euler_cgd)
% legend('gd x','gd y','gd z','cgd x','cgd y','cgd z')

%% DOE
tauAcc= 1;
tauMag= 5;
zeta= 0.1;
accRating= 0;

out =EMDI(Accelerometer,Gyroscope,Magnetic,sample_freq, tauAcc, tauMag, zeta, accRating);
quat_doe=out.q;
for i=1:length(quat_doe)
Quat_doe(i)=quaternion(quat_doe(i,:));
end
euler_doe=eulerd(Quat_doe,'ZXY','frame');


%% CDOE
Sigma_acc= 0.0591 ;   % rad_acc / sample_freq is the maximum rotation angle in a sampling intervel, 0.05 is the gyroscope drift percentage 
Sigma_mag= 0.1092; 
Sigma_mag= 0.05; 

outc =EMDI_C(Accelerometer,Gyroscope,Magnetic,sample_freq, tauAcc, tauMag, zeta, accRating,Sigma_acc,Sigma_mag);
quat_emdc=outc.q;
for i=1:length(quat_emdc)
Quat_EMDC(i)=quaternion(quat_emdc(i,:));
end
euler_cdoe=eulerd(Quat_EMDC,'ZXY','frame');

%%
data.imu=IMU;
data.euler_eskf=euler_eskf-mean(euler_eskf(1:2000,:));
data.euler_mkmc=euler_mkmc-mean(euler_mkmc(1:2000,:));
data.euler_gd=euler_gd-mean(euler_gd(1:2000,:));
data.euler_cgd=euler_cgd-mean(euler_cgd(2:1000,:));
data.euler_doe=euler_doe-mean(euler_doe(2:1000,:));
data.euler_cdoe=euler_cdoe-mean(euler_cdoe(2:1000,:));

%% 
enc=load('enc_gait_0.5_magd.mat');
t_s=24.181+7.130; COR=[53.94,0.9274,-91.82];
Ang.t=enc.gait(:,17)-enc.gait(1,17)+t_s;
Ang.ang=enc.gait(:,1)+enc.gait(:,2);
Ang.ang=-(Ang.ang-Ang.ang(1))/pi*180;
encindex=find(Ang.t<time(end));
Ang.t=Ang.t(encindex);
Ang.ang=Ang.ang(encindex);


TrueAng=[];ZeroM=[];
TrueAng=zeros(length(Ang.t),3);
TrueAng(:,2)=TrueAng(:,2)+Ang.ang;

ZeroM=zeros(round(t_s*1000),3);
TrueAng=[ZeroM;TrueAng];
TrueT=0:0.001:Ang.t(end);
TrueT=TrueT';

TrueT_fs=0:1/fs:TrueT(end);
lenData=length(TrueT_fs);
TrueAng_fs=zeros(lenData,3);
TrueAng_fs(:,1)=spline(TrueT, TrueAng(:,1),TrueT_fs);
TrueAng_fs(:,2)=spline(TrueT, TrueAng(:,2),TrueT_fs);
TrueAng_fs(:,3)=spline(TrueT, TrueAng(:,3),TrueT_fs);


figure
plot(Ang.t,Ang.ang)
hold on
plot(time,data.euler_mkmc(:,2),'LineWidth',0.8,'Color','red','LineStyle','--')
plot(time,data.euler_cgd(:,2),'LineWidth',0.8,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
plot(time,data.euler_cdoe(:,2),'LineWidth',0.8,'Color','black','LineStyle','-.')
legend('enc','MKMC','CGD','CDOE')
%% error calculation

data.euler_eskf(lenData+1:end,:)=[];
data.euler_mkmc(lenData+1:end,:)=[];
data.euler_gd(lenData+1:end,:)=[];
data.euler_cgd(lenData+1:end,:)=[];
data.euler_doe(lenData+1:end,:)=[];
data.euler_cdoe(lenData+1:end,:)=[];

data.t=TrueT_fs;
data.TrueAng_fs=TrueAng_fs;
data.euler_eskf_e=data.euler_eskf-TrueAng_fs;
data.euler_mkmc_e=data.euler_mkmc-TrueAng_fs;
data.euler_gd_e=data.euler_gd-TrueAng_fs;
data.euler_cgd_e=data.euler_cgd-TrueAng_fs;
data.euler_doe_e=data.euler_doe-TrueAng_fs;
data.euler_cdoe_e=data.euler_cdoe-TrueAng_fs;


data.euler_eskf_e_rmse=rms(data.euler_eskf_e);
data.euler_gd_e_rmse=rms(data.euler_gd_e);
data.euler_doe_e_rmse=rms(data.euler_doe_e);
data.euler_mkmc_e_rmse=rms(data.euler_mkmc_e);
data.euler_cgd_e_rmse=rms(data.euler_cgd_e);
data.euler_cdoe_e_rmse=rms(data.euler_cdoe_e);



data.euler_eskf_e_me=max(abs(data.euler_eskf_e));
data.euler_gd_e_me=max(abs(data.euler_gd_e));
data.euler_doe_e_me=max(abs(data.euler_doe_e));
data.euler_mkmc_e_me=max(abs(data.euler_mkmc_e));
data.euler_cgd_e_me=max(abs(data.euler_cgd_e));
data.euler_cdoe_e_me=max(abs(data.euler_cdoe_e));

%% 
md_s=12.1; md_e=52.7;
ad_s=31.47; ad_e=46.68;

figure
x1=subplot(3,1,1);
hold on
box on
plot(data.t,data.euler_eskf_e(:,1),'LineWidth',0.5,'Color','red','LineStyle','-')
plot(data.t,data.euler_gd_e(:,1),'LineWidth',0.5,'Color','blue','LineStyle','-')
plot(data.t,data.euler_doe_e(:,1),'LineWidth',0.5,'Color','black','LineStyle','-')
plot(data.t,data.euler_mkmc_e(:,1),'LineWidth',0.5,'Color','red','LineStyle','--')
plot(data.t,data.euler_cgd_e(:,1),'LineWidth',0.5,'Color','blue','LineStyle','--')
plot(data.t,data.euler_cdoe_e(:,1),'LineWidth',0.5,'Color','black','LineStyle','-.')
ylabel('yaw error','interpreter','latex')
x= [md_s ad_s ad_s md_s];
y= [-29 -29 3.5 3.5];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
hold on
xx= [ad_s ad_e ad_e ad_s];
yy= [-29 -29 3.5 3.5];
patch(xx,yy,[0.9290 0.6940 0.1250],'LineStyle','none','FaceAlpha',0.2)
x= [ad_e md_e md_e ad_e];
y= [-29 -29 3.5 3.5];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
ylabel('yaw error ($\deg$)','interpreter','latex')
xticks([])
ylim([-29 3.5])
set(gca,'FontSize',14)
x2=subplot(3,1,2);
hold on
box on
plot(data.t,data.euler_eskf_e(:,2),'LineWidth',0.5,'Color','red','LineStyle','-')
plot(data.t,data.euler_gd_e(:,2),'LineWidth',0.5,'Color','blue','LineStyle','-')
plot(data.t,data.euler_doe_e(:,2),'LineWidth',0.5,'Color','black','LineStyle','-')
plot(data.t,data.euler_mkmc_e(:,2),'LineWidth',0.5,'Color','red','LineStyle','--')
plot(data.t,data.euler_cgd_e(:,2),'LineWidth',0.5,'Color','blue','LineStyle','--')
plot(data.t,data.euler_cdoe_e(:,2),'LineWidth',0.5,'Color','black','LineStyle','-.')
x= [md_s ad_s ad_s md_s];
y= [-8.5 -8.5 6.5 6.5];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
hold on
xx= [ad_s ad_e ad_e ad_s];
yy= [-8.5 -8.5 6.5 6.5];
patch(xx,yy,[0.9290 0.6940 0.1250],'LineStyle','none','FaceAlpha',0.2)
x= [ad_e md_e md_e ad_e];
y= [-8.5 -8.5 6.5 6.5];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
xticks([])
ylim([-8.5 6.5])
ylabel('roll error ($\deg$)','interpreter','latex')
set(gca,'FontSize',14)
box on

x3=subplot(3,1,3);
hold on
box on
plot(data.t,data.euler_eskf_e(:,3),'LineWidth',0.5,'Color','red','LineStyle','-')
plot(data.t,data.euler_gd_e(:,3),'LineWidth',0.5,'Color','blue','LineStyle','-')
plot(data.t,data.euler_doe_e(:,3),'LineWidth',0.5,'Color','black','LineStyle','-')
plot(data.t,data.euler_mkmc_e(:,3),'LineWidth',0.5,'Color','red','LineStyle','--')
plot(data.t,data.euler_cgd_e(:,3),'LineWidth',0.5,'Color','blue','LineStyle','--')
plot(data.t,data.euler_cdoe_e(:,3),'LineWidth',0.5,'Color','black','LineStyle','-.')
x= [md_s ad_s ad_s md_s];
y= [-5 -5 4 4];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
hold on
xx= [ad_s ad_e ad_e ad_s];
yy= [-5 -5 4 4];
patch(xx,yy,[0.9290 0.6940 0.1250],'LineStyle','none','FaceAlpha',0.2)
x= [ad_e md_e md_e ad_e];
y= [-5 -5 4 4];
patch(x,y,[0.3010 0.7450 0.9330],'LineStyle','none','FaceAlpha',0.2)
ylim([-5 4])
lg=legend('ESKF','GD','DOE','MKMC','CGD','CDOE','$M_d$','$A_d,M_d$','interpreter','latex','Orientation','horizontal');
set(lg,'FontSize',10)
ylabel('pitch error ($\deg$)','interpreter','latex')
set(gca,'FontSize',14)
linkaxes([x1,x2,x3],'x')
xlabel('tims (s)','interpreter','latex')
xlim([10 65])
set(gcf,'position',[100 100 750 600])




