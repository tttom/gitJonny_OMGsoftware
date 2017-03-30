load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube','yRange','zRange');
scan2nocomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube; yRange=yRange*1e6; zRange=zRange*1e6;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_ASG0.71_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan2highint=restoredDataCube(11:end-10,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan2fullcomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;

scan2nocomp = scan2nocomp .* (scan2nocomp>0);
scan2highint = scan2highint .* (scan2highint>0);
scan2fullcomp = scan2fullcomp .* (scan2fullcomp>0);

yRange_start = 601;
yRange_end = 1732;
zRange_start = 1;
zRange_end = 451;

h=figure;
subplot(3,1,1); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan2nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,1,2); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan2highint(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compnesation - increased power');
subplot(3,1,3); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan2fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('with compensation');
saveas(h,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-highint-fullcomp2','fig');

clear scan2nocomp scan2highint scan2fullcomp;

load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube','yRange','zRange');
scan3nocomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;yRange=yRange*1e6; zRange=zRange*1e6;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan3halfcomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan3fullcomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;

scan3nocomp = scan3nocomp .* (scan3nocomp>0);
scan3halfcomp = scan3halfcomp .* (scan3halfcomp>0);
scan3fullcomp = scan3fullcomp .* (scan3fullcomp>0);

yRange_start = 799;
yRange_end = 1534;
zRange_start = 1;
zRange_end = 351;

g=figure;
%big figures
subplot(3,4,[1 2]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,[5 6]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,[9 10]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
%small figures (left)
yRange_start = 884;
yRange_end = 1025;
zRange_start = 76;
zRange_end = 201;
subplot(3,4,3); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,7); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,11); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
%small figures (right)
yRange_start = 1392;
yRange_end = 1534;
zRange_start = 1;
zRange_end = 126;
subplot(3,4,4); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,8); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,12); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan3fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
saveas(g,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\sample3_nocomp-halfcomp-fullcomp2','fig');

clear scan3nocomp scan3halfcomp scan3fullcomp;

% load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4nocomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
% load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4halfcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
% load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4fullcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
% 
% scan4nocomp = scan4nocomp .* (scan4nocomp>0);
% scan4halfcomp = scan4halfcomp .* (scan4halfcomp>0);
% scan4fullcomp = scan4fullcomp .* (scan4fullcomp>0);

% g=figure; subplot(3,2,1); imagesc(yRange,zRange,squeeze(max(scan3nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,3); imagesc(yRange,zRange,squeeze(max(scan3halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,5); imagesc(yRange,zRange,squeeze(max(scan3fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,2); imagesc(yRange,zRange,squeeze(max(scan4nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,4); imagesc(yRange,zRange,squeeze(max(scan4halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,6); imagesc(yRange,zRange,squeeze(max(scan4fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% saveas(g,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-halfcomp-fullcomp_v2','fig');
% 
% clear scan4nocomp scan4halfcomp scan4fullcomp;

load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube','yRange','zRange');
scan5nocomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;yRange=yRange*1e6; zRange=zRange*1e6;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan5halfcomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube');
scan5fullcomp=restoredDataCube(11:end-10,:,:); clear restoredDataCube;

scan5nocomp = scan5nocomp .* (scan5nocomp>0);
scan5halfcomp = scan5halfcomp .* (scan5halfcomp>0);
scan5fullcomp = scan5fullcomp .* (scan5fullcomp>0);

% m=figure; subplot(3,2,1); imagesc(yRange,zRange,squeeze(max(scan3nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,3); imagesc(yRange,zRange,squeeze(max(scan3halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,5); imagesc(yRange,zRange,squeeze(max(scan3fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,2); imagesc(yRange,zRange,squeeze(max(scan5nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,4); imagesc(yRange,zRange,squeeze(max(scan5halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% subplot(3,2,6); imagesc(yRange,zRange,squeeze(max(scan5fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
% saveas(m,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-halfcomp-fullcomp_v1','fig');

yRange_start = 742;
yRange_end = 1760;
zRange_start = 1;
zRange_end = 500;

m=figure;
%big figures
subplot(3,4,[1 2]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,[5 6]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,[9 10]); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
%small figures (left)
yRange_start = 827;
yRange_end = 968;
zRange_start = 276;
zRange_end = 476;
subplot(3,4,3); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,7); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,11); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
%small figures (right)
yRange_start = 1336;
yRange_end = 1562;
zRange_start = 76;
zRange_end = 276;
subplot(3,4,4); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5nocomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('no compensation');
subplot(3,4,8); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5halfcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('half compensation');
subplot(3,4,12); imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),squeeze(max(scan5fullcomp(:,yRange_start:yRange_end,zRange_start:zRange_end),[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');axis image;
title('full compensation');
saveas(m,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\sample5_nocomp-halfcomp-fullcomp2','fig');

clear scan5nocomp scan5halfcomp scan5fullcomp;

clear all;