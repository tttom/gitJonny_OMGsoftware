load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube','yRange','zRange'); scan2nocomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube; yRange=yRange*1e6; zRange=zRange*1e6;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_ASG0.71_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan2highint=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan2\2017-03-27 17_44_32.811\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan2fullcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
scan2nocomp = scan2nocomp .* (scan2nocomp>0);
scan2highint = scan2highint .* (scan2highint>0);
scan2fullcomp = scan2fullcomp .* (scan2fullcomp>0);

h=figure; subplot(3,1,1); imagesc(yRange,zRange,squeeze(max(scan2nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,1,2); imagesc(yRange,zRange,squeeze(max(scan2highint,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,1,3); imagesc(yRange,zRange,squeeze(max(scan2fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
saveas(h,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-highint-fullcomp','fig');

clear scan2nocomp scan2highint scan2fullcomp;

load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan3nocomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan3halfcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan3\2017-03-27 17_53_13.354\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan3fullcomp=restoredDataCube; clear restoredDataCube;

scan3nocomp = scan3nocomp .* (scan3nocomp>0);
scan3halfcomp = scan3halfcomp .* (scan3halfcomp>0);
scan3fullcomp = scan3fullcomp .* (scan3fullcomp>0);

load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4nocomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4halfcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan4\2017-03-27 18_06_25.646\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan4fullcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;

scan4nocomp = scan4nocomp .* (scan4nocomp>0);
scan4halfcomp = scan4halfcomp .* (scan4halfcomp>0);
scan4fullcomp = scan4fullcomp .* (scan4fullcomp>0);

g=figure; subplot(3,2,1); imagesc(yRange,zRange,squeeze(max(scan3nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,3); imagesc(yRange,zRange,squeeze(max(scan3halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,5); imagesc(yRange,zRange,squeeze(max(scan3fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,2); imagesc(yRange,zRange,squeeze(max(scan4nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,4); imagesc(yRange,zRange,squeeze(max(scan4halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,6); imagesc(yRange,zRange,squeeze(max(scan4fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
saveas(g,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-halfcomp-fullcomp_v2','fig');

clear scan4nocomp scan4halfcomp scan4fullcomp;

load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan5nocomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan5halfcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;
load('E:\RESULTS\20170327_DFopercula_attnComp\sample2\scan5\2017-03-27 18_15_27.067\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat','restoredDataCube'); scan5fullcomp=restoredDataCube(2:end-1,:,:); clear restoredDataCube;

scan5nocomp = scan5nocomp .* (scan5nocomp>0);
scan5halfcomp = scan5halfcomp .* (scan5halfcomp>0);
scan5fullcomp = scan5fullcomp .* (scan5fullcomp>0);

m=figure; subplot(3,2,1); imagesc(yRange,zRange,squeeze(max(scan3nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,3); imagesc(yRange,zRange,squeeze(max(scan3halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,5); imagesc(yRange,zRange,squeeze(max(scan3fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,2); imagesc(yRange,zRange,squeeze(max(scan5nocomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,4); imagesc(yRange,zRange,squeeze(max(scan5halfcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
subplot(3,2,6); imagesc(yRange,zRange,squeeze(max(scan5fullcomp,[],1))'); colormap hot; xlabel('x-axis [µm]'); ylabel('z-axis [µm]');
saveas(m,'C:\Users\Jonathan\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\Deconvolution\AttenuationCompensation_WaterImmersionSystemUpdate_November2016\nocomp-halfcomp-fullcomp_v1','fig');
