clear all;

fileName='F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_surface\2015-08-13 16_09_37.889\recording0_lambda532nm_alpha7_beta100.mat';
load(fileName,'restoredDataCube','xRange','yRange','zRange');

% fileName='F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_surface\2015-08-13 16_09_37.889\recording0_lambda532nm_alpha0_beta100.mat';
% load(fileName,'recordedImageStack','xRange','yRange','zRange');
% restoredDataCube=recordedImageStack;
% clear recordedImageStack;

% xRange is y-dimension
% yRange is x-dimension
% zRange is actually z-dimension
restoredDataCube(1:10,:,:)=0;
restoredDataCube(end-9:end,:,:)=0;
restoredDataCube(:,1:10,:)=0;
restoredDataCube(:,end-9:end,:)=0;
restoredDataCube(:,:,1:10)=0;
restoredDataCube(:,:,end-9:end)=0;

% Gaussian filter to (smoothly) crop in-focus region of Gaussian LS
[Y,X]=meshgrid(yRange,xRange);
dataWidth=8e-6;
dataFilter=exp(-(Y/2/dataWidth).^2);

for n=1:size(restoredDataCube,3)
    restoredDataCube(:,:,n)=restoredDataCube(:,:,n).*dataFilter;
end

kxRange=([1:length(xRange)]-floor(length(xRange)/2)-1)/length(xRange)/(xRange(2)-xRange(1))*532e-9/2/0.42*2;

% figure();
for n=1:length(xRange)
    [restoredDataCube(n,:,:),rotYRange,rotZRange]=rotate2DArray(squeeze(restoredDataCube(n,:,:)),-45/360*2*pi,yRange,zRange,0,0);
%     imagesc(rotYRange,rotZRange,squeeze(restoredDataCube(n,:,:)).');axis image;
%     title(strcat('Frame: (',num2str(n),') of (',num2str(length(xRange)),')'));
%     drawnow;shg;
    disp(n)
end

figure();
imagesc(rotYRange,rotZRange,squeeze(max(restoredDataCube,[],1)).');axis image;
title('x-z projection');
drawnow;shg;


kyRange=([1:length(rotYRange)]-floor(length(rotYRange)/2)-1)/length(rotYRange)/(rotYRange(2)-rotYRange(1))*532e-9/2/0.42*2;
kzRange=([1:length(rotZRange)]-floor(length(rotZRange)/2)-1)/length(rotZRange)/(rotZRange(2)-rotZRange(1))*532e-9/2/0.42*2;

[kY,kX]=meshgrid(kyRange,kxRange);
fftFilter=exp(-(sqrt(kY.^2+kX.^2).^8/1));
fftFilterSignal=exp(-(sqrt(kY.^2+kX.^2)/1).^8).*(1-exp(-(sqrt(kY.^2+kX.^2)/0.2).^8));


totalLevel=zeros(1,size(restoredDataCube,3));
noiseReduced=zeros(1,size(restoredDataCube,3));
signalLevel=zeros(1,size(restoredDataCube,3));

fft00_01=zeros(1,size(restoredDataCube,3));
fft01_02=zeros(1,size(restoredDataCube,3));
fft02_03=zeros(1,size(restoredDataCube,3));
fft03_04=zeros(1,size(restoredDataCube,3));
fft04_05=zeros(1,size(restoredDataCube,3));
fft05_06=zeros(1,size(restoredDataCube,3));
fft06_07=zeros(1,size(restoredDataCube,3));
fft07_08=zeros(1,size(restoredDataCube,3));
fft08_09=zeros(1,size(restoredDataCube,3));
fft09_10=zeros(1,size(restoredDataCube,3));

figure(99);
for n=1:size(restoredDataCube,3)
    fftImage=fftshift(fft2(squeeze(restoredDataCube(:,:,n))));
%     filteredFFT=fftFilter.*fftImage;
%     filteredImage=ifft2(ifftshift(filteredFFT));
%     
%     filteredFFTSignal=fftFilterSignal.*fftImage;
%     filteredImageSignal=ifft2(ifftshift(filteredFFTSignal));
%     
%     totalLevel(n)=sum(sum(abs(fftImage)));
%     noiseReduced(n)=sum(sum(abs(filteredFFT)));
%     signalLevel(n)=sum(sum(abs(filteredFFTSignal)));

    fft00_01(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.0) & (sqrt(kY.^2+kX.^2)<=0.1))./sum(sum(((sqrt(kY.^2+kX.^2)>0) & (sqrt(kY.^2+kX.^2)<=0.1)))))));
    fft01_02(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.1) & (sqrt(kY.^2+kX.^2)<=0.2))/sum(sum(((sqrt(kY.^2+kX.^2)>0.1) & (sqrt(kY.^2+kX.^2)<=0.2)))))));
    fft02_03(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.2) & (sqrt(kY.^2+kX.^2)<=0.3))/sum(sum(((sqrt(kY.^2+kX.^2)>0.2) & (sqrt(kY.^2+kX.^2)<=0.3)))))));
    fft03_04(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.3) & (sqrt(kY.^2+kX.^2)<=0.4))/sum(sum(((sqrt(kY.^2+kX.^2)>0.3) & (sqrt(kY.^2+kX.^2)<=0.4)))))));
    fft04_05(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.4) & (sqrt(kY.^2+kX.^2)<=0.5))/sum(sum(((sqrt(kY.^2+kX.^2)>0.4) & (sqrt(kY.^2+kX.^2)<=0.5)))))));
    fft05_06(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.5) & (sqrt(kY.^2+kX.^2)<=0.6))/sum(sum(((sqrt(kY.^2+kX.^2)>0.5) & (sqrt(kY.^2+kX.^2)<=0.6)))))));
    fft06_07(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.6) & (sqrt(kY.^2+kX.^2)<=0.7))/sum(sum(((sqrt(kY.^2+kX.^2)>0.6) & (sqrt(kY.^2+kX.^2)<=0.7)))))));
    fft07_08(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.7) & (sqrt(kY.^2+kX.^2)<=0.8))/sum(sum(((sqrt(kY.^2+kX.^2)>0.7) & (sqrt(kY.^2+kX.^2)<=0.8)))))));
    fft08_09(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.8) & (sqrt(kY.^2+kX.^2)<=0.9))/sum(sum(((sqrt(kY.^2+kX.^2)>0.8) & (sqrt(kY.^2+kX.^2)<=0.9)))))));
    fft09_10(n)=sum(sum(abs(fftImage.*((sqrt(kY.^2+kX.^2)>0.9) & (sqrt(kY.^2+kX.^2)<=1.0))/sum(sum(((sqrt(kY.^2+kX.^2)>0.9) & (sqrt(kY.^2+kX.^2)<=1.0)))))));
    
%     subplot(4,6,[1 8]);imagesc(yRange,xRange,squeeze(restoredDataCube(:,:,n)));axis image;
%     title(strcat('base image:',num2str(n)));
%     subplot(4,6,9);imagesc(kyRange,kxRange,log(abs(fftImage)));axis image;xlim([-2 2]);ylim([-2 2]);
%     fourierScaling=get(gca,'clim');
%     
%     subplot(4,6,15);imagesc(kyRange,kxRange,fftFilter);axis image;xlim([-2 2]);ylim([-2 2]);
%     subplot(4,6,21);imagesc(kyRange,kxRange,log(abs(filteredFFT)));axis image;xlim([-2 2]);ylim([-2 2]);
%     set(gca,'clim',fourierScaling);
%     subplot(4,6,[13 20]);imagesc(yRange,xRange,abs(filteredImage));axis image;
%     title('Noise filtering');
%     
%     subplot(4,6,6);imagesc(kyRange,kxRange,fftFilterSignal);axis image;xlim([-2 2]);ylim([-2 2]);
%     subplot(4,6,12);imagesc(kyRange,kxRange,log(abs(filteredFFTSignal)));axis image;xlim([-2 2]);ylim([-2 2]);
%     set(gca,'clim',fourierScaling);
%     subplot(4,6,[4 11]);imagesc(yRange,xRange,abs(filteredImageSignal));axis image;
%     title('Signal isolation');
    
%     subplot(4,6,[16 18]);plot([1:size(restoredDataCube,3)],totalLevel,[1:size(restoredDataCube,3)],noiseReduced,[1:size(restoredDataCube,3)],signalLevel);
%     legend('total level','noise reduced','signal level');
%     title('Fourier content');
%     subplot(4,6,[22 24]);plot([1:size(restoredDataCube,3)],totalLevel*0,[1:size(restoredDataCube,3)],noiseReduced./totalLevel,[1:size(restoredDataCube,3)],signalLevel./totalLevel);
%     title('Fourier content - normalised');


%     subplot(4,6,[16 24]);plot([1:size(restoredDataCube,3)],fft00_01...
%         ,[1:size(restoredDataCube,3)],fft01_02...
%         ,[1:size(restoredDataCube,3)],fft02_03...
%         ,[1:size(restoredDataCube,3)],fft03_04...
%         ,[1:size(restoredDataCube,3)],fft04_05...
%         ,[1:size(restoredDataCube,3)],fft05_06...
%         ,[1:size(restoredDataCube,3)],fft06_07...
%         ,[1:size(restoredDataCube,3)],fft07_08...
%         ,[1:size(restoredDataCube,3)],fft08_09...
%         ,[1:size(restoredDataCube,3)],fft09_10);
%     legend('0-10%'...
%         ,'10-20%'...
%         ,'20-30%'...
%         ,'30-40%'...
%         ,'40-50%'...
%         ,'50-60%'...
%         ,'60-70%'...
%         ,'70-80%'...
%         ,'80-90%'...
%         ,'90-100%')
    disp(n)

end

    %uncomment for Airy (blue)
    subplot(2,5,1);plot(rotZRange*1e6,fft00_01,'b','Linewidth',3);title('0-10% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,2);plot(rotZRange*1e6,fft01_02,'b','Linewidth',3);title('10-20% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,3);plot(rotZRange*1e6,fft02_03,'b','Linewidth',3);title('20-30% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,4);plot(rotZRange*1e6,fft03_04,'b','Linewidth',3);title('30-40% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,5);plot(rotZRange*1e6,fft04_05,'b','Linewidth',3);title('40-50% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,6);plot(rotZRange*1e6,fft05_06,'b','Linewidth',3);title('50-60% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,7);plot(rotZRange*1e6,fft06_07,'b','Linewidth',3);title('60-70% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,8);plot(rotZRange*1e6,fft07_08,'b','Linewidth',3);title('70-80% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,9);plot(rotZRange*1e6,fft08_09,'b','Linewidth',3);title('80-90% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    subplot(2,5,10);plot(rotZRange*1e6,fft09_10,'b','Linewidth',3);title('90-100% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
    
%     %uncomment for Gaussian (red)
%     subplot(2,5,1);plot(rotZRange*1e6,fft00_01,'r','Linewidth',3);title('0-10% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,2);plot(rotZRange*1e6,fft01_02,'r','Linewidth',3);title('10-20% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,3);plot(rotZRange*1e6,fft02_03,'r','Linewidth',3);title('20-30% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,4);plot(rotZRange*1e6,fft03_04,'r','Linewidth',3);title('30-40% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,5);plot(rotZRange*1e6,fft04_05,'r','Linewidth',3);title('40-50% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,6);plot(rotZRange*1e6,fft05_06,'r','Linewidth',3);title('50-60% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,7);plot(rotZRange*1e6,fft06_07,'r','Linewidth',3);title('60-70% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,8);plot(rotZRange*1e6,fft07_08,'r','Linewidth',3);title('70-80% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,9);plot(rotZRange*1e6,fft08_09,'r','Linewidth',3);title('80-90% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);
%     subplot(2,5,10);plot(rotZRange*1e6,fft09_10,'r','Linewidth',3);title('90-100% k_m_a_x');xlabel('Depth [um]');ylabel('Magnitude [a.u.]');xlim([min(rotZRange) max(rotZRange)]*1e6);

    
    drawnow;shg;