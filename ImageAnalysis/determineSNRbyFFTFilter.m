%%% name:           determineSNRbyFFTFilter
%%% author:         Jonathan Nylk
%%% date created:   22/07/2015
%%% description:    This function performs low and high pass Fourier
%%%                 filtering to determine the signal-to-noise (SNR) ratio
%%%                 of an image.
%%%                 The image is high-pass filtered to remove low frequency
%%%                 noise, and low-pass filtered to remove high frequency
%%%                 noise. The "signal" is the sum of the frequency
%%%                 components remaining in this passband. The "noise" is
%%%                 the sum of all components rejected by the band-pass
%%%                 filter.
%%%
%%% updates (latest first):
%%%
%%%
%%% END %%%

function determineSNRbyFFTFilter(inputImage,normaliseImage,pixelSize,numericalAperture,lambda,lowCutOff,highCutOff,filterStrength)

    % default image
    if nargin<1
%         inputImage=imread('C:\Users\Jonathan Nylk\Downloads\NoiseImages_domar\FluoCells_SFG_pos1_z203_170759.tiff');
%         inputImage=imread('C:\Users\Jonathan Nylk\Downloads\NoiseImages_domar\FluoCells_V5_pos1_z197_140101.tiff');
%         inputImage=imread('C:\Users\Jonathan Nylk\Downloads\NoiseImages_domar\Hela_SFG_pos2_z290_154955.tiff');
        inputImage=imread('C:\Users\Jonathan Nylk\Downloads\NoiseImages_domar\Hela_V5_pos2_z289_163455.tiff');

    end
    % default values
    if nargin<2
        normaliseImage=1;
    end
    if nargin<3
        pixelSize=0.5e-6; %m
    end
    if nargin<4
        numericalAperture=1;
    end
    if nargin<5
        lambda=800e-9; %m (illumiation wavelength)
    end
    if nargin<6
        lowCutOff=0.004; %rising edge of signal bandpass filter
    end
    if nargin<7
        highCutOff=0.45; %falling edge of signal bandpass filter
    end
    if nargin<8
        filterStrength=4; %steepness of filter edge
    end
    
    inputImage=single(inputImage);
    %normalisation
    if normaliseImage
        inputImage=inputImage/max(inputImage(:));
    end
    
    %generate real-space coordinate system
    yRange=([1:size(inputImage,1)]-floor(size(inputImage,1)/2)-1)*pixelSize;
    xRange=([1:size(inputImage,2)]-floor(size(inputImage,2)/2)-1)*pixelSize;
    
    %perform Fourier transform and generate fourier-space coordinates
    fftImage=fftshift(fft2(inputImage));
    k_xRange=([1:size(inputImage,2)]-floor(size(inputImage,2)/2)-1)/(size(inputImage,2)*pixelSize/2);
    k_yRange=([1:size(inputImage,1)]-floor(size(inputImage,1)/2)-1)/(size(inputImage,2)*pixelSize/2);
    %normalise to NA and lambda
    k_xRange=k_xRange*lambda/2/numericalAperture;
    k_yRange=k_yRange*lambda/2/numericalAperture;
    %create fourier-space coordinate grid
    [k_x,k_y]=meshgrid(k_xRange,k_yRange);
    
    %generate low-frequency noise removal (high-pass), high-frequency
    %noise removal (low-pass), and signal selection (band-pass) filters.
    lowFreqNoiseFilter=exp(-((sqrt(k_x.^2+k_y.^2))/lowCutOff).^filterStrength);
    highFreqNoiseFilter=1-exp(-((sqrt(k_x.^2+k_y.^2))/highCutOff).^filterStrength);
    signalBandFilter=(1-lowFreqNoiseFilter).*(1-highFreqNoiseFilter);
    
    %apply filters
    fftSignal=fftImage.*signalBandFilter;
    fftLowFreqNoise=fftImage.*lowFreqNoiseFilter;
    fftHighFreqNoise=fftImage.*highFreqNoiseFilter;
    
    Signal=real(ifft2(ifftshift(fftSignal)));
    LowFreqNoise=real(ifft2(ifftshift(fftLowFreqNoise)));
    HighFreqNoise=real(ifft2(ifftshift(fftHighFreqNoise)));
    
    %sum Fourier magnitudes
    signalMag=sum(abs(fftSignal(:)));
    lowFreqNoiseMag=sum(abs(fftLowFreqNoise(:)));
    highFreqNoiseMag=sum(abs(fftHighFreqNoise(:)));
    
    disp(strcat('Signal magnitude = (',num2str(signalMag),')'));
    disp(strcat('Low-Freq Noise magnitude = (',num2str(lowFreqNoiseMag),')'));
    disp(strcat('High-Freq Noise magnitude = (',num2str(highFreqNoiseMag),')'));
    disp(strcat('Signal-to-Low-Freq Noise Ratio = (',num2str(signalMag/lowFreqNoiseMag),')'));
    disp(strcat('Signal-to-High-Freq Noise Ratio = (',num2str(signalMag/highFreqNoiseMag),')'));
    disp(strcat('Signal-to-Noise Ratio = (',num2str(signalMag/(lowFreqNoiseMag+highFreqNoiseMag)),')'));
    
    %plot results
    figure(1);
    %plot input image
        subplot(3,4,9);
        imagesc(xRange*1e6,yRange*1e6,inputImage);axis image;
        xlabel('x [um]');ylabel('y [um]');
        title('Input Image');
        subplot(3,4,5);
        imagesc(k_xRange,k_yRange,log10(abs(fftImage)));axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('|fft(Input Image)|');
    %plot filters
        subplot(3,4,2);
        imagesc(k_xRange,k_yRange,signalBandFilter);axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('Signal Filter');
        subplot(3,4,3);
        imagesc(k_xRange,k_yRange,lowFreqNoiseFilter);axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('Low-Freq Noise Filter');
        subplot(3,4,4);
        imagesc(k_xRange,k_yRange,highFreqNoiseFilter);axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('High-Freq Noise Filter');
    %plot filtered Fourier results
        subplot(3,4,6);
        imagesc(k_xRange,k_yRange,log10(abs(fftSignal)));axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('log_1_0(|fft(Signal)|)');
        subplot(3,4,7);
        imagesc(k_xRange,k_yRange,log10(abs(fftLowFreqNoise)));axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('log_1_0(|fft(Low-Freq Noise)|)');
        subplot(3,4,8);
        imagesc(k_xRange,k_yRange,log10(abs(fftHighFreqNoise)));axis image;
        xlabel('k_x [m^-^1*lambda/2/NA]');ylabel('k_y [m^-^1*lambda/2/NA]');
        title('log_1_0(|fft(High-Freq Noise)|)');
    %plot filtered real space results
        subplot(3,4,10);
        imagesc(xRange*1e6,yRange*1e6,Signal);axis image;
        xlabel('x [um]');ylabel('y [um]');
        title('Signal');
        subplot(3,4,11);
        imagesc(xRange*1e6,yRange*1e6,LowFreqNoise);axis image;
        xlabel('x [um]');ylabel('y [um]');
        title('Low-Freq Noise');
        subplot(3,4,12);
        imagesc(xRange*1e6,yRange*1e6,HighFreqNoise);axis image;
        xlabel('x [um]');ylabel('y [um]');
        title('High Freq Noise');
        

    figure(2);
    subplot(2,1,1);
    imagesc(xRange*1e6,yRange*1e6,inputImage);axis image;
    xlabel('x [um]');ylabel('y [um]');
    title('Input Image');
    subplot(2,1,2);
    imagesc(xRange*1e6,yRange*1e6,Signal);axis image;
    xlabel('x [um]');ylabel('y [um]');
    title('Filtered Image');
        
        
end