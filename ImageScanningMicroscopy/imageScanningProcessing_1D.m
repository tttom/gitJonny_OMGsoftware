%%% function is work in progress, currently for testing and simulations,
%%% will be fully usable for real data processing when finished. Full
%%% documentation as function evolves.

function imageScanningProcessing_1D

    lambda = 532e-9; % [metres]
    NA = 1.49;
    
    % Define beam/stage scan coordinates
    rDelta = 10e-9; % scan step size [metres]
    rRange = [-5e-6:rDelta:5e-6] ; % scan coordinates [metres]
    
    % Define sample coordinates in camera plane
    cameraPixelSize = 6.5e-6; % pixel size/pitch [metres]
    cameraPixels = 51; % number of pixels to image (always make odd)
    cameraMagnification = 100; % magnification of optical system to camera
    sDelta = cameraPixelSize / cameraMagnification; % sample step size [metres]
    sRange = ([1:cameraPixels] - ceil(cameraPixels / 2)) * sDelta; % sample CCD coords [metres]

    % Define target to image
    sampleTarget = zeros(size(rRange));
    sampleTarget(1,ceil(length(rRange) / 2) - 1:ceil(length(rRange) / 2) + 1) = 1;
    sampleTarget(1,ceil(length(rRange) / 4) - 1:ceil(length(rRange) / 4) + 1) = 1;
    sampleTarget(1,3 * ceil(length(rRange) / 4) - 1:3 * ceil(length(rRange) / 4) + 1) = 1;
    sampleTarget(1,3 * ceil(length(rRange) / 4) + floor(300e-9/rDelta) - 1:3 * ceil(length(rRange) / 4) + floor(300e-9/rDelta) + 1) = 1;
    
    % Define illumination/detection PSF
    illuminationPSF = exp(-1 .* rRange.^2 * NA.^2 / 2 / 0.42.^2 / lambda.^2);
    illuminationOTF = fft(ifftshift(illuminationPSF));
    illuminationOTF = illuminationOTF ./ illuminationOTF(1);
    detectionPSF = illuminationPSF; % same as illumination as through same objective
    detectionOTF = illuminationOTF; % same as illumination as through same objective
    
    % Simulate confocal imaging
    fluorescenceIntensity = squeeze(zeros([size(rRange,2) 1]));
    detectedImage = squeeze(zeros([size(rRange,2) 1]));
    imageArray = zeros([size(rRange,2) size(sRange,2)]);
    for rIdx = 1:length(rRange)
        scanningIlluminationPSF = exp(-1 .* (rRange - rRange(rIdx)).^2 * NA.^2 / 2 / 0.42.^2 / lambda.^2);
        fluorescenceIntensity = scanningIlluminationPSF .* sampleTarget;
        detectedImage = ifft(fft(fluorescenceIntensity) .* detectionOTF,'symmetric');
        imageArray(rIdx,:) = interp1(squeeze(rRange),squeeze(detectedImage),squeeze(rRange(rIdx) - sRange),'pchip',0);
    end
    
    % Plot target
    figure();
    plot(rRange * 1e6,sampleTarget);
    xlabel('rRange [um]');
    ylabel('Intensity [a.u.]');
    title('Target');
    
    % Plot illumination PSF
    figure();
    plot(rRange * 1e6,illuminationPSF);
    xlabel('rRange [um]');
    ylabel('Intensity [a.u.]');
    title('Illumination/Detection PSF');
    
    % Plot image scanning matrix
    figure();
    imagesc(sRange * 1e6,rRange * 1e6,imageArray);
    xlabel('sRange [um]');
    ylabel('rRange [um]');
    title('Image scanning matrix');

    confocalImage_smallPinhole = imageArray(:,ceil(cameraPixels / 2));
    confocalImage_largePinhole = sum(imageArray,2);
    confocal_mediumPinholeHalfSize = 2; % [pixels]
    confocalImage_mediumPinhole = sum(imageArray(:,ceil(cameraPixels / 2) - confocal_mediumPinholeHalfSize:ceil(cameraPixels / 2) + confocal_mediumPinholeHalfSize),2);

    % Plot confocal images
    % Plot target
    figure();
    plot(rRange * 1e6,sampleTarget,rRange * 1e6,confocalImage_smallPinhole,rRange * 1e6,confocalImage_mediumPinhole,rRange * 1e6,confocalImage_largePinhole);
    xlabel('rRange [um]');
    ylabel('Intensity [a.u.]');
    title('Confocal images (non-normalised)');

    % Plot confocal images
    % Plot target
    figure();
    plot(rRange * 1e6,sampleTarget,rRange * 1e6,confocalImage_smallPinhole / max(confocalImage_smallPinhole(:)),rRange * 1e6,confocalImage_mediumPinhole / max(confocalImage_mediumPinhole(:)),rRange * 1e6,confocalImage_largePinhole / max(confocalImage_largePinhole(:)));
    xlabel('rRange [um]');
    ylabel('Intensity [a.u.]');
    title('Confocal images (normalised)');


end