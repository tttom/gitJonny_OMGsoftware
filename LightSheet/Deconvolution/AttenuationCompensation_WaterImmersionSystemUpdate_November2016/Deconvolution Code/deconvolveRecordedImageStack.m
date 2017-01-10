function [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,config,showFigures,sampleAttenuation)
% recordedImageStack format: The recorded values [x y z]=[down right back] on camera, or [swipe propagation scan]
% config: a light sheet microscope set-up configuration file (.json), including the wavelengths and the pupil functions.

    cubicInterpolation=true;
    
    % Sample grid specification
    stageTranslationStepSize=norm(median(diff(config.stagePositions.target)));
    
    if isempty(config.detector.center)
        config.detector.center=[0 0];
    end
    
    if isfield(config.detection,'calibratedSystemMagnification')
        realMagnification = config.detection.calibratedSystemMagnification; % use known magnification of system if available
    else
        realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength; % calculate theoretical perfect magnification if system magnification not known
    end
    
    xRange=-config.detector.center(1)+config.detector.pixelSize(1)*([1:size(recordedImageStack,1)]-floor(size(recordedImageStack,1)/2)-1)/realMagnification; % up/down on camera
    yRange=-config.detector.center(2)+config.detector.pixelSize(2)*([1:size(recordedImageStack,2)]-floor(size(recordedImageStack,2)/2)-1)/realMagnification; % left/right on camera
    tRange=stageTranslationStepSize*([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2+1)); %Translation range (along z-axis)
    zRange=tRange; % detection axis, zero centered
    
    inversePropagation=isfield(config.excitation,'inversePropagation') && config.excitation.inversePropagation;
        
    tilt=0;
    logMessage('Calculating light sheet...');
    %Call function calcLightSheetPsf
    lightSheetPsf=calcLightSheetPsf(single(xRange),single(yRange*((-1).^inversePropagation)),single(zRange),tilt,config.excitation,config.modulation,config.sample.refractiveIndex);
    
    % Modify lightSheetPsf to account for depth-of-focus of detection lens.
    sigma=16e-6; %width of Gaussian
    GaussEnvelope=zeros(size(lightSheetPsf));
    for n=1:length(yRange)
        GaussEnvelope(1,n,:)=exp(-(zRange.^2)/2/(sigma)^2);
    end
    lightSheetPsf=lightSheetPsf.*GaussEnvelope;
    
    % normalise
    lightSheetPsf = lightSheetPsf ./ max(lightSheetPsf(:));
    
    % Implement attenuation of the light-sheet
    % simulate "perfect" absorption
    absorption_decay = repmat(exp(-sampleAttenuation * yRange),[length(zRange) 1]);
    lightSheetPsf = lightSheetPsf .* absorption_decay;
    % re-normalise
    lightSheetPsf = lightSheetPsf ./ max(lightSheetPsf(:));
    
    
    
    %Import scan shear and perspective scaling parameters from the config
    %file or manually enter the results of the
    %determineScanShearPerspectiveScalingParameters script
    if (isfield(config.detector,'scanShear'))
        scanShear=config.detector.scanShear;
    else
        scanShear=[0.9790 0.3829]*1e-3;
    end
    if (isfield(config.detector,'perspectiveScaling'))
        scaling=config.detector.perspectiveScaling;
    else
        scaling=[0.1270 0.1176]*1e-5;
    end
    
    if showFigures
        % Display XY, XZ, and YZ projections of the recorded data cube
        % (before geometric correction).
        % and the (static) light sheet PSF.
        figure();
        subplot(4,2,[1 3]);imagesc(yRange,zRange,squeeze(lightSheetPsf).');axis image;title('light sheet - recorded data (before geometric correction)');
        subplot(4,2,2);imagesc(yRange,xRange,squeeze(max(recordedImageStack,[],3)));axis image;title('xy-proj');xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(4,2,[5 7]);imagesc(yRange,zRange,squeeze(max(recordedImageStack,[],1)).');axis image;title('xz-proj');xlabel('x-axis [um]');ylabel('z-axis [um]');
        subplot(4,2,[4 8]);imagesc(xRange,zRange,squeeze(max(recordedImageStack,[],2)).');axis image;title('yz-proj');xlabel('y-axis [um]');ylabel('z-axis [um]');
        colormap gray
        figure();
        subplot(2,1,1);imagesc(yRange*1e6,xRange*1e6,squeeze(max(recordedImageStack,[],3)));axis image;title('xy-proj - before geometric correction');
        xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(2,1,2);imagesc(yRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],1)).');axis image;title('xz-proj');
        xlabel('x-axis [um]');ylabel('z-axis [um]');
        colormap gray
    end
    
    lightSheetPsfOrig=lightSheetPsf;
    if (~all(scanShear==0) || ~all(scaling==0))
        
        % Geometrically correcting recorded data cube and light sheet
        logMessage('Geometrically shifting and deforming recorded date cube and light sheet by [%0.3f%%,%0.3f%%] and a magnification change of [%0.3f%%,%0.3f%%] per micrometer...',-[scanShear scaling*1e-6]*100);
        for (zIdx=1:size(recordedImageStack,3))         
            zPos=zRange(zIdx);
            sampleXRange = scanShear(1)*zPos + xRange*(1-scaling(1)*zPos);
            sampleYRange = scanShear(2)*zPos + yRange*(1-scaling(2)*zPos);
            if (cubicInterpolation)
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'cubic',0);
            else
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'nearest',0);
            end
            recordedImageStack(:,:,zIdx)=interpolatedSlice;
            lightSheetPsf(:,:,zIdx)=interp1(yRange,lightSheetPsf(:,:,zIdx),sampleYRange,'*pchip',0);
        end
    end
 
    if showFigures
        % Display XY, XZ, and YZ projections of the recorded data cube
        % (after geometric correction).
        % and the (static) light sheet PSF.
        figure();
        subplot(4,2,[1 3]);imagesc(yRange,zRange,squeeze(lightSheetPsf).');axis image;title('light sheet - recorded data (after geometric correction)');
        subplot(4,2,2);imagesc(yRange,xRange,squeeze(max(recordedImageStack,[],3)));axis image;title('xy-proj');xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(4,2,[5 7]);imagesc(yRange,zRange,squeeze(max(recordedImageStack,[],1)).');axis image;title('xz-proj');xlabel('x-axis [um]');ylabel('z-axis [um]');
        subplot(4,2,[4 8]);imagesc(xRange,zRange,squeeze(max(recordedImageStack,[],2)).');axis image;title('yz-proj');xlabel('y-axis [um]');ylabel('z-axis [um]');
        colormap gray
        figure();
        subplot(2,1,1);imagesc(yRange*1e6,xRange*1e6,squeeze(max(recordedImageStack,[],3)));axis image;title('xy-proj - after geometric correction');
        xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(2,1,2);imagesc(yRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],1)).');axis image;title('xz-proj');
        xlabel('x-axis [um]');ylabel('z-axis [um]');
        colormap gray
    end
    
    %Call subfunction deconvolveLightSheet to perform deconvolution of the
    %image and the Airy illumination PSF.
    logMessage('Reconstructing convolved data set...');
    [recordedImageStack lightSheetDeconvFilter lightSheetOtf ZOtf tRange]=deconvolveLightSheet(xRange,yRange,zRange,tRange,recordedImageStack,config,lightSheetPsf);
    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab
    
    %Undo the geometrical deformation of the data cube
    if (~all(scanShear==0) || ~all(scaling==0))
        logMessage('Undoing the geometrically shifting and deformation of the data cube and light sheet by [%0.3f%%,%0.3f%%] and a magnification change of [%0.3f%%,%0.3f%%] per micrometer...',-[scanShear scaling*1e-6]*100);
        for (tIdx=1:size(restoredDataCube,3))
            tPos=tRange(tIdx);
            sampleXRange = scanShear(1)*tPos + xRange*(1-scaling(1)*tPos);
            sampleYRange = scanShear(2)*tPos + yRange*(1-scaling(2)*tPos);
            if (cubicInterpolation)
                interpolatedSlice=interp2(sampleYRange.',sampleXRange,restoredDataCube(:,:,tIdx),yRange.',xRange,'cubic',0);
            else
                interpolatedSlice=interp2(sampleYRange.',sampleXRange,restoredDataCube(:,:,tIdx),yRange.',xRange,'nearest',0);

            end
            restoredDataCube(:,:,tIdx)=interpolatedSlice;
        end
        lightSheetPsf=lightSheetPsfOrig;
    end  
    
    if showFigures
        % Display XY, XZ, and YZ projections of the deconvolved data cube.
        % and the (static) light sheet PSF.
        figure();
        subplot(4,2,[1 3]);imagesc(yRange,zRange,squeeze(lightSheetPsf).');axis image;title('light sheet - deconvolved data');
        subplot(4,2,2);imagesc(yRange,xRange,squeeze(max(restoredDataCube,[],3)));axis image;title('xy-proj');xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(4,2,[5 7]);imagesc(yRange,zRange,squeeze(max(restoredDataCube,[],1)).');axis image;title('xz-proj');xlabel('x-axis [um]');ylabel('z-axis [um]');
        subplot(4,2,[4 8]);imagesc(xRange,zRange,squeeze(max(restoredDataCube,[],2)).');axis image;title('yz-proj');xlabel('y-axis [um]');ylabel('z-axis [um]');
        colormap gray
        figure();
        subplot(2,1,1);imagesc(yRange*1e6,xRange*1e6,squeeze(max(restoredDataCube,[],3)));axis image;title('xy-proj - after deconvolution');
        xlabel('x-axis [um]');ylabel('y-axis [um]');
        subplot(2,1,2);imagesc(yRange*1e6,zRange*1e6,squeeze(max(restoredDataCube,[],1)).');axis image;title('xz-proj');
        xlabel('x-axis [um]');ylabel('y-axis [um]');
        colormap gray
    end
    
end

function [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf tRange]=deconvolveLightSheet(xRange,yRange,zRange,tRange,recordedImageStack,config,lightSheetPsf)
    
% Calculate an extended range over which the deconvolved image can be
    % spread out due to light sheet diffraction.
    if (length(tRange)>1)
        dt=diff(tRange([1:2]));
    else dt=1;
    end
    tRangePadded=tRange(floor(end/2)+1) + dt*([1:length(tRange)*2]-length(tRange)-1);
    
    inputSize=size(recordedImageStack);
    if (length(inputSize)<4)
        if (length(inputSize)<3)
            inputSize(3)=1;
        end
        inputSize(4)=1;
    end
    
    %(Sub-pixel) shift the PSF to the center and pad to double the size 
    if (length(zRange)>1)
        lightSheetDeconvPsf=interp1(-zRange,squeeze(lightSheetPsf(floor(end/2)+1,:,:)).',tRangePadded,'*pchip',0).';
        lightSheetDeconvPsf=permute(lightSheetDeconvPsf,[3 1 2]);
    else
        lightSheetDeconvPsf=lightSheetPsf;
        lightSheetDeconvPsf(:,:,2)=0;
    end
    deconvSize=size(lightSheetDeconvPsf); deconvSize(2)=size(recordedImageStack,2);
    
    %Deconvolve the light-sheet in Z
    
    %Fourier transform the OTF of the light sheet
    ZOtf=([1:deconvSize(3)]-floor(deconvSize(3)/2)-1)/(dt*deconvSize(3));
    
    actualNumericalAperture=config.excitation.fractionOfNumericalApertureUsed*config.excitation.objective.numericalAperture;
    excitationOpticalCutOffSpFreq=2*actualNumericalAperture/config.excitation.wavelength; % Independent of refractive index of medium
    excitationNoiseToSignalRatio=config.sample.backgroundLevel*(ZOtf/excitationOpticalCutOffSpFreq)/config.sample.signalLevel;

    lightSheetOtf=fftshift(fft(ifftshift(lightSheetDeconvPsf(1,:,:),3),[],3),3);
    clear lightSheetDeconvPsf;
    lightSheetOtf=lightSheetOtf./max(abs(lightSheetOtf(:))); % Maintain the mean brightness
    
    % Scale OTF for attenuation compensation parameters
    % Check for attenuation compensation parameters
        if isfield(config.modulation,'sigmaV')
            sigmaV = config.modulation.sigmaV;
        else
            sigmaV = 0;
        end
    compensationScalingParameter = 10.478 * (sigmaV.^2) + (1.8838 * sigmaV) + 1;
    lightSheetOtf = lightSheetOtf * compensationScalingParameter;
    
    
    %Construct the deconvolution filter
    lightSheetDeconvFilter=conj(lightSheetOtf)./(abs(lightSheetOtf).^2+repmat(permute(excitationNoiseToSignalRatio.^2,[1 3 2]),[size(lightSheetOtf,1) size(lightSheetOtf,2) 1]));
    
    %Extend edges to double the size and convolve by multiplying the
    %deconvolution filter with the fft of each slice of the image stack,
    %then ifft the product.
    for (xIdx=1:size(recordedImageStack,1))
        recordedImageStackSliceFft=fft(recordedImageStack(xIdx,:,[1:end, end*ones(1,floor(end/2)), ones(1,floor((end+1)/2))],:),[],3);
        restoredDataCubeSlice=ifft(recordedImageStackSliceFft.*repmat(ifftshift(lightSheetDeconvFilter,3),[1 1 1 inputSize(4:end)]),[],3,'symmetric');

        %Drop the bit that might have overlap from a wrapped kernel
        recordedImageStack(xIdx,:,:,:)=restoredDataCubeSlice(:,:,1:end/2,:); % overwrite recordedImageStack to save memory
    end
    
    %Output
    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab
end
