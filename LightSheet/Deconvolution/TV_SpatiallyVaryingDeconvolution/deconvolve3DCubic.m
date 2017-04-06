%[restoredDataCube, xRange,yRange,zRange, tRange, lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,config)
%
% recordedImageStack: The recorded values [x y z]=[down right back] on camera, or [swipe propagation scan]
% config: a light sheet microscope set-up configuration file, including the wavelengths and the pupil functions. 
%
function [restoredDataCube, xRange,yRange,zRange, tRange, lightSheetPsf] = deconvolve3DCubic(recordedImageStack,config)
   cubicInterpolation=true;
    
    % Sample grid specification
    stageTranslationStepSize=norm(median(diff(config.stagePositions.target)));
    
    if isempty(config.detector.center)     %%% If in the config file the detector's centre isn't defined then set it to [0 0]
        config.detector.center=[0 0];
    end
    
         %%% Defining the magnification, and ranges of axis using
         %%% information from the config file 

    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    pixelPitchInSample=[config.detector.pixelSize./realMagnification stageTranslationStepSize];     
    xRange=-config.detector.center(1)+config.detector.pixelSize(1)*([1:size(recordedImageStack,1)]-floor(size(recordedImageStack,1)/2)-1)/realMagnification; % up/down on camera
    yRange=-config.detector.center(2)+config.detector.pixelSize(2)*([1:size(recordedImageStack,2)]-floor(size(recordedImageStack,2)/2)-1)/realMagnification; % left/right on camera
    tRange=stageTranslationStepSize*([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2+1)); %Translation range (along z-axis)
    zRange=tRange; % detection axis, zero centered
    
    inversePropagation=isfield(config.excitation,'inversePropagation') && config.excitation.inversePropagation;
    
    %%% Here is where the Light Sheet's PSF is calcuated in another function
    
    tilt=0;
    logMessage('Calculating light sheet...');
    lightSheetPsf=calcLightSheetPsf(single(xRange),single(yRange*((-1).^inversePropagation)),single(zRange),tilt,config.excitation,config.modulation,config.sample.refractiveIndex);
         
    %%% Shear & Scale parameters are searched for in the detctor config file  if
    %%% not there they are given default values
    
    if (isfield(config.detector,'scanShear'))
        scanShear=config.detector.scanShear;
    else
        scanShear=[0 0];
    end
    if (isfield(config.detector,'perspectiveScaling'))
        scaling=config.detector.perspectiveScaling;
    else
        scaling=[0 0];
    end
    
    if (~all(scanShear==0) || ~all(scaling==0))     %%% checking that there are no scan scaling values equal to zero
        % Geometrically correcting recorded data cube and light sheet
        logMessage('Geometrically shifting and deforming recorded data cube and light sheet by [%0.3f%%,%0.3f%%] and a magnification change of [%0.3f%%,%0.3f%%] per micrometer...',-[scanShear scaling*1e-6]*100);
%         YRangeMatrix=repmat(yRange,[size(xRange,2) 1]);XRangeMatrix=repmat(xRange.',[1 size(yRange,2)]);

        for (zIdx=1:size(recordedImageStack,3))   %%% Going through each slice in the z plane
%             [~, sV] = memory();
%             logMessage('zIdx=%d, memory available: %0.0f MB',[zIdx sV.PhysicalMemory.Available/2^20]);            
            zPos=zRange(zIdx);
            sampleXRange = scanShear(1)*zPos + xRange*(1-scaling(1)*zPos);      %%% This sorts any skew in the tails of the airy beam and straightens them
            sampleYRange = scanShear(2)*zPos + yRange*(1-scaling(2)*zPos);
            if (cubicInterpolation)
%                 interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'*cubic',0);
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'nearest',0);   %%% Creating an interpolated image of the slice to the specified range required
            else
%                 interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'*linear',0);
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'nearest',0);
%                 interpolatedSlice=qinterp2(YRangeMatrix,XRangeMatrix,recordedImageStack(:,:,zIdx),repmat(sampleYRange,[size(sampleXRange,2) 1]),repmat(sampleXRange.',[1 size(sampleYRange,2)]), 2);
%                 interpolatedSlice(isnan(interpolatedSlice))=0;
            end
            recordedImageStack(:,:,zIdx)=interpolatedSlice;    %%% Then saving the interpolated slices into the stack
            lightSheetPsf(:,:,zIdx)=interp1(yRange,lightSheetPsf(:,:,zIdx),sampleYRange,'*pchip',0);   %%% Saving the 1D interpolation as the PSF
        end
    end
      
    %%% This is the line of code that takes the data to the
    %%% deconvolution function at the bottom of this tab
        
    
    
 % recordedImageStack = deconvolveCubic(recordedImageStack,config.detector.pixelSize./config.detection.objective.magnification,config.excitation.wavelength,...
                  %  config.detection.objective.numericalAperture,config.sample.refractiveIndex,config.detection.modulation.alpha,config.detection.modulation.beta);  
    
    
    kSNR=config.sample.signalLevel/config.sample.backgroundLevel;
        
    detectionPupilFunctor=@(U,V) exp(config.detection.modulation.alpha*2i*pi*(U.^3+V.^3)); % cubic
    %detectionPupilFunctor=@(U,V) exp(alpha*2i*pi*(U.^3+V.^3-3*(U.^2.*V+U.*V.^2))); % gen. cubic
    
    % define the sample grid
    imgSize=size(recordedImageStack);
    xRange=pixelPitchInSample(1)*([1:imgSize(1)]-floor(imgSize(1)./2)-1);
    yRange=pixelPitchInSample(2)*([1:imgSize(2)]-floor(imgSize(2)./2)-1);
    
    % determine the PSF and the corresponding Wiener filter
    detectionPsf=calcVectorialPsf(xRange,yRange,zRange,config.detection.wavelength,...
        detectionPupilFunctor,@(U,V) 1i*detectionPupilFunctor(U,V),... % circular polarization
        config.detection.objective.numericalAperture,config.sample.refractiveIndex);
    psf=detectionPsf.*repmat(lightSheetPsf,[imgSize(1) 1 1]);
    otf=fftn(ifftshift(psf));
    cutOffSpatialFrequency=2*config.detection.objective.numericalAperture./config.excitation.wavelength;
    otfSteps=1./pixelPitchInSample./imgSize;
    fRel=calcRadius(otfSteps,imgSize)./cutOffSpatialFrequency;
    NSR=ifftshift(fRel)./kSNR;
    clear fRel;
    filter=conj(otf)./(abs(otf).^2+NSR.^2); % Wiener filter
    
    % Deconvolve 3D
    restoredDataCube=fftn(recordedImageStack); % in place calculation to save memory
    restoredDataCube=ifftn(restoredDataCube.*filter,'symmetric');
    
    if nargout<1,
        close all;
        
        fig=figure('NumberTitle','off','Position',[100 100 800 600]);
        for zIdx=1:imgSize(3),
            set(fig,'Name',sprintf('frame idx %d',zIdx));
            axs(1)=subplot(1,2,1);
            imagesc(yRange*1e6,xRange*1e6,max(0,recordedImageStack(:,:,zIdx))); axis equal tight; title('pre-deconv');
            xlabel('x [\mum]'); ylabel('y [\mum]');
            colorbar();
            axs(2)=subplot(1,2,2);
            imagesc(yRange*1e6,xRange*1e6,max(0,restoredDataCube(:,:,zIdx))); title('deconvolved');
            xlabel('x [\mum]'); ylabel('y [\mum]'); axis equal tight;
            colorbar();

            linkaxes(axs);
            
            drawnow();
            
            pause(0.1);
        end
        
        clear deconvolvedImg;
    end
end
