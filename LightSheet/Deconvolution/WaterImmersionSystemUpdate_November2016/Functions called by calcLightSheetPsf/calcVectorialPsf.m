function [psf, psfField, varargout]=calcVectorialPsf(xRange,yRange,zRange,wavelength,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture,refractiveIndexOfSample,objectiveMagnification,objectiveTubeLength,projectionDimensions)
% Calculates the 3D point spread function at the grid specified by
% x/y/zRange for a pupil function given by pupilFunctorH/V(U,V) where U and V
% are normalized Cartesian pupil coordinates, in a medium with
% refractive index refractiveIndexOfSample and an objective. The horizontal
% polarization given by pupilFunctorH is along the first dimension, and the
% vertical given by pupilFunctorV is along the second dimension of the
% output.

    debug=false;
    
    %Manual definition of inputs
    if (nargin<1 || isempty(xRange))
        xRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<2 || isempty(yRange))
        yRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<3 || isempty(zRange))
        zRange=[-5000:500:5000]*1e-9;
    end
    if (nargin<4 || isempty(wavelength))
        wavelength=532e-9;
    end
    if (nargin<5 || isempty(pupilFunctorH))
        % Don't change the following, it is a sensible default!
        pupilFunctorH=@(normalU,normalV) 0.0; % Along the first dimension
    end
    if (nargin<6)
        pupilFunctorV=[]; %Along the second dimension
    end
    if (nargin<7 || isempty(objectiveNumericalAperture))
        objectiveNumericalAperture=1.0;
    end
    if (nargin<8 || isempty(refractiveIndexOfSample))
        refractiveIndexOfSample=1.0;
    end
    if (nargin<9 || isempty(objectiveMagnification))
        objectiveMagnification=1;
    end
    if (nargin<10 || isempty(objectiveTubeLength))
        objectiveTubeLength=200e-3;
    end
    if (nargin<11)
        projectionDimensions=[];
    end
    
    %When scalars are specified instead of function, assume constant input fields
    if (~isa(pupilFunctorH,'function_handle') && isscalar(pupilFunctorH))
        pupilFunctorH=@(normalU,normalV) pupilFunctorH;
    end
    if (~isa(pupilFunctorV,'function_handle') && isscalar(pupilFunctorV))
        pupilFunctorV=@(normalU,normalV) pupilFunctorV;
    end
    vectorialCalculation=~isempty(pupilFunctorV);
    if (~vectorialCalculation)
        logMessage('Starting a scalar calculation of the PSF...');
    end
    
    % Check how many multi-photon orders of the intensity have to be calculated
    highestOrderIntensityRequired=max(1,nargout-1);
    
    focalLengthInSample=objectiveTubeLength/objectiveMagnification; %TODO: Check for correctness
    
    objectiveSinMaxHalfAngleInSample=objectiveNumericalAperture/refractiveIndexOfSample;
    
    % Determine the requested step size for each non-singleton dimension
    sampleDelta=zeros(1,3);
    if (length(xRange)>1)
        sampleDelta(1)=diff(xRange(1:2));
    end
    if (length(yRange)>1)
        sampleDelta(2)=diff(yRange(1:2));
    end
    if (length(zRange)>1)
        sampleDelta(3)=diff(zRange(1:2));
    end
    
    minPupilSize=[1 1]*256; % To ensure that rapid pupil changes are properly sampled
    maxPupilSize=[1 1]*1024; % To avoid memory problems
    
    %The minimum pupil size to avoid PSF replication
    requiredPupilSize=2*ceil([length(xRange) length(yRange)].*sampleDelta(1:2)/(wavelength/(objectiveSinMaxHalfAngleInSample*refractiveIndexOfSample)));
    if (requiredPupilSize>maxPupilSize)
        logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle.\nThis will cause replicas.',[maxPupilSize,requiredPupilSize]);
    end
    
    %Check if pupil sampling not too sparse for the defocus we intend to simulate
    maxSampleNA=objectiveSinMaxHalfAngleInSample*(1-0.25/(max(maxPupilSize)/2)); % Use of 'max' because the sampling rate near the edge for NA=1 diverges
    minPupilSizeToHandleDefocus=minPupilSize+[1 1]*max(abs(zRange/(wavelength/refractiveIndexOfSample)))*4*objectiveSinMaxHalfAngleInSample*maxSampleNA/sqrt(1-maxSampleNA^2);
    if (minPupilSize<minPupilSizeToHandleDefocus)
        if debug,
            logMessage('A minimum pupil size of (%0.0f,%0.0f) is required to handle the specified defocus.',minPupilSizeToHandleDefocus);
        end
        requiredPupilSize=max(requiredPupilSize,minPupilSizeToHandleDefocus);
    end
    if (minPupilSizeToHandleDefocus>maxPupilSize)
        logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle. This can cause aliasing artefacts.',[maxPupilSize,minPupilSizeToHandleDefocus]);
    end
    
    pupilSize=ceil(min(maxPupilSize,max(minPupilSize,requiredPupilSize)));
    if debug,
        logMessage('Setting the pupil sampling grid size to (%0.0f,%0.0f)',pupilSize);
    end
    
    %Choose the pupil grid
    wavelengthInSample=wavelength/refractiveIndexOfSample;
    uRange=objectiveSinMaxHalfAngleInSample*2*[-floor(pupilSize(1)/2):floor((pupilSize(1)-1)/2)]/pupilSize(1);
	vRange=objectiveSinMaxHalfAngleInSample*2*[-floor(pupilSize(2)/2):floor((pupilSize(2)-1)/2)]/pupilSize(2);
    
    [U,V]=ndgrid(uRange,vRange);
    sinApAngle2=U.^2+V.^2;
    apertureFieldTransmission=double(sinApAngle2<objectiveSinMaxHalfAngleInSample^2);
    apertureArea=numel(U)*pi/4;
    sinApAngle=sqrt(apertureFieldTransmission.*sinApAngle2);
    cosApAngle=apertureFieldTransmission.*sqrt(1-apertureFieldTransmission.*sinApAngle2);
    clear sinApAngle2;
    
    %Scale so that the total intensity is 1 for a unity uniform
    %illumination
    apertureFieldTransmission=apertureFieldTransmission./sqrt(apertureArea);
    
    %Calculate the 2D pupil function
    pupilFunctionX=apertureFieldTransmission.*pupilFunctorH(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
    if (vectorialCalculation)
        pupilFunctionY=apertureFieldTransmission.*pupilFunctorV(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
        %Convert pupil function to polar coordinates
        T=atan2(V,U); CT=cos(T); ST=sin(T);
        pupilFunctionR =  CT.*pupilFunctionX+ST.*pupilFunctionY; % Radial component is rotated by the focusing
        pupilFunctionA = -ST.*pupilFunctionX+CT.*pupilFunctionY; % Azimutal component is unaffected by the focusing
        %Calculate the polarization change due to focussing
        pupilFunctionZ = sinApAngle.*pupilFunctionR;
        pupilFunctionR = cosApAngle.*pupilFunctionR;
        %Convert back to carthesian coordinates
        pupilFunctionX = CT.*pupilFunctionR-ST.*pupilFunctionA;
        pupilFunctionY = ST.*pupilFunctionR+CT.*pupilFunctionA;
        clear pupilFunctionR pupilFunctionA CT ST T apertureFieldTransmission sinApAngle cosApAngle;

        pupilFunction2D=cat(3,pupilFunctionX,pupilFunctionY,pupilFunctionZ);
        clear pupilFunctionX pupilFunctionY pupilFunctionZ;
    else
        pupilFunction2D=pupilFunctionX;
        clear pupilFunctionX;
    end
    
    %Calculate the focal plain fields using subroutine czt2andDefocus
    [psfField psfIntensities]=czt2andDefocus(pupilFunction2D,objectiveSinMaxHalfAngleInSample,xRange/wavelengthInSample,yRange/wavelengthInSample,zRange/wavelengthInSample, focalLengthInSample/wavelengthInSample, projectionDimensions, highestOrderIntensityRequired);

    % Rename the output
    psf=psfIntensities(:,:,:,1);
    if (size(psfIntensities,4)>2)
        varargout=mat2cell(psfIntensities(:,:,:,2:end),size(psfIntensities,1),size(psfIntensities,2),size(psfIntensities,3),ones(1,size(psfIntensities,4)-1));
    else
        if (size(psfIntensities,4)==2)
            varargout={psfIntensities(:,:,:,2)};
        end
    end
    clear psfIntensities;
    
    if (nargout==0)
        close all;
        %Display results
        maxNormalization=1./max(abs(psf(:)));
        for (zIdx=1:size(psf,3))
            subplot(2,2,1);
            showImage(psf(:,:,zIdx).'.*maxNormalization,[],xRange*1e6,yRange*1e6);
            title(sprintf('Total intensity for z=%0.3f \\mu m',zRange(zIdx)*1e6));
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,2);
            showImage(psfField(:,:,zIdx,3).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ez for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,3);
            showImage(psfField(:,:,zIdx,1).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ex for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,4);
            showImage(psfField(:,:,zIdx,2).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ey for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            drawnow();
            
            logMessage('Total intensity = %0.3f%%',100*sum(sum(psf(:,:,zIdx))));
            
            pause(1/30);
        end
           
        clear psf; % Don't litter on the command prompt
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, psfIntensities]=czt2andDefocus(x,objectiveSinMaxHalfAngleInSample,nxRange,nyRange,nzRange, nfocalLengthInSample, projectionDimensions, highestOrderIntensityRequired)
% Calculate the partial spectrum of x using the chirp z transform.
% This returns the complex field.
%
% x is the pupil and should not be ifftshifted, implicit zero padding to the right!
% objectiveSinMaxHalfAngleInSample: the input matrix must cover this disk exactly.
% the following arguments, also specifyable as a list, are:
%     nxRange and nyRange: the sample points in the units of wavelength in the sample medium.
%     nzRange: the sample points in the z dimension in units of wavelength in the sample.
%     nfocalLengthInSample: (optional, default infinite) The focal length specified in units of wavelength.
%     projectionDimensions: (optional, default none) The dimension along which an integration is done
%     highestOrderIntensityRequired: (optional, default depends on nargout) If specified, (higher order) intensities upto this number are returned as well.
%
% If more than one output argument is specified, the first and higher order
% intensities will be returned as 3D arrays stacked into a single 4D array.

    %Manual definition of inputs
    if (nargin<6 || isempty(nfocalLengthInSample))
        nfocalLengthInSample=Inf; %Assuming that focusLength >> z
    end
    if (nargin<7)
        projectionDimensions=[];
    end
    if (nargin<8)
        highestOrderIntensityRequired=max(0,nargout-1);
    end
    
    %Prepare the output matrix with zeros
    inputSize=size(x);
    outputSize=[length(nxRange) length(nyRange) length(nzRange), inputSize(3:end)]; % Add dimension
    if (~isempty(projectionDimensions))
        outputSize(projectionDimensions)=1;
    end
    f=zeros(outputSize,class(x));
    psfIntensities=zeros([outputSize(1:3) highestOrderIntensityRequired],class(x));
    
    uRange=objectiveSinMaxHalfAngleInSample*2*[-floor(inputSize(1)/2):floor((inputSize(1)-1)/2)]/inputSize(1);
	vRange=objectiveSinMaxHalfAngleInSample*2*[-floor(inputSize(2)/2):floor((inputSize(2)-1)/2)]/inputSize(2);
    [U,V]=ndgrid(uRange,vRange);
    R2=min(1.0,U.^2+V.^2);
    clear U V;
    cosHalfAngleInSampleMatrix=sqrt(1-R2);
    clear R2;
    
    %Loop through the z-stack, calculating the PSF at each slice
    for zIdx=1:length(nzRange)
        normalizedZ=nzRange(zIdx);
        phaseOffsetOfAxialWaveletInRad=2*pi*mod(normalizedZ,1); %Center on focal point
        % Calculate the phase delay due to z-displacement with respect to
        % the wavelet leaving the center of the pupil. This causes the Gouy
        % phase shift.
        if (~isinf(nfocalLengthInSample))
            % Geometrical difference between the axial and the off-axis ray
            relativeDefocusPhaseDelayAcrossPupilInRad=2*pi*(...
                sqrt(nfocalLengthInSample^2+2*nfocalLengthInSample*cosHalfAngleInSampleMatrix*normalizedZ+normalizedZ^2)...
                -(nfocalLengthInSample+normalizedZ));
        else
            unityDefocusInRad=2*pi*(cosHalfAngleInSampleMatrix-1); %Assuming that focalLength >> z
            % Approximate the above equation for nfocalLengthInSample >> normalizedZ
            relativeDefocusPhaseDelayAcrossPupilInRad=normalizedZ*unityDefocusInRad;
        end
        pupil=x.*repmat(exp(1i*(phaseOffsetOfAxialWaveletInRad+relativeDefocusPhaseDelayAcrossPupilInRad)),[1 1 inputSize(3:end)]);
        
        %Call function czt2fromRanges, which frames the pupil function (the
        %single-slice PSF) appropriately to to calculate the partial
        %spectrum using the built-in function czt.
        psfSlice=czt2fromRanges(pupil,nxRange*2*objectiveSinMaxHalfAngleInSample,nyRange*2*objectiveSinMaxHalfAngleInSample);
        
        % Project the output before continuing to save memory
        psfSliceIntensity=sum(abs(psfSlice).^2,3);
        psfSliceIntensity=repmat(psfSliceIntensity,[1 1 highestOrderIntensityRequired]);
        for (photonNb=1:highestOrderIntensityRequired)
            psfSliceIntensity(:,:,photonNb)=psfSliceIntensity(:,:,photonNb).^photonNb;
        end
        if (~isempty(projectionDimensions))
            for (projIdx=1:size(projectionDimensions,2))
                projectionDimension=projectionDimensions(projIdx);
                if (projectionDimension>=3)
                    projectionDimension=projectionDimension-1;
                end
                psfSlice=sum(psfSlice,projectionDimension);
                psfSliceIntensity=sum(psfSliceIntensity,projectionDimension);
            end
        end
        if (any(projectionDimensions==3))
            f(:,:,1,:)=f(:,:,1,:)+psfSlice;
            psfIntensities(:,:,1,:)=psfIntensities(:,:,1,:)+psfSliceIntensity;
        else
            f(:,:,zIdx,:)=psfSlice;
            psfIntensities(:,:,zIdx,:)=psfSliceIntensity;
        end
    end
end
