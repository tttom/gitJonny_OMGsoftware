function psf=calcLightSheetPsf(xRange,yRange,zRange,tilt,excitation,modulation,refractiveIndexOfSample)
% Calculates the 2D intensity profile created by a swiped light-sheet.
%
% Inputs:
%     tilt: a scalar indicating the tilt of the pupil, or thus the lateral
%           position of the light sheet (defined in
%           deconvolveRecordedImageStack).
%     excitation: struct with fields wavelength and objective, the latter
%                 must contain the fields magnification, tubeLength, and optionally fractionOfNumericalApertureUsed
%     modulation: a struct containing the fields alpha and beta or more
%     refractiveIndexOfSample: the refractive index of the sample medium
%Excitation, modulation, and refractive index are stored in the .json file
%associated with the raw data and are imported in processWaterImmersion....
%
% Returns:
%     psf = the single photon point-spread function (PSF)

%Manual definition of input values
    if (nargin<5 || isempty(excitation))
        excitation=struct('wavelength',532e-9,...
            'objective',struct('numericalAperture',0.80,'magnification',40,'tubeLength',200e-3),...
            'fractionOfNumericalApertureUsed',1.0);
    end
    
    if (isfield(excitation,'fractionOfNumericalApertureUsed'))
        numericalAperture=excitation.objective.numericalAperture*excitation.fractionOfNumericalApertureUsed;
    end
    if (nargin<7 || isempty(refractiveIndexOfSample))
        refractiveIndexOfSample=1.0;
    end
    if (nargin<1 || isempty(xRange))
        xRange=(7.4*1e-6/excitation.objective.magnification)*[-200:199];
        xRange=single(xRange);
    end
    if (nargin<2 || isempty(yRange))
        yRange=(7.4*1e-6/excitation.objective.magnification)*[-200:199];
        yRange=single(yRange);
    end
    if (nargin<3 || isempty(yRange))
        stageTranslationStepSize=0.1*1e-6;
        zRange=(stageTranslationStepSize*refractiveIndexOfSample)*[-250:249]; %Translation range (along z-axis)
        zRange=single(zRange);
    end
    if (nargin<4 || isempty(tilt))
        tilt=0;
    end
    if (nargin<6 || isempty(modulation))
        modulation=struct('alpha',0,'beta',1);
    end
    if (~isfield(excitation,'illuminationClippingFactors'))
        illuminationClippingFactors=0*[1 1; 1 1];
    else
        illuminationClippingFactors=excitation.illuminationClippingFactors;
    end
    if (~isfield(excitation,'gaussianIlluminationStd'))
        gaussianIlluminationStd=[];
    else
        gaussianIlluminationStd=excitation.gaussianIlluminationStd;
    end
    if (~isfield(excitation,'beamAngle'))
        beamAngle=[];
    else
        beamAngle=excitation.beamAngle;
    end
    
    % Check for attenuation compensation parameters and modify pupil
    % function if needed
    if isfield(modulation,'sigmaU')
        sigmaU = modulation.sigmaU;
        sigmaV = modulation.sigmaV;
    else
        sigmaU = 0;
        sigmaV = 0;
    end
    
    % Check for super-Gaussian window parameters and modify pupil
    % function if needed
    if isfield(modulation,'A_SG')
        A_SG = modulation.A_SG;
        k_SG = modulation.k_SG;
        sigma_SG = modulation.sigma_SG;
    else
        A_SG = 0;
        k_SG = 0;
        sigma_SG = 0;
    end
    
    %Define cubic phase mask for Airy beam PSF
    if ~isa(modulation,'function_handle')
        if ~isempty(beamAngle)
            UR=cos(beamAngle)*U-sin(beamAngle)*V;
            V=sin(beamAngle)*U+cos(beamAngle)*V;
            U=UR;
            clear UR;
        end
        modulation=@(U,V) (sqrt(U.^2+V.^2)>=(1-modulation.beta)).*exp(2i*pi*modulation.alpha*(U.^3+V.^3));
    end    
    
    projectionDimension=1;
    if (length(xRange)<=1)
        projRangeLength=512; %Nyquist sampling
        xRangeForProj=([1:projRangeLength]-floor(projRangeLength/2)-1)*0.5*0.5*excitation.wavelength/numericalAperture;
    else
        xRangeForProj=xRange;
    end
   
    crop=@(U,V) 1.0*(U>=-(1-illuminationClippingFactors(1,1))&U<=(1-illuminationClippingFactors(1,2)) & V>=-(1-illuminationClippingFactors(2,1))&V<=(1-illuminationClippingFactors(2,2)));
    if (~isempty(gaussianIlluminationStd))
        apodization=@(U,V) exp(-(U.^2+V.^2)./(2*gaussianIlluminationStd^2));
    else
        apodization=@(U,V) 1;
    end

    %Rotate the coordinate system so that X and Z are interchanged.
    pupilFunctor=@(U,V) crop(U,V).*apodization(U,V).*exp(2i*pi*tilt*V).*modulation(U,V);
    
    % Check for attenuation compensation parameters and modify pupil
    % function if needed
    if sigmaU ~= 0 || sigmaV ~= 0
        pupilFunctor = @(U,V) pupilFunctor(U,V)...
            .* exp(sigmaU .* (U - 1)) .* exp(sigmaV .* (V-1));
    end
    
    % Check for super-Gaussian window parameters and modify pupil
    % function if needed
    if A_SG ~= 0 && k_SG ~= 0 && sigma_SG ~=0
%         pupilFunctor = @(U,V) pupilFunctor(U,V)...
%             .* A_SG .* exp(-1 .* (k_SG .* (U + V)).^sigma_SG);
        pupilFunctor = @(U,V) pupilFunctor(U,V)...
            .* A_SG .* exp(-1 .* ((k_SG .* U).^sigma_SG + (k_SG .* V).^sigma_SG));
    end
    
    %Assume circular polarization propagating along the x-axis and
    %call the function calcVectorialPsf to calculate the PSF
    cache=Cache(); % Get the default cache
    key={xRangeForProj,zRange,yRange,{excitation.wavelength,1/sqrt(2),1i/sqrt(2),illuminationClippingFactors,gaussianIlluminationStd,beamAngle,pupilFunctor,tilt},numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension};
    if isfield(cache,key)
        logMessage('Loading PSF projection from cache...');
        %If available, loads a previously-calculated PSF with identical parameters and
        %uses it to save memory.
        psf=cache(key);
    else
        stopWatch=tic();
        psf=calcVectorialPsf(xRangeForProj,zRange,yRange,excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension);
        if toc(stopWatch)>10, % Only cache slow calculations
            cache(key)=psf;
        end
    end
        
    %Return to original coordinate system and define output
    psf=permute(psf,[1 3 2]);
    
    %display PSF if there are no outputs (i.e. if the function was not
    %called as part of a deconvolution; if it was, the PSF will be
    %displayed at the end of the calling function,
    %deconvolveRecordedImageStack).
    if (nargout==0 && numel(psf)>100)
        figure;
        tmp=squeeze(psf(:,1,:)).';
        tmp=tmp./repmat(max(tmp),[size(tmp,1) 1]);
        imagesc(xRange*1e6,zRange*1e6,tmp);axis equal
        clear psf;
    end
end