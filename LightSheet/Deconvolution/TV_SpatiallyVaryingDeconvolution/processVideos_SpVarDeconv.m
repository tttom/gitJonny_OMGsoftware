%
%
%
function processVideos_SpVarDeconv(folderNames,detNA,PMalpha)
%Tom Vettenburg's deconvolveCubicLightSheet code wrapped in the nested
%folder search and recursion from processWaterImmersionLightSheetVideos
%(the normal, spatially-invariant deconvolution code). It contains none of
%the other frills of the processWaterImmersion code, like the option to
%delete .avi files as you go, or to skip files that have already been
%processed. (So on that last point, take heed: if you run this code on
%already-processed movies, it WILL overwrite your .mat files. If you want
%multiple deconvolutions, move the .mat files to another folder first.)
    if nargin<1 || isempty(folderNames)
        folderNames={'E:\FatmaZohra-LightSheetData_2016-09-08\2016-09-08 18_47_55.289'};
        %experimentFilePath=fullfile(rootFolder,'recording0_lambda532nm_alpha-7_beta100'); % NA 0.42
    end
    
     if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        logMessage('Checking folder %s for recorded videos.',folderName);
        % and process all videos
        videoFileNameList=dir(strcat(folderName,'/*.avi'));
        for (fileName={videoFileNameList.name})
            fileName=fileName{1}(1:end-4);
            experimentFilePath=strcat(folderName,'/',fileName);
            logMessage('Processing %s...',experimentFilePath);

    %-------------------
    [experimentPath,fileName]=fileparts(experimentFilePath);
    dataFilePath=fullfile(experimentPath,[fileName,'.avi']);
    configFilePath=fullfile(experimentPath,[fileName,'.json']);
    outputFilePath=fullfile(experimentPath,[fileName,'.mat']);
  
    %Load the default configuration
    functionName=mfilename();
    defaultConfigPath=mfilename('fullpath');
    defaultConfigPath=defaultConfigPath(1:end-length(functionName));
    defaultConfig=loadjson([defaultConfigPath,'/waterImmersion.json']);
    experimentConfig=loadjson(configFilePath);
    config=structUnion(defaultConfig,experimentConfig);
    
    % Fix and extend the configuration with useful information
    % (could be written by LabVIEW so the json files it generates contain
    % all information needed)
    %config.excitation.objective.numericalAperture = 0.70;
    config.excitation.gaussianIlluminationStd = 0.50; %0.50, an estimate % not specified in json file
    config.modulation.alpha = -config.modulation.alpha; % for some reason this is inverted now
    config.detection.tubeLength = 190e-3; % Set-up changed from the expected 200e-3, otherwise the magnification doesn't match
    config.detection.objective.irisDiameter = detNA; % in meters
    backApertureDiameterOfDetection = ...
        2*tan(asin(config.detection.objective.numericalAperture/config.detection.objective.refractiveIndex))*...
        config.detection.objective.tubeLength/config.detection.objective.magnification;
    config.detection.objective.fractionOfNumericalApertureUsed = config.detection.objective.irisDiameter/backApertureDiameterOfDetection;
    
    config.detection.mask.diameter = 10e-3; % what you told me
    config.detection.mask.alpha = 20.65; % in wavelengths, over its diameter (i.e. not considering the iris here)
%     config.detection.mask.alpha = 0.0; % no mask in detection path
    config.detection.wavelength = 532e-9; % the approximate emission wavelength
%     config.detection.mask.beta=0; % use this for the classic cubic
    config.detection.mask.beta = -3*config.detection.mask.alpha;% generalized cubic
    config.detection.mask.rotationAngle = 0; % update this if you place the mask differently
    config.detection.mask.fractionOfDiameterUsed = config.detection.objective.irisDiameter / config.detection.mask.diameter;
    
    lightSheetOffset = [0 0 0]*1e-6; %[0 32e-6 0]*1e-6; % in case of no misalignment between the light sheet focus and the center of the FOV (x,y) and focus (z).
    shearFactors = [0.0002 0.0302]; %[0.0008 -0.0019]; %[0.0673 -0.0133]; % it appears that the stage translates at an angle slightly non-orthogonal. That is ok as long as it doesn't change
    viewSize=[128 1024 500]; % in pixels, a subset makes processing faster and doesn't overload the output image
%     viewSize=[128 256 128]; %[256, 128, 256]; % in pixels, a subset makes processing faster and doesn't overload the output image
    
    % The number of point-spread functions (PSFs) to use for spatially-variant deconvolution
    % affects. This should be chosen wisely as it also defines in how many sections the 
    % field-of-view is divided. For each section the PSF is calculated at the center of 
    % its corresponding section, the deconvolution kernel is calculated and the section is
    % deconvolved with the linear interpolation of all neighboring kernels so that it is smooth.
    % Pick the number of PSFs the sections size (viewSize./nbPsfs) is about twice as large 
    % as the PSF intensity spread. Too many PSFs will be  inaccurate due to potential cropping 
    % and overlap. Too few PSFs may fail to model some of the spatial variation.
    % Note also that a larger number of PSFs is more memory efficient because of smallar sub section sizes. 
    nbPsfs = [1 15 1]; % spatially variant, 15 different psfs along the light sheet propagation axis
%     nbPsfs = [1 1 1]; % spatially invariant, a single psf = also far more memory required because all work is done at once.

psfSize = floor(viewSize./nbPsfs); 
    
    % Extract the information we need for the deconvolution
    actualMagnification = config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    pixelPitchInSample = config.detector.pixelSize/actualMagnification;
    pixelPitchInSample(3)=norm(median(diff(config.stagePositions.target)));
    detectionObjectiveNumericalAperture = config.detection.objective.numericalAperture*config.detection.objective.fractionOfNumericalApertureUsed;
    kSNR = config.sample.signalLevel/config.sample.backgroundLevel; % could be estimated from data in principle
    detectionAlpha = config.detection.mask.alpha*config.detection.mask.fractionOfDiameterUsed^3; % this only works because the masks are all third order terms
    detectionBeta = config.detection.mask.beta*config.detection.mask.fractionOfDiameterUsed^3; % this only works because the masks are all third order terms
    % The following is the same as:
    %     pupilFunction=@(U,V) exp(2i*pi*(detectionAlpha*(U.^3+V.^3)+detectionBeta*(U.^2.*V+U.*V.^2)));
    % but rotated by config.detection.mask.rotationAngle
    cRA=cos(config.detection.mask.rotationAngle);
    sRA=sin(config.detection.mask.rotationAngle);
    %Generalized Cubic pupil function
    pupilFunctionGeneral = @(U,V) exp(2i*pi*(detectionAlpha*((U*cRA-V*sRA).^3+(U*sRA+V*cRA).^3)+detectionBeta*((U*cRA-V*sRA).^2.*(U*sRA+V*cRA)+(U*cRA-V*sRA).*(U*sRA+V*cRA).^2)));
    %2D Cubic pupil function
    pupilFunction2D = @(U,V) exp(2i*pi*(detectionAlpha*((U*cRA-V*sRA).^3+(U*sRA+V*cRA).^3)));
    pixelPitchInSample(3) = norm(median(diff(config.stagePositions.target)));
    if PMalpha == 20.65
        pupilFunction = pupilFunction2D;
    elseif PMalpha == 15.31
        pupilFunction = pupilFunctionGeneral;
    elseif PMalpha == 0
        pupilFunction = 1;
    end
    
    % Load the image data
    logMessage('Loading recorded data...');
    recordedImageStack = readDataCubeFromAviFile(dataFilePath,[],1);
    dataSize = [size(recordedImageStack) size(config.stagePositions.actual,1)];
    % define the sample grid
    xRange = pixelPitchInSample(1)*([1:dataSize(1)]-floor(dataSize(1)/2)-1);
    yRange = pixelPitchInSample(2)*([1:dataSize(2)]-floor(dataSize(2)/2)-1);
    zRange = pixelPitchInSample(3)*([1:dataSize(3)]-floor(dataSize(3)/2)-1);
    % center the view
    viewCenter = 1+floor(dataSize./2);
    
    if isempty(viewSize),
        viewSize = dataSize;
    else
        viewSize=min(viewSize,dataSize);
    end
    for dimIdx=1:3,
        rangeIndexes{dimIdx}=viewCenter(dimIdx)+[-floor(viewSize(dimIdx)/2):floor((viewSize(dimIdx)-1)/2)];
    end
    recordedImageStack = readDataCubeFromAviFile(dataFilePath,[],rangeIndexes);
    % Adapt the ranges
    dataSize = cellfun(@(r)numel(r),rangeIndexes);
    xRange = xRange(rangeIndexes{1});
    yRange = yRange(rangeIndexes{2});
    zRange = zRange(rangeIndexes{3});
    
    % determine the Point-spread Function
    [xRangePsfCalc, yRangePsfCalc, zRangePsfCalc] = calcRanges(psfSize, pixelPitchInSample);
    logMessage('Calculating the detection point-spread function...');
    detectionPsfs = calcVectorialPsf(single(xRangePsfCalc),single(yRangePsfCalc),single(zRangePsfCalc),config.detection.wavelength,...
        pupilFunction,@(U,V) 1i*pupilFunction(U,V),... % circular polarization
        detectionObjectiveNumericalAperture,config.sample.refractiveIndex);
    
    logMessage('Calculating light sheet...');
    lightSheetPsf = calcLightSheetPsf(single(xRange),single(yRange-lightSheetOffset(2)),single(zRange-lightSheetOffset(3)),0,config.excitation,config.modulation,config.sample.refractiveIndex);
    
    logMessage('Calculating the PSF at %d points...',prod(nbPsfs));
    [psfs, tiledPsfs] = calcCombinedPsf(lightSheetPsf,detectionPsfs,nbPsfs,dataSize,pixelPitchInSample,shearFactors);
    
    % Align all psfs so they can be interpolated more accurately
    cellSize = dataSize./nbPsfs; % floating point number
    for psfIdx=1:prod(nbPsfs),
        cellCenter = round(1+(ind2subs(nbPsfs,psfIdx)-0.5).*cellSize);
        psfs(:,:,:,psfIdx) = straightenAiryBend(psfs(:,:,:,psfIdx),...
            yRange(cellCenter(2))+yRangePsfCalc*0,zRange(cellCenter(3))+zRangePsfCalc,config); % *0: shift each psf as a block
    end
    
%     centerPsfIdx = sub2ind(nbPsfs,1+floor(nbPsfs(1)/2),1+floor(nbPsfs(2)/2),1+floor(nbPsfs(3)/2));
    
    logMessage('Calculating the deconvolution filter...');
    otfs = fft(fft(fft(ifftshift(ifftshift(ifftshift(psfs,1),2),3),[],1),[],2),[],3);
    %otfs = otfs./repmat(otfs(1,1,1,:,:,:));
    otfs = otfs./max(abs(otfs(:)));
    clear psfs;
    otfGridPitch = 1./pixelPitchInSample./psfSize;
    cutOffSpatialFrequency = 2*detectionObjectiveNumericalAperture./config.detection.wavelength;
    fRel = calcDistanceToCenterOfGrid(single(psfSize),single(otfGridPitch))./cutOffSpatialFrequency;
    PNSR = (ifftshift(fRel)./kSNR).^2;
    clear fRel;
    filters = conj(otfs)./(abs(otfs).^2+repmat(PNSR,[1 1 1 nbPsfs]));
    clear PNSR otfs;
    kernels = fftshift(fftshift(fftshift(ifft(ifft(ifft(filters,[],1),[],2),[],3),1),2),3);
%     squeeze(real(sum(sum(sum(kernels))))).'./squeeze(real(sum(sum(sum(psfs))))).'
    clear filters;
 
    restoreData = @spVarDeconvolution;

    % Deconvolve 2D
    logMessage('Restoring the data...');
    restoredDataCube = restoreData(recordedImageStack,kernels,yRange,zRange,config);
    logMessage('Deconvolving the PSFs too...');
    restoredTiledPsfs = restoreData(tiledPsfs,kernels,yRange,zRange,config);
    
    if nargout<1,
        logMessage('Displaying...');
        close all;
        
        figs(1) = displayMultiView(xRange,yRange,zRange,{'detection','lightsheet','system','deconvolved'},...
            detectionPsfs,lightSheetPsf,tiledPsfs,restoredTiledPsfs);
        
        figs(2) = displayMultiView(xRange,yRange,zRange,{'psf','recorded','restored'},...
            tiledPsfs,recordedImageStack,restoredDataCube);
        
        for figIdx = 1:numel(figs),
            outputFigurePath = strcat(outputFilePath(1:end-4),sprintf('_%dx%dx%d_%d.fig',[nbPsfs figIdx]));
            logMessage('Saving figure to %s...',outputFigurePath);
            saveas(figs(figIdx),outputFigurePath);
        end
        
        %clear restoredDataCube;
    end
     logMessage('Saving to %s...',outputFilePath);
     save(outputFilePath,'restoredDataCube','xRange','yRange','zRange','tiledPsfs','recordedImageStack');
     
     %---------------------------------
     end
        
        % Check for subfolders and handles these recursively
        directoryList=dir(folderName);
        for listIdx=1:length(directoryList),
            if directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.'
                expandedFolderName=strcat(folderName,'/',directoryList(listIdx).name);
                processVideos_SpVarDeconv(expandedFolderName,detNA,PMalpha);
            end
        end
        
    end
end

% Remove the distortion
function img = straightenAiryBend(img,yRange,zRange,config)
    %logMessage('Straightening data...');
    effectiveNA = config.excitation.objective.numericalAperture*config.excitation.fractionOfNumericalApertureUsed;
    effectiveNA = effectiveNA*min(1,config.excitation.gaussianIlluminationStd);
    effectiveAlpha = config.modulation.alpha*min(1,config.excitation.gaussianIlluminationStd).^3;
    if effectiveAlpha~=0, % a division by zero would invalidate the result
        sinTheta = effectiveNA/config.excitation.objective.refractiveIndex;
        cosTheta = sqrt(1-sinTheta^2);

        lambdaMedium = config.excitation.wavelength/config.sample.refractiveIndex;
        wRange = -yRange.*(1-cosTheta)/lambdaMedium; % in waves
        zShift = -wRange.^2*(lambdaMedium/(3*effectiveAlpha*sinTheta/cosTheta)); % metric
        zShift = 2*zShift; % TODO: not sure why this is off by a factor of two
        freqSlopeZ = zRange./((pixelPitchInSample(3).^2).*size(img,3));
        for yIdx=1:size(img,2),
            % sub-pixel shift plane by plane using the spatial frequency domain
            shiftFt = repmat(exp(2i*pi*zShift(yIdx)*reshape(freqSlopeZ,1,1,[])),[size(img,1) 1 1]);
            sliceFt = fft(img(:,yIdx,:),[],3);
            shiftedSlice = ifft(sliceFt.*ifftshift(shiftFt),[],3,'symmetric');
            img(:,yIdx,:) = shiftedSlice;
        end
    end
end

function restoredDataCube = spVarDeconvolution(recordedImageStack,kernels,yRange,zRange,config)
    restoredDataCube = spatiallyVariantConvolution(recordedImageStack,kernels);
    restoredDataCube = real(restoredDataCube);
    restoredDataCube = straightenAiryBend(restoredDataCube,yRange,zRange,config); % straighten result to undo the psf alignment
end
    
function fig = displayMultiView(xRange,yRange,zRange,dataNames,varargin)
%         colorMap=interpolatedColorMap(1024,[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 1 1],[[0:5]*.05 1]);
    colorMap=hot(1024);
    
    dataSize = [numel(xRange), numel(yRange), numel(zRange)];
    
%     function extendedProj = extend(proj,dim)
%         newSize = dataSize;
%         newSize(dim) = 1;
%         for dimIdx = 1:3,
%             paddedRanges{dimIdx} = 1+floor(size(proj,dimIdx)/2) + [1:newSize(dimIdx)]-1-floor(newSize(dimIdx)/2);
%             paddedRanges{dimIdx} = min(size(proj,dimIdx),max(1,paddedRanges{dimIdx}));
%         end
%         extendedProj = proj(paddedRanges{:});
%     end
    function paddedProj = pad(proj,projDim)
        oldSize = size(proj);
        oldSize((end+1):3) = 1;
        newSize = dataSize;
        newSize(projDim) = 1;
        % Initialize with zeros
        paddedProj = zeros(newSize);
        offset = floor(newSize/2)-floor(oldSize/2);
        for dimIdx = 1:3,
            projLimits = [1 oldSize(dimIdx)];
            projLimits = [max(1-offset(dimIdx),projLimits(1)) min(newSize(dimIdx)-offset(dimIdx),projLimits(2))];
            projRanges{dimIdx} = [projLimits(1):projLimits(2)];
            dataRanges{dimIdx} = offset(dimIdx)+projRanges{dimIdx};
        end
        % copy
        paddedProj(dataRanges{:}) = proj(projRanges{:});
    end

    function p = centerSlice(data,dim)
        ranges = {[1:size(data,1)],[1:size(data,2)],[1:size(data,3)]};
        ranges{dim} = ranges{dim}(1+floor(end/2));
        p = data(ranges{:});
        p = pad(p,dim);
    end
    function p = maxProj(data,dim)
        p = max(data,[],dim);
        p = pad(p,dim);
    end
    proj = @(data,dim) maxProj(data,dim);

    fig = figure('Position',[100 100 1024 768]);
    nbDataSets = numel(varargin);
    for dataIdx=1:nbDataSets,
        background = max(0,quantile(varargin{dataIdx}(:),0.01));
        normalization = max(varargin{dataIdx}(:))-background;
        
        axs(dataIdx,1) = subplot(3,nbDataSets,dataIdx);
        img = mapColor((proj(varargin{dataIdx},3)-background)./normalization,colorMap);
        if size(img,1) == 1,
            img = repmat(img,[numel(xRange) 1 1]);
        end
        showImage(img,[],yRange*1e6,xRange*1e6);
        xlabel('y [\mum]'); ylabel('x [\mum]'); axis equal tight;
        title(dataNames{dataIdx});
        
        axs(dataIdx,2) = subplot(3,nbDataSets,nbDataSets+dataIdx);
        showImage(mapColor(squeeze(proj(varargin{dataIdx},1)-background).'./normalization,colorMap),[],yRange*1e6,zRange*1e6);
        xlabel('y [\mum]'); ylabel('z [\mum]'); axis equal tight;
        
        axs(dataIdx,3) = subplot(3,nbDataSets,2*nbDataSets+dataIdx);
        img = mapColor(squeeze(proj(varargin{dataIdx},2)-background).'./normalization,colorMap);
        if size(img,1) == 1,
            img = repmat(permute(img,[2 1 3]),[1 numel(xRange) 1]);
        end
        showImage(img,[],xRange*1e6,zRange*1e6);
        xlabel('x [\mum]'); ylabel('z [\mum]'); axis equal tight;
    end

    for viewIdx=1:size(axs,2),
        linkaxes(axs(:,viewIdx));
    end
    
    drawnow();
end

%
% Multiply the light sheet with the detection psfs
%
function [psfs, tiledPsfs] = calcCombinedPsf(lightSheetPsf,detectionPsfs,nbPsfs,dataSize,pixelPitchInSample,shearFactors)
    if nargin<3,
        nbPsfs = ones(1,ndims(detectionPsfs));
    end
    if nargin<6 || isempty(shearFactors),
        shearFactors = [0 0];
    end
    
    psfSize = size(detectionPsfs);
    [xRangePsfCalc, yRangePsfCalc, zRangePsfCalc] = calcRanges(psfSize, pixelPitchInSample);
    cellSize = dataSize./nbPsfs; % floating point number
    if any(cellSize<psfSize),
        logMessage('The PSF size is larger than the deconvolution cells!');
    end
    for dimIdx = 1:3,
        cellPosRanges{dimIdx} = 1+cellSize(dimIdx)*([1:nbPsfs(dimIdx)]-0.5);
    end
    psfs = repmat(detectionPsfs,[1 1 1 nbPsfs]); % safe space by in-place multiplication
    if nargout>=2,
        tiledPsfs = zeros(dataSize);
    end
    % Handle one local PSF and its cell at a time
    for psfIdx = 1:prod(nbPsfs),
        cellIndexes=ind2subs(nbPsfs,psfIdx);
        for dimIdx = 1:3,
            cellPos(dimIdx) = cellPosRanges{dimIdx}(cellIndexes(dimIdx)); % floating point number
        end
        psfOffset = round(cellPos) -1-floor(psfSize./2); % offset of psf with respect to the output data cube
        for dimIdx = 1:3,
            psfLimits = [1 psfSize(dimIdx)];
            psfLimits = [max(1-psfOffset(dimIdx),psfLimits(1)) min(dataSize(dimIdx)-psfOffset(dimIdx),psfLimits(2))];
            psfRangesInSubSection{dimIdx} = [psfLimits(1):psfLimits(2)];
            psfRanges{dimIdx} = psfOffset(dimIdx)+[psfLimits(1):psfLimits(2)];
        end
        % ranges = arrayfun(@(s) [1:s]-1-floor(s/2), psfSize, 'UniformOutput',false);
        lightSheetSubSection = zeros([1 psfSize(2:end)]);
        psfRangeY = psfRanges{2}(1+floor(end/2))*ones(1,numel(psfRanges{2}));
        lightSheetSubSection(1,psfRangesInSubSection{2:3}) = lightSheetPsf(1,psfRangeY,psfRanges{3});

        % Do the actual multiplication (in a memory efficient way)
        for xIdx=1:psfSize(1),
            psfs(xIdx,:,:,psfIdx) = psfs(xIdx,:,:,psfIdx).*lightSheetSubSection;
        end
        
        % Shear the psfs if requested
        if any(shearFactors~=0),
            %logMessage('Shearing PSF (%0.2f%%,%0.2f%%) in (x/z,y/z).',shearFactors*100);
            [shiftX,shiftY] = ndgrid(xRangePsfCalc*shearFactors(1)./((pixelPitchInSample(1).^2).*psfSize(1)),...
                            yRangePsfCalc*shearFactors(2)./((pixelPitchInSample(2).^2).*psfSize(2)));
            minus2iPiXpY = -2i*pi*(shiftX+shiftY);
            clear shiftX shiftY;
            for zIdx=1:psfSize(3),
                % sub-pixel shift plane by plane using the spatial frequency domain
                shiftFt = exp(minus2iPiXpY*zRangePsfCalc(zIdx));
                slice = psfs(:,:,zIdx,psfIdx);
                slice = ifft2(fft2(slice).*ifftshift(shiftFt),'symmetric');
                psfs(:,:,zIdx,psfIdx) = slice;
            end
        end
        if nargout>=2,
            % Combining into a single psf image for backward compatibility
            tiledPsfs(psfRanges{:}) = tiledPsfs(psfRanges{:})+psfs(psfRangesInSubSection{:},psfIdx);
        end
    end
end
