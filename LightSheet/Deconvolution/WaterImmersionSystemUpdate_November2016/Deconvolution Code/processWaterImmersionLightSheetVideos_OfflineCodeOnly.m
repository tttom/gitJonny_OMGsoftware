function processWaterImmersionLightSheetVideos_OfflineCodeOnly(folderNames,reprocessData,deleteAvi,centerOffset,scanShear,perspectiveScaling)
% Reads the recorded data, converts it to matrices and deconvolves.
%
% Input:
%      folderNames: a string or cell array of strings indicating folders with avi and json files.
%      reprocessData: When false, date which already has a mat file with
%                     a restoredDataCube matrix in it will be skipped.
%      centerOffset: The offset of the light sheet beam waist in meters given
%                    as a two-element vector containing the vertical (swipe direction, Y) and the
%                    horizontal (propagation direction, X) offset, respectivelly.
%      scanShear: The fraction change in x and y (vertical and horizontal) when moving in z (axially)
%      perspectiveScaling: The linear change in [x,y]=[vertical,horizontal] magnification when increasing z by one meter, units m^-1.

    %Manual input
    if (nargin<1 || isempty(folderNames))
%         folderNames={'E:\2015-12-08_JavierTello_BeadInjectedMouseBrain_Cleared'};
%         folderNames={'G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\__Acquired Data'};
%         folderNames = {'H:\Stored Files\NEW_LSM_SYSTEM_RESULTS\2016-12-07_CalibrationMeasurements'};
        folderNames = {'H:\Stored Files\NEW_LSM_SYSTEM_RESULTS\2016-12-09_pixelreassignment\Airystep200_scan100_int25'};

    end
    if (nargin<2) 
        reprocessData=true;
        %If true, the program skips .avi files that have the _PROCESSED
        %suffix (see below) and performs the deconvolution using the
        %outputted .mat file.  (This is good for, e.g., testing different
        %shear and scaling parameters, etc.)
    end
    if (nargin<3)
        deleteAvi=false;
        %If true, the raw data files are DELETED and replaced by a smaller
        %.avi file of the same name plus the suffix _PROCESSED.  Set to 
        %FALSE unless there is a copy of the raw data stored elsewhere.
    end
    if (nargin<4)
%         centerOffset=[0,1e-05];
        centerOffset=[0 0]*1e-6; % [y x] in metres
    end
    if (nargin<6)
        %perspectiveScaling: two-element vector in m^-1, components ordered vertical horizontal
        perspectiveScaling=[0 0];
%         perspectiveScaling=[0.5796 0.7339]*1e-3;
%         perspectiveScaling=[0.0008344,-0.0002917];
    end
    if (nargin<5)
        %scanShear also a two-element vector in m^-1 ([vertical
        %horizontal])
        scanShear=[0 0];
%         scanShear=[-0.0594,-0.0044];
    end
 
    if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    %Load the default configuration from the .json file
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'/waterImmersion.json');
    defaultConfig=loadjson(defaultConfigFileName);
    
    % Go through all the specified folders
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        logMessage('Checking folder %s for recorded videos.',folderName);
        experimentConfigFileName=fullfile(folderName,'experimentConfig.json');
        if (exist(experimentConfigFileName,'file'))
            experimentConfig=loadjson(experimentConfigFileName);
        else
            experimentConfig=struct();
            experimentConfig.detector=struct();
            experimentConfig.detector.center=[0 0];
            experimentConfig.detector.scanShear=[0 0];
            experimentConfig.detector.perspectiveScaling=[0 0];
        end
        if (~isempty(centerOffset))
            experimentConfig.detector.center=centerOffset;
        end
        if (~isempty(scanShear))
            experimentConfig.detector.scanShear=scanShear;
        end
        if (~isempty(perspectiveScaling))
            experimentConfig.detector.perspectiveScaling=perspectiveScaling;
        end
        savejson([],experimentConfig,experimentConfigFileName);
        experimentConfig=structUnion(defaultConfig,experimentConfig);
        
        % and process all videos
        videoFileNameList=dir(strcat(folderName,'/*.avi'));
        for (fileName={videoFileNameList.name})
            fileName=fileName{1}(1:end-4);
            
            %check if the file has "_PROCESSED" suffix (appended if
            %deleteAvi was set to 'true' during a previous analysis of the
            %same data sets).
            if ~isempty(strfind(fileName,'_PROCESSED'))
                fileName=fileName(1:end-10);
                skipAvi=true;
            else
                skipAvi=false;
            end
            
            filePathAndName=strcat(folderName,'/',fileName);
            logMessage('Processing %s...',filePathAndName);
            
            % Load specific configuration description
            configFile=strcat(filePathAndName,'.json');
            inputFileName=strcat(filePathAndName,'.avi');
            outputFileName=strcat(filePathAndName,'.mat');
            reprocessThisFile=true;
            if (~reprocessData && exist(outputFileName,'file'))
                storedVariables = whos('-file',outputFileName);
                if (ismember('restoredDataCube', {storedVariables.name}))
                    reprocessThisFile=false;
                    logMessage('Already done %s, skipping it!',outputFileName);
                end
            end
            if (reprocessThisFile)
                if (exist(configFile,'file'))
                    try
                        specificConfig=loadjson(configFile);
                        setupConfig=structUnion(experimentConfig,specificConfig);
                    catch Exc
                        logMessage('Could not read json config file %s, assuming defaults!',configFile);
                        setupConfig=experimentConfig;
                    end
                else
                    logMessage('Description file with extension .json is missing, assuming defaults!');
                end
                
                % Load recorded data
                if skipAvi
                    logMessage('Loading %s...',outputFileName);
                    try
                        load(outputFileName,'recordedImageStack');
                    catch Exc
                        logMessage('Failed to load data stack from %s!',outputFileName);
                        recordedImageStack=[];
                    end
                else  
                    logMessage('Loading %s...',inputFileName);
                    try
                        recordedImageStack=readDataCubeFromAviFile(inputFileName);
                    catch Exc
                        logMessage('Failed to load data stack from %s!',inputFileName);
                        recordedImageStack=[];
                    end
                end

                if (~isempty(recordedImageStack))
                    if skipAvi
                        % Save time by not overwriting with identical data
                        logMessage('Saving recorded data to %s...',outputFileName);
                        save(outputFileName,'setupConfig', '-append');
                    else
                        % Store partial results
                        logMessage('Saving recorded data to %s...',outputFileName);
                        save(outputFileName,'recordedImageStack','setupConfig', '-v7.3');
                    end

                    % Call the function deconvolveRecordedDataStack to
                    % perform the deconvolution
                    logMessage('Starting image reconstruction...');
                    [recordedImageStack lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,setupConfig);
                    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab

                    % Append the rest of the results
                    logMessage('Saving restored data cube to %s...',outputFileName);
                    save(outputFileName,'restoredDataCube','xRange','yRange','zRange','tRange','ZOtf','lightSheetPsf','lightSheetOtf','lightSheetDeconvFilter','-append');
                    clear restoredDataCube;
                    
                    %.avi file deletion (if deleteAvi==true)
                    if (deleteAvi && ~skipAvi)
                        % Delete original .avi file
                        delete(inputFileName);
                        % Replace with small placeholder file
                        placeHolderFileName=strcat(filePathAndName,'_PROCESSED.avi');
                        writerObj=VideoWriter(placeHolderFileName);
                        open(writerObj);
                        writeVideo(writerObj,zeros(2));
                        close(writerObj);
                    end
                    
                end
            else
                logMessage('Skipping file %s, already done.',inputFileName);
            end
        end
        
        % Checks for subfolders and handles these recursively
        directoryList=dir(folderName);
        for listIdx=1:length(directoryList),
            if directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.'
                expandedFolderName=strcat(folderName,'/',directoryList(listIdx).name);
                processWaterImmersionLightSheetVideos_OfflineCodeOnly(expandedFolderName,reprocessData,deleteAvi);
            end
        end
        
    end
end
