%%% Function:           wavefrontAberrationMeasurement_CizmarMethod2010.m
%%% Author:             Jonathan Nylk (University of St Andrews)
%%% Created:            20/10/2016
%%% Description:        This function performs a wavefront aberration
%%%                     measurement as described in [UPDATED PAPER
%%%                     REFERENCE]. Naturally this function is heavily
%%%                     influenced by the above stated manuscript.
%%%                     A number of cameras and SLM control options may be
%%%                     implemented.
%%%
%%%
%%% Inputs:             Camera_Type:    A string indicating the camera that
%%%                                     will be used.
%%%                     Camera_Settings:A structure containing relevant
%%%                                     parameters for setting the camera.
%%%                     SLM_Type:       A string indicating the SLM that
%%%                                     will be used.
%%%                     SLM_Settings:   A structure containing relevant
%%%                                     parameters for settign the SLM.
%%%                     Probe_Mode_Size:A 1x2 array detailing the desired
%%%                                     size of the modes to be used for
%%%                                     the measurement.
%%%                     
%%%
%%% Updates (latest first):
%%%         
%%%
%%% END


function wavefrontAberrationMeasurement_CizmarMethod2010(Camera_Type,Camera_Settings,SLM_Type,SLM_Settings,Probe_Mode_Size,Phase_Periods,Phase_Samples,Correction_HalfROI)

    if nargin < 1
        % Hamamatsu Orca Flash 4.0 v2
        Camera_Type = 'Hamamatsu_Orca_Flash4.0_v2_C11440';
    end
    
    if nargin < 2
       Camera_Settings = [];
       % Video input settings
       Camera_Settings.videoinput.adaptorname = 'hamamatsu'; % Adaptor name
       Camera_Settings.videoinput.deviceID = 1; % Device ID
       Camera_Settings.videoinput.format = 'MONO16_BIN4x4_512x512_FastMode'; % Format
       Camera_Settings.pixelDepth = 16; % Pixel bit depth
       Camera_Settings.pixelFormat = 'uint16'; % Pixel bit format
       % Triggering
       Camera_Settings.Triggering.triggerConfig = 'manual'; % Trigger Config
       Camera_Settings.Triggering.framesPerTrigger = 1; % Frames per trigger
       Camera_Settings.Triggering.triggerRepeat = Inf; % Trigger repeat
       % Exposure
       Camera_Settings.exposureTime = 0.02; %Exposure time (in seconds - I think)
    end
    
    if nargin < 3
        % Hamamatsu LCOS
        SLM_Type = 'Hamamatsu_LCOS';
    end
    
    if nargin < 4
        SLM_Settings = [];
        screenSize = get(0,'MonitorPositions');
        SLM_Settings.monitorNumber = 1; % SLM monitor number
        SLM_Settings.size.slmHeight = screenSize(1 + SLM_Settings.monitorNumber,4); % SLM height (pixels)
        SLM_Settings.size.slmWidth = screenSize(1 + SLM_Settings.monitorNumber,3); % SLM width (pixels)
        SLM_Settings.location.verticalStart = screenSize(1 + SLM_Settings.monitorNumber,2); % SLM vertical offset (pixels)
        SLM_Settings.location.horizontalStart = screenSize(1 + SLM_Settings.monitorNumber,1); % SLM horizontal offset (pixels)
        SLM_Settings.pixelDepth = 8; % Pixel bit depth
        SLM_Settings.pixelFormat = 'uint8'; % Pixel bit format
        SLM_Settings.modulation.horizontalGratingPeriod = 0.12; % pixels
        SLM_Settings.modulation.verticalGratingPeriod = 0; % pixels
        SLM_Settings.modulation.Range2Pi = 1; % Fraction of pixel bit range that equals 2pi
    end
    
    if nargin < 5
        % Must divide 100 with no remainder
        Probe_Mode_Size = 50; % Pixels
%         Probe_Mode_Size = 20; % Pixels
    end
    
    if nargin < 6
        Phase_Periods = 2;
    end
    
    if nargin < 7
        Phase_Samples = 16; % Samples per period
    end
    
    if nargin < 8
        Correction_HalfROI = 1; % Half-width of box to collect signal from (pixels)
    end
    
    Flag_Error = 0; % 0 if no (predicted) errors, 1 if (predicted) error.
    
    
    %%% Initialise camera
    switch Camera_Type
        case 'Hamamatsu_Orca_Flash4.0_v2_C11440'
            disp('Initialising camera: Hamamatsu Orca Flash 4.0 v2.')
            vidObj = videoinput(Camera_Settings.videoinput.adaptorname,Camera_Settings.videoinput.deviceID,Camera_Settings.videoinput.format);
            srcObj = getselectedsource(vidObj);
            triggerconfig(vidObj,Camera_Settings.Triggering.triggerConfig);
            vidObj.FramesPerTrigger = Camera_Settings.Triggering.framesPerTrigger;
            vidObj.TriggerRepeat = Camera_Settings.Triggering.triggerRepeat;
            srcObj.ExposureTime = Camera_Settings.exposureTime; %in seconds (I think)
            imSize = vidObj.VideoResolution;
            imHeight = imSize(1);
            imWidth = imSize(2);
            start(vidObj);
            disp('WARNING: Camera object not started for debugging purposes.')
            pause(0.5);
            disp('Camera initialisation complete')
        otherwise
            disp('No recognised camera type input. Exiting function.')
            Flag_Error = 1;
    end
    
    %%% Initialise SLM
    switch SLM_Type
        case 'Hamamatsu_LCOS'
            disp('Initialising SLM: Hamamatsu LCOS.')
            SLMWindow = figure('Name','SLMWindow','Position',[SLM_Settings.location.horizontalStart SLM_Settings.location.verticalStart SLM_Settings.size.slmWidth SLM_Settings.size.slmHeight]);
            SLMObj = image(ones(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'uint8') .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi);
            SLMAxes = get(SLMObj,'Parent');
            set(SLMAxes,'Position',[0 0 1 1]);
            colormap(gray(2.^SLM_Settings.pixelDepth));
            drawnow;shg;
            % Define SLM Coords
            [xSLM,ySLM] = meshgrid([1:SLM_Settings.size.slmWidth] - floor(SLM_Settings.size.slmWidth/2),[1:SLM_Settings.size.slmHeight] - floor(SLM_Settings.size.slmHeight/2));
            xSLM = single(xSLM);
            ySLM = single(ySLM);
            % Default Blazed Grating
            SLM_BlazedGrating = mod(angle(exp(2 * pi * 1i * SLM_Settings.modulation.horizontalGratingPeriod .* xSLM)...
                .* exp(2 * pi * 1i * SLM_Settings.modulation.verticalGratingPeriod .* ySLM)) / 2 / pi,1)...
                 .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi;
            SLMObj = image(SLM_BlazedGrating);
        otherwise
            disp('No recognised SLM type input. Exiting function.')
            Flag_Error = 1;
    end
  
    if Flag_Error == 0
        %%% Select pixel of interest
        disp('Use previewFrame Window to select point of image to perfrom correction on.')
%         previewFrame = zeros(imHeight,imWidth,Camera_Settings.pixelFormat); % Pre-allocation of memory not crucial here.
        trigger(vidObj);
        previewFrame = getdata(vidObj,1);
        previewFig = figure();
        imagesc(previewFrame);
        axis image;
        title('Click on point in image to perform correction on to continue');
        drawnow;shg;
        [camXCoord,camYCoord] = ginput(1);
        close(previewFig);
        camXCoord = round(camXCoord); % Closest pixel
        camYCoord = round(camYCoord); % Closest pixel
        
        %%% Set probe parameters
        noVerticalModes = SLM_Settings.size.slmHeight / Probe_Mode_Size;
        noHorizontalModes = SLM_Settings.size.slmWidth / Probe_Mode_Size;
        % Check integer number of modes in each dimension
        if mod(noVerticalModes,1) ~= 0 || mod(noHorizontalModes,1) ~= 0
            Probe_Mode_Size = 50;
            noVerticalModes = SLM_Settings.size.slmHeight / Probe_Mode_Size;
            noHorizontalModes = SLM_Settings.size.slmWidth / Probe_Mode_Size;
            disp('Probe_Mode_Size does not fit in integer tiling onto SLM, setting default Probe_Mode_Size.');
        end
        
        
    else
        disp('Error encountered during initialisation. Wavefront measurement not made.')
    end
    
    
    keyboard


    %%% De-initialise camera
    switch Camera_Type
        case 'Hamamatsu_Orca_Flash4.0_v2_C11440'
            delete(vidObj);
            clear vidObj srcObj;
        otherwise
    end
    
    %%% De-initialise SLM
    switch SLM_Type
        case 'Hamamatsu_LCOS'
            close SLMWindow;
            clear SLMWindow SLMAxes SLMObj;
        otherwise
    end

end