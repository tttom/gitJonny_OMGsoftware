%%% Function:           wavefrontCorrection_CizmarMethod2010.m
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
%%%                     Phase_Periods:  Integer indicating the number of
%%%                                     2Pi periods to sample.
%%%                     Phase_Samples:  Integer indicating the number of
%%%                                     samples to take per 2Pi period.
%%%                     Correction_HalfROI:
%%%                                     Integer indicating the half-width
%%%                                     of ROI to be used for correction.
%%%                                     Area of ROI is:
%%%                                       [centreV - Correction_HalfROI:
%%%                                       centreV + Correction_HalfROI,
%%%                                       centreH - Correction_HalfROI:
%%%                                       centreH + Correction_HalfROI].
%%%                     Delay_Timer:    Sets a delay timer which pauses the
%%%                                     program to allow the user to leave
%%%                                     the room before the correction
%%%                                     process begins.
%%%                     
%%%
%%% Updates (latest first):
%%%         
%%%
%%% END


function [SLM_Amplitude_Map,SLM_Phase_Map] = wavefrontCorrection_CizmarMethod2010(Camera_Type,Camera_Settings,SLM_Type,SLM_Settings,Probe_Mode_Size,Phase_Periods,Phase_Samples,Correction_HalfROI,Delay_Timer)

    if nargin < 1
        % Hamamatsu Orca Flash 4.0 v2
        Camera_Type = 'Hamamatsu_Orca_Flash4.0_v2_C11440';
%         Camera_Type = 'Test_No_Camera';
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
        Probe_Mode_Size = [50,50]; % [vert,horiz] Pixels
%         Probe_Mode_Size = [20,20]; % [vert,horiz] Pixels
    end
    
    if nargin < 6
        Phase_Periods = 2;
    end
    
    if nargin < 7
        Phase_Samples = 8; % Samples per period
    end
    
    if nargin < 8
        Correction_HalfROI = 1; % Half-width of box to collect signal from (pixels)
    end
    
    if nargin < 9
%         Delay_Timer = 60; % seconds (time to wait before starting)
        Delay_Timer = 0; % seconds
    end
    
    Flag_Error = 0; % 0 if no (predicted) errors, 1 if (predicted) error.
    Flag_No_Camera = 0; % 0 if camera set, 1 if no camera set.
    
    live_Cam = 1;
    
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
            pause(0.5);
            disp('Camera initialisation complete.')
        case 'Test_No_Camera'
            Flag_No_Camera = 1;
            disp('No camera set. Running without camera acquisition.')
        otherwise
            disp('No recognised camera type input. Exiting function.')
            Flag_Error = 1;
    end
    
    %%% Initialise SLM
    switch SLM_Type
        case 'Hamamatsu_LCOS'
            disp('Initialising SLM: Hamamatsu LCOS.')
            SLMWindow = figure('Name','SLMWindow','Position',[SLM_Settings.location.horizontalStart SLM_Settings.location.verticalStart SLM_Settings.size.slmWidth SLM_Settings.size.slmHeight]);
            SLMObj = image(ones(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'uint8')...
                .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi);
            SLMAxes = get(SLMObj,'Parent');
            set(SLMAxes,'Position',[0 0 1 1],'XTick',[],'YTick',[],'Box','Off');
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
            drawnow;shg;
        otherwise
            disp('No recognised SLM type input. Exiting function.')
            Flag_Error = 1;
    end
    
    %%% Other initialisation
    if Flag_No_Camera == 0
        previewFrame = zeros(imHeight,imWidth,Camera_Settings.pixelFormat);
        image_Frame = zeros(imHeight,imWidth,Camera_Settings.pixelFormat);
    end
    %%% Set probe parameters
    noVerticalModes = SLM_Settings.size.slmHeight / Probe_Mode_Size(1);
    noHorizontalModes = SLM_Settings.size.slmWidth / Probe_Mode_Size(2);
    % Check integer number of modes in each dimension
    if mod(noVerticalModes,1) ~= 0 || mod(noHorizontalModes,1) ~= 0
        Probe_Mode_Size = [50,50];
        noVerticalModes = SLM_Settings.size.slmHeight / Probe_Mode_Size(1);
        noHorizontalModes = SLM_Settings.size.slmWidth / Probe_Mode_Size(2);
        disp('Probe_Mode_Size does not fit in integer tiling onto SLM, setting default Probe_Mode_Size.');
    end
    ref_Mode_Vertical = round(noVerticalModes/2);
    ref_Mode_Horizontal = round(noHorizontalModes/2);
    ref_Amplitude_Mask = zeros(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'single');
    ref_Amplitude_Mask(Probe_Mode_Size(1) * (ref_Mode_Vertical - 1) + 1:Probe_Mode_Size(1) * ref_Mode_Vertical...
        ,Probe_Mode_Size(2) * (ref_Mode_Horizontal - 1) + 1:Probe_Mode_Size(2) * ref_Mode_Horizontal)...
        = 1;
    ref_Mask = ref_Amplitude_Mask .* SLM_BlazedGrating;
    % Set fitting parameters
    intensity_Matrix = zeros(noVerticalModes,noHorizontalModes,Phase_Periods * Phase_Samples,'single');
    cosine_Fitting_Matrix = zeros(noVerticalModes,noHorizontalModes,Phase_Periods * Phase_Samples,'single');
    cosine_Function = @(cosine_Params,t) cosine_Params(1) .* cos(cosine_Params(2) .* t + cosine_Params(3)) + cosine_Params(4);
    cosine_Params0 = [0,0,0,0];
    cosine_Parameter_Minimisation_Function = @(data,cosine_Params,t) sum((squeeze(data) - squeeze(cosine_Function(cosine_Params,t)).').^2);
    cosine_Fitting_Parameters = zeros(noVerticalModes,noHorizontalModes,4,'single'); %y = A * cos(B * x + C) + D (4 parameters to fit).
  
    if Flag_Error == 0
        %%% Preview
        if Flag_No_Camera == 0
            stop(vidObj);
            preview(vidObj);
            disp('check beam is present then continue.')
            keyboard
            stoppreview(vidObj);
            closepreview(vidObj);
            start(vidObj);
        else
            disp('No camera set. Skipping stage: Preview image with reference mode only and set intensity.')
        end
        
        %%% Select pixel of interest
        if Flag_No_Camera ==0
            disp('Use previewFrame Window to select point of image to perfrom correction on.')
            trigger(vidObj);
            previewFrame = getdata(vidObj,1);
            previewFig = figure();
            imagesc(round(imHeight / 4):round(imHeight * 3 / 4),round(imWidth / 4):round(imWidth * 3 / 4)...
                ,previewFrame(round(imHeight / 4):round(imHeight * 3 / 4),round(imWidth / 4):round(imWidth * 3 / 4)));
            axis image;
            title('Click on point in image to perform correction on to continue');
            drawnow;shg;
            [camXCoord,camYCoord] = ginput(1);
            close(previewFig);
            camXCoord = round(camXCoord); % Closest pixel
            camYCoord = round(camYCoord); % Closest pixel
        else
            disp('No camera set. Skipping stage: Select pixel of interest.')
        end
        
        %%% Preview image with reference mode only and set intensity
        if Flag_No_Camera == 0
            figure(SLMWindow);
            axes(SLMAxes);
            SLMObj = image(ref_Mask);
            drawnow;shg;
            pause(1/60); % Short delay to allow SLM to update
            stop(vidObj);
            preview(vidObj);
            disp('Set laser power then continue.')
            keyboard
            stoppreview(vidObj);
            closepreview(vidObj);
            start(vidObj);
        else
            disp('No camera set. Skipping stage: Preview image with reference mode only and set intensity.')
        end
        
        %%% Optional delay before starting correction (lets me leave the room)
        disp(strcat('Waiting for (',num2str(Delay_Timer),') seconds before starting wavefront measurement.'))
        pause(Delay_Timer); % delay in seconds
        
        %%% Wavefront aberration measurement loop
        disp('Beginning wavefront measurement procedure.')
        measurement_Timer = tic;
        % Loop over all vertical mode coords
        for V_idx = 1:noVerticalModes
            % Loop over all horizontal mode coords
            for H_idx = 1:noHorizontalModes
                if V_idx == ref_Mode_Vertical && H_idx == ref_Mode_Horizontal
                    % No need to probe the reference mode itself
                    intensity_Matrix(V_idx,H_idx,:) = 0;
                    intensity_Matrix(V_idx,H_idx,Phase_Periods * Phase_Samples) = 1; % 0 phase difference
                else
                    % Loop over all phase samples and periods
                    for P_idx = 1: Phase_Samples * Phase_Periods
                        % Set probe mask
                        probe_Amplitude_Mask = zeros(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'single');
                        probe_Amplitude_Mask(Probe_Mode_Size(1) * (V_idx - 1) + 1:Probe_Mode_Size(1) * V_idx...
                            ,Probe_Mode_Size(2) * (H_idx - 1) + 1:Probe_Mode_Size(2) * H_idx)...
                            = 1;
                        probe_Mask = probe_Amplitude_Mask .* mod(SLM_BlazedGrating...
                            + mod(angle(exp(2 * pi * 1i * P_idx / Phase_Samples)) / 2 / pi,1)...
                            .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi,255);
                        full_Probe_Mask = ref_Mask + probe_Mask;
                        figure(SLMWindow);
                        axes(SLMAxes);
                        SLMObj = image(full_Probe_Mask);
                        drawnow;shg;
                        pause(1/60); % Short delay to allow SLM to update
                        % Acquire frame
                        if Flag_No_Camera == 0
                            trigger(vidObj);
                            image_Frame = getdata(vidObj,1);
                            intensity_Matrix(V_idx,H_idx,P_idx)...
                                = sum(sum(image_Frame(camYCoord - Correction_HalfROI:camYCoord + Correction_HalfROI...
                                ,camXCoord - Correction_HalfROI:camXCoord + Correction_HalfROI)));
                            if max(image_Frame(:)) >= 2^(Camera_Settings.pixelDepth - 1)
                                disp('Saturation warning.')
                            end
                            if live_Cam == 1
                                figure(40);
                                imagesc(image_Frame(camYCoord - 50:camYCoord + 50 ...
                                ,camXCoord - 50:camXCoord + 50));
                                axis image;
                                drawnow;shg;
                            end
                        else
                            intensity_Matrix(V_idx,H_idx,P_idx) = cos(2 * pi * P_idx / Phase_Samples) + 1 + 0.5 * rand(1); % sin function with random noise (used to test fitting algorithm)
                        end
                    end
                end
            end
        end
        measurement_Time = toc(measurement_Timer);
        disp(strcat('Wavefront correction measurement complete. Data acquisition took (',num2str(measurement_Time),') seconds.'))
        
        SLMObj = image(ones(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'uint8')...
            .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi);
        
        %%% Analysis of acquired data
        % Loop over all vertical mode coords
        for V_idx = 1:noVerticalModes
            % Loop over all horizontal mode coords
            for H_idx = 1:noHorizontalModes
                if V_idx == ref_Mode_Vertical && H_idx == ref_Mode_Horizontal
                    % No need to probe the reference mode itself
                    cosine_Fitting_Parameters(V_idx,H_idx,1) = max(max(cosine_Fitting_Parameters(:,:,1)));
                    cosine_Fitting_Parameters(V_idx,H_idx,2) = (2.^SLM_Settings.pixelDepth - 1) / SLM_Settings.modulation.Range2Pi;
                    cosine_Fitting_Parameters(V_idx,H_idx,3) = 2 * pi;
                    cosine_Fitting_Parameters(V_idx,H_idx,4) = cosine_Fitting_Parameters(V_idx,H_idx,1);
%                     cosine_Fitting_Parameters = zeros(noVerticalModes,noHorizontalModes,4,'single');
                        %y = A * cos(B * x + C) + D (4 parameters to fit).
                else
                    [max_val,max_idx] = max(intensity_Matrix(V_idx,H_idx,:));
                    min_val = min(intensity_Matrix(V_idx,H_idx,:));
                    cosine_Params0(1) = (max_val - min_val) / 2;
                    cosine_Params0(2) = (2.^SLM_Settings.pixelDepth - 1) / SLM_Settings.modulation.Range2Pi;
                    cosine_Params0(3) = max_idx / Phase_Samples * 2 * pi;
                    cosine_Params0(4) = max_val - cosine_Params0(1);
                    cosine_Fitting_Parameters(V_idx,H_idx,:) = fminsearch(@(cosine_Params) cosine_Parameter_Minimisation_Function(intensity_Matrix(V_idx,H_idx,:),cosine_Params,[1:Phase_Samples * Phase_Periods] * 2 * pi / Phase_Samples),cosine_Params0);
                end
                cosine_Fitting_Matrix(V_idx,H_idx,:) = cosine_Function(cosine_Fitting_Parameters(V_idx,H_idx,:),[1:Phase_Samples * Phase_Periods] * 2 * pi / Phase_Samples);
                
                Flag_Plot_Results = 0;
                if Flag_Plot_Results
                    figure(50);
                    plot([1:Phase_Samples * Phase_Periods] * 2 * pi / Phase_Samples,squeeze(intensity_Matrix(V_idx,H_idx,:)),[1:Phase_Samples * Phase_Periods] * 2 * pi / Phase_Samples,squeeze(cosine_Fitting_Matrix(V_idx,H_idx,:)));
                    title(strcat('V_idx = (',num2str(V_idx),') of (',num2str(noVerticalModes),'), H_idx = (',num2str(H_idx),') of (',num2str(noHorizontalModes),').'));
                    disp('Fitted function: y(t) = A * cos (B * t + C) + D. Parameters:')
                    disp(strcat('A = (',num2str(cosine_Fitting_Parameters(V_idx,H_idx,1)),').'))
                    disp(strcat('B = (',num2str(cosine_Fitting_Parameters(V_idx,H_idx,2)),').'))
                    disp(strcat('C = (',num2str(cosine_Fitting_Parameters(V_idx,H_idx,3)),').'))
                    disp(strcat('D = (',num2str(cosine_Fitting_Parameters(V_idx,H_idx,4)),').'))
                    pause();
                end
            end
        end
        
        %%% Process results
        SLM_Amplitude_Map = zeros(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'single');
        SLM_Phase_Map = zeros(SLM_Settings.size.slmHeight,SLM_Settings.size.slmWidth,'single');
        % Loop over all vertical mode coords
        for V_idx = 1:noVerticalModes
            % Loop over all horizontal mode coords
            for H_idx = 1:noHorizontalModes
                    SLM_Amplitude_Map(Probe_Mode_Size(1) * (V_idx - 1) + 1:Probe_Mode_Size(1) * V_idx...
                        ,Probe_Mode_Size(2) * (H_idx - 1) + 1:Probe_Mode_Size(2) * H_idx)...
                        = cosine_Fitting_Parameters(V_idx,H_idx,1);
                    SLM_Phase_Map(Probe_Mode_Size(1) * (V_idx - 1) + 1:Probe_Mode_Size(1) * V_idx...
                        ,Probe_Mode_Size(2) * (H_idx - 1) + 1:Probe_Mode_Size(2) * H_idx)...
                        = cosine_Fitting_Parameters(V_idx,H_idx,3);
            end
        end
        SLM_Phase_Map = mod(SLM_Phase_Map / 2 / pi,1)...
            .* (2.^SLM_Settings.pixelDepth - 1) .* SLM_Settings.modulation.Range2Pi;
        figure(51);
        imagesc(SLM_Amplitude_Map);
        axis image;
        title('amplitude map');
        figure(52);
        imagesc(SLM_Phase_Map);
        axis image;
        title('phase map');
        disp(strcat('2 Pi dynamic range at correction wavelength = (',num2str(round(mean(mean(cosine_Fitting_Parameters(:,:,2))))),').'))
        
        %%% 
        
    else
        disp('Error encountered during initialisation. Wavefront measurement not made.')
    end


    %%% De-initialise camera
    switch Camera_Type
        case 'Hamamatsu_Orca_Flash4.0_v2_C11440'
            stop(vidObj);
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