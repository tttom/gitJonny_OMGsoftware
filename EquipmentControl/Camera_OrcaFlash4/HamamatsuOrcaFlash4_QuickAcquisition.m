%%% Initialise camera
    disp('Initialising camera...')
    vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode');
    triggerconfig(vid,'manual')
    src = getselectedsource(vid);
    src.ExposureTime = 0.01; %in seconds (I think)
    vid.FramesPerTrigger = 1;
    vid.TriggerRepeat = Inf;
    imageSize = vid.VideoResolution;
    imHeight = imageSize(1);
    imWidth = imageSize(2);
    start(vid)
    pause(1)
    disp('Camera initialisation complete')
    %%%
    
    disp('Setting other parameters...')
    %%% Preview camera
    % Acquire frame
    previewFrame = zeros(imHeight,imWidth,'uint16');
    trigger(vid);
    previewFrame = getdata(vid,1);
    previewFig = figure();imagesc(previewFrame);axis image;
    title('Click on point in image to perform correction on to continue')
    [~,~] = ginput(1);
    close(previewFig)
    %%%
    
    nFrames = 1000;
    frames = zeros(imHeight,imWidth,nFrames,'uint16');
    
    tic
    for n = 1:nFrames
        % Acquire frame
        trigger(vid);
        frames(:,:,n) = getdata(vid,1);
    end
    timeElapsed = toc
    averageFrameRate = nFrames / timeElapsed
    
%     for n = 1:nFrames
%         %%% Camera display
%         figure(100);
%         imagesc(frames(:,:,n));axis image;
%         title(num2str(n));
%         drawnow;shg;
%         %%%
%     end
    
    delete(vid);
    clear vid src;