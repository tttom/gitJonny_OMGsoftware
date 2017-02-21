%%% File: CalciumSignallingAnalysis.m
%%% Author: Jonathan Nylk
%%% Last Updated: 12/05/2015
%%%
%%% This program is part of a body of work that is published and described 
%%% in Nylk, J. et al, "Development of a graded index microlens based fiber
%%% optical trap and its characterization using principal component
%%% analysis", Biomedical Optics Express 6(4) 1512-1519 (2015) doi:
%%% 10.1364/BOE.6.001512.
%%%
%%% This program determines the total intensity resulting from 2 cells
%%% (designated "Signal" and "Control") from image frames aquired over time
%%% during an experiment bringing a B3Z T-cell hybridoma into contact with
%%% an immobilised antibody-coated bead with a GRIN lens based fibre
%%% optical trap. Photobleaching is ocrrected for by monitoring the
%%% intensity of the Control cell.
%%% Due to unwanted motion during the experiment, the centre of mass (CoM)
%%% of each cell is tracked and used to dynamically update the image
%%% regions associated with the Signal and Control cells.
%%% This program outputs two figures. The first contains a plot of intensity
%%% for the Signal and Control cells with an exponential decay curve fitted
%%% to Control. The second contains 3 cropped image frames from the
%%% experiment (one before cell-bead contact, one immediately after
%%% contact, and one much later after contact) which are also shown in the
%%% paper.

%% Initialisation
    % Folder location of image files and base filename with time-index removed
    baseFileName='...\ActivationImageFrames\bead_b3z_activation';

    % Define starting frame (frame when cell-bead contact is initiated)
    totalFrames=3256;
    frameStart=250;
    nFrames=totalFrames-frameStart;

    % Allocate memory for variables and results
    [x,y]=meshgrid(1:512,1:512); %coordinate system (size of image in pixels)
    signal=zeros(1,nFrames);
    control=zeros(1,nFrames);
    background=zeros(1,nFrames);
    signalCoMTracker_x=zeros(1,nFrames);
    signalCoMTracker_y=zeros(1,nFrames);
    controlCoMTracker_x=zeros(1,nFrames);
    controlCoMTracker_y=zeros(1,nFrames);
    frame=zeros(512,512);
    selectedFrames=zeros(512,512,3);
    m=1; %counter for "selectedFrames" selection

%% Track CoM coordinates to acocunt for sample movement
    for n=1:nFrames
        % Generate full filename for a time-indexed image and load
        fileNum=10000+n+frameStart;
        fileNum=num2str(fileNum);
        fileNum=fileNum(2:end);
        fileName=strcat(baseFileName,fileNum,'.tif');
        frame=double(imread(fileName));

        % Determine CoM for Signal
            % Select area of image only containing Signal
            cropRegion=frame(281:400,ceil(end/2):end);
            % rescale and threshold image to remove background
            cropRegion=cropRegion./max(cropRegion(:));
            cropRegion=cropRegion.*(cropRegion>0.5);
            %determine COM x- and y-coordinates for Signal
            signalCoMTracker_x(1,n)=sum(sum(cropRegion.*x(281:400,ceil(end/2):end)))./sum(cropRegion(:));
            signalCoMTracker_y(1,n)=sum(sum(cropRegion.*y(281:400,ceil(end/2):end)))./sum(cropRegion(:));
        % Determine CoM for Control (same as above)
            cropRegion=frame(281:400,1:floor(end/2));
            cropRegion=cropRegion./max(cropRegion(:));
            cropRegion=cropRegion.*(cropRegion>0.5);
            controlCoMTracker_x(1,n)=sum(sum(cropRegion.*x(281:400,1:floor(end/2))))./sum(cropRegion(:));
            controlCoMTracker_y(1,n)=sum(sum(cropRegion.*y(281:400,1:floor(end/2))))./sum(cropRegion(:));

        % Select frames to be displayed in paper to show the cases:
        % 1: before cell-bead contact
        % 2: immediately after cell-bead contact
        % 3: much time after cell-bead contact
        if max(n==[10,200,1000])==1
            selectedFrames(:,:,m)=frame;
            m=m+1; %increase counter
        end
        
        % Progress update to user
        disp(strcat('frame: (',num2str(n),') of (',num2str(nFrames),')'))
    end

    % Round CoM coordinates to nearest pixel
    signalCoMTracker_x=round(signalCoMTracker_x);
    signalCoMTracker_y=round(signalCoMTracker_y);
    controlCoMTracker_x=round(controlCoMTracker_x);
    controlCoMTracker_y=round(controlCoMTracker_y);

%% Integrate intensity values within Signal, Control, and Background regions
    for n=1:nFrames
        % Generate full filename for a time-indexed image and load
        fileNum=10000+n+frameStart;
        fileNum=num2str(fileNum);
        fileNum=fileNum(2:end);
        fileName=strcat(baseFileName,fileNum,'.tif');
        frame=double(imread(fileName));    

        % Define spatial regions for Signal, Control, and Background (all
        % pixels not within Signal or Control regions), dynamically updated
        % based on CoM tracking.
        signalFilter=(((x-signalCoMTracker_x(1,n)).^2)+((y-signalCoMTracker_y(1,n)).^2))<=12^2;
        controlFilter=(((x-controlCoMTracker_x(1,n)+4).^2)+((y-controlCoMTracker_y(1,n)-2).^2))<=8^2;
        backgroundFilter=min((((x-controlCoMTracker_x(1,n)).^2)+((y-controlCoMTracker_y(1,n)).^2)>=50^2),(((x-signalCoMTracker_x(1,n)).^2)+((y-signalCoMTracker_y(1,n)).^2)>=50^2));

        %Isolate Signal region and integrate over all pixels
            signalFrame=frame.*signalFilter;
            signal(1,n)=sum(signalFrame(:));
        %Isolate Control region and integrate over all pixels
            controlFrame=frame.*controlFilter;
            control(1,n)=sum(controlFrame(:));
        %Isolate Background region and integrate over all pixels
        backgroundFrame=frame.*backgroundFilter;
        background(1,n)=sum(backgroundFrame(:));
        
        % Progress update to user
        disp(strcat('frame: (',num2str(n),') of (',num2str(nFrames),')'))

    end
    
    % Normalise values by intensity per pixel for comparison
    signal=signal./sum(signalFilter(:));
    control=control./sum(controlFilter(:));
    background=background./sum(backgroundFilter(:));

    % Subtract background level
    signal=signal-background;
    control=control-background;

    % Normalise control (in arbitrary units)
    control=control./control(1,1);
    % Normalise Signal to Control value at corresponding time-index to
    % account for loss of intensity by photobleaching
    signal=signal./control;

    % Average frames to smooth traces
        frames2Average=3; %number of frames to average
        %allocate memory for temporary Signal and Control variables
        signal2=zeros(1,nFrames-(frames2Average-1));
        control2=zeros(1,nFrames-(frames2Average-1));
        for n=1:nFrames-(frames2Average-1);
            signal2(n)=mean(signal(n:n+frames2Average-1));
            control2(n)=mean(control(n:n+frames2Average-1));
        end
    % Overwrite Signal and Control variables with smoothed traces
    signal=signal2;
    control=control2;
    nFrames=nFrames-(frames2Average-1); %reduce number of elements by frames2Average-1 when averaging

%% Exponential Fitting
    % Fit a single exponential function to Control
        fittedFunct=fit((1:nFrames)'*50e-3,(control./max(control(:)))','exp1');
        fittedData=fittedFunct((1:nFrames)*50e-3); %50ms (50e-3s) is framerate for experiment

%% Plot Results
    % Plot intensity of Signal and Control as function of time
    figure(1);[Ax,H1,H2]=plotyy((1:nFrames)*50e-3,signal./max(signal(:)),(1:nFrames)*50e-3,control./max(control(:))); %50ms (50e-3s) is framerate for experiment
    drawnow;axes(Ax(1));
    xlim([0 150]);ylim([.1 1.1]);
    ylabel('Calcium Signal [a.u.]');
    drawnow;axes(Ax(2));
    xlim([0 150]);ylim([.7 1.1]);
    xlabel('Time [s]');
    ylabel('Calcium Control [a.u.]');
    legend('Control','Signal','location','northwest'); %Exponential fit is added to legend manually after exporting to inkscape
    hold on
    plot((1:nFrames)*50e-3,fittedData,'r');
    hold off

    % Display cropped versions of the selected frames for use in
    % publication
    % Normalise and scale images to use full dynamic range
        selectedFrames=selectedFrames./max(selectedFrames(:));
        scaleFactor=64;
    figure(2);
    subplot(1,3,1);image((1:512)*.11,(1:512)*.11,scaleFactor*squeeze(selectedFrames(:,:,1))); %.11 is for .11um/pixel scale factor
    axis image;xlim([20 400]*.11);ylim([250 450]*.11);colormap gray;
    subplot(1,3,2);image((1:512)*.11,(1:512)*.11,scaleFactor*squeeze(selectedFrames(:,:,2)));
    axis image;xlim([20 400]*.11);ylim([250 450]*.11);colormap gray;
    subplot(1,3,3);image((1:512)*.11,(1:512)*.11,scaleFactor*squeeze(selectedFrames(:,:,3)));
    axis image;xlim([20 400]*.11);ylim([250 450]*.11);colormap gray;

    clear all %clear all variables from memory

%%% End