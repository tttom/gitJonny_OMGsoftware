%%% name:           equipartitionTheoremFromAvi
%%% author:         Jonathan Nylk
%%% date created:   05/07/2016
%%% description:    This function determines optical trap stiffness by
%%%                 equipartition theorem method. Give "videoName" of avi
%%%                 data file and set "user inputs" to tweak tracking
%%%                 algorithm parameters.
%%%                 Function used 3 methods to determine trap stiffness,
%%%                 "imfindcircles" function, centre-of-mass determination
%%%                 of thresholded image, and centre-of-mass determination
%%%                 of thresholded image that has been passed through a
%%%                 spatial low-pass filter.
%%%                 All traces and histograms are output figures, trap
%%%                 stiffness by each method is printed in the command
%%%                 window.
%%%
%%% updates (latest first):
%%%                 11/07/2016: Made code more robust to errors in circle
%%%                 finding algorithm. If no circles are found, this entry
%%%                 is filled as 'Nan', then as the mean value of the
%%%                 non-NaN entries.
%%%                 05/07/2016: Added drift correction by low-order
%%%                 polynomial fitting of particle traces.
%%%
%%%
%%% END %%%


% constants
    k_b = 1.38e-23; % Boltzmann constant in m^2*kg*s^-2*K^-1
    T = 23 + 273; % (kelvin)

% read video info
    videoName = 'E:\Stored Files\Frances Trapping Data\2016-06-20_IlluminationLevelSetting\Trapped bead\24mW trans beam a.avi';
    vidObj = VideoReader(videoName);
    imHeight = vidObj.height;
    imWidth = vidObj.width;
    noFrames = vidObj.NumberOfFrames - 1;
%     noFrames = 400; % for testing
    showFrames = randi(noFrames,[40,1]);
    showFramesCounter = 1;

% user inputs
    imageScalingFactor = 2;
    camera_pixel_size = [7.4e-6,7.4e-6]; %[y,x] (metres)
    magnification = 100;
    bead_diameter = 2.19e-6; %(metres)
    trap_power = 24e-3; %(watts)
    % imfindcircles
        min_radius_limit = floor(6 * imageScalingFactor);
        max_radius_limit = ceil(15 * imageScalingFactor);
    % centre-of-mass
        image_threshold_value = 0.6;
    % fourier filtered centre-of-mass
        filter_type = 'Circle';
    %     filter_type = 'Gaussian';
        filter_strength = 0.15 * round(imWidth * imageScalingFactor /2);
        filtered_image_threshold_value = 0.6;
    % drift compensation
        flag_correct_drift = 1; % 1: yes, 0: no
        drift_polynomial_order = 3;


% allocate memory for results
currentFrame = zeros([imHeight,imWidth,3],'uint8');
currentImage = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
% imfindcircles
    xPosition_pixels = zeros([noFrames,1]);
    yPosition_pixels = zeros([noFrames,1]);
    xPosition = zeros([noFrames,1]);
    yPosition = zeros([noFrames,1]);
    radius_pixels = zeros([noFrames,1]);
    number_circles_found = zeros([noFrames,1]);
% centre-of-mass
    currentImageInverted = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    thresholded_image = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    [xPixels,yPixels] = meshgrid([1:imWidth * imageScalingFactor] - round(imWidth * imageScalingFactor /2)...
        ,[1:imHeight * imageScalingFactor] - round(imHeight * imageScalingFactor /2));
    xPosition_pixels_COM = zeros([noFrames,1]);
    yPosition_pixels_COM = zeros([noFrames,1]);
    xPosition_COM = zeros([noFrames,1]);
    yPosition_COM = zeros([noFrames,1]);
% fourier filtered centre-of-mass
    fourier_filter = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    if strcmp(filter_type,'Circle')
        fourier_filter = xPixels.^2 + yPixels.^2 <= filter_strength^2;
    elseif strcmp(filter_type,'Gaussian')
        fourier_filter = exp(-(x.^2 + yPixels.^2) / 2 / filter_strength^2);
    end
    fourier_image = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    filtered_fourier_image = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    filtered_image = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    thresholded_filtered_image = zeros([imHeight * imageScalingFactor,imWidth * imageScalingFactor],'single');
    xPosition_pixels_fCOM = zeros([noFrames,1]);
    yPosition_pixels_fCOM = zeros([noFrames,1]);
    xPosition_fCOM = zeros([noFrames,1]);
    yPosition_fCOM = zeros([noFrames,1]);

start_timer = tic;

for fIdx = 1:noFrames
    currentFrame = read(vidObj,fIdx+1); %1st frame is spacer
    currentImage = imresize(single(currentFrame(:,:,1)),imageScalingFactor);
    currentImage = currentImage/max(currentImage(:)); %normalise
    
    % imfindcircles
        [centres, radii] = imfindcircles(currentImage,[min_radius_limit max_radius_limit]...
            ,'Sensitivity', 0.9, 'ObjectPolarity', 'bright'); %finding the circle fit

        if size(radii,1)>0
            centre = centres(1,:);
            radius = radii(1,:);
        else
            centre = [NaN,NaN];
            radius = NaN;
        end

        xPosition_pixels(fIdx,1) = centre(1,1);
        yPosition_pixels(fIdx,1) = centre(1,2);
        radius_pixels(fIdx,1) = radius(1,1);
        number_circles_found(fIdx,1) = size(radii,1);

%         figure(96);
%         imagesc(currentImage);axis image;colormap gray;
%         viscircles(centres,radii,'EdgeColor','r');
%         viscircles(centre,radius,'EdgeColor','b');
%         title(sprintf('imfindcircles: Frame %d of %d',[fIdx,noFrames]));
%         drawnow;shg;
% %         pause(0.1);

    % centre-of-mass
        currentImageInverted = 1 - currentImage;
        thresholded_image = currentImageInverted .* (currentImageInverted >= image_threshold_value);
        xPosition_pixels_COM(fIdx,1) = sum(sum(thresholded_image .* xPixels)) / sum(thresholded_image(:));
        yPosition_pixels_COM(fIdx,1) = sum(sum(thresholded_image .* yPixels)) / sum(thresholded_image(:));
%         
%         figure(97);
%         subplot(1,3,1);
%         imagesc(currentImageInverted);axis image;colormap gray;
%         title('starting image');
%         subplot(1,3,2);
%         imagesc(thresholded_image);axis image;colormap gray;
%         title('thresholded image');
%         subplot(1,3,3);
%         scatter(xPosition_pixels_COM(fIdx,1),yPosition_pixels_COM(fIdx,1),5);axis square;
%         xlim([xPixels(1) xPixels(end)]);
%         ylim([yPixels(1) yPixels(end)]);
%         title(sprintf('CoM: Frame %d of %d',[fIdx,noFrames]));
%         drawnow;shg;
% %         pause(0.1);
        
    % fourier filtered centre-of-mass
        fourier_image = fftshift(fft2(currentImageInverted));
        filtered_fourier_image = fourier_image .* fourier_filter;
        filtered_image = real(ifft2(ifftshift(filtered_fourier_image)));
        filtered_image = filtered_image / max(filtered_image(:));
        thresholded_filtered_image = filtered_image .* (filtered_image >= filtered_image_threshold_value);
        xPosition_pixels_fCOM(fIdx,1) = sum(sum(thresholded_filtered_image .* xPixels)) / sum(thresholded_filtered_image(:));
        yPosition_pixels_fCOM(fIdx,1) = sum(sum(thresholded_filtered_image .* yPixels)) / sum(thresholded_filtered_image(:));
    
%         figure(98);
%         subplot(2,3,1);
%         imagesc(currentImageInverted);axis image;colormap gray;
%         title('starting image');
%         subplot(2,3,2);
%         imagesc(filtered_image);axis image;colormap gray;
%         title('filtered image');
%         subplot(2,3,5);
%         imagesc(log(abs(fourier_image)));axis image;colormap gray;
%         title('fourier image (log)');
%         subplot(2,3,4);
%         imagesc(log(abs(filtered_fourier_image)));axis image;colormap gray;
%         title('filtered fourier image (log)');
%         subplot(2,3,3);
%         imagesc(thresholded_filtered_image);axis image;colormap gray;
%         title('thresholded filtered image');
%         subplot(2,3,6);
%         scatter(xPosition_pixels_fCOM(fIdx,1),yPosition_pixels_fCOM(fIdx,1),5);axis square;
%         xlim([xPixels(1) xPixels(end)]);
%         ylim([yPixels(1) yPixels(end)]);
%         title(sprintf('filtered CoM: Frame %d of %d',[fIdx,noFrames]));
%         drawnow;shg;
% %         pause(0.1);
        
    % loop progress
    if rem(fIdx,round(noFrames/10))==0
        fprintf('Processed %d of %d frames \n',[fIdx,noFrames]);
    end

%     % show 40 random video frames and their circle fits
%     if max(fIdx==showFrames)
%         figure(99);
%         subplot(5,8,showFramesCounter);
%         imagesc(currentImage);axis image;colormap gray;
%         viscircles(centres,radii,'EdgeColor','r');
%         viscircles(centre,radius,'EdgeColor','b');
%         title(sprintf('Frame %d',fIdx));
%         drawnow;shg;
%         showFramesCounter = showFramesCounter + 1;
%     end
end


NaN_counter = 0;
if max(isnan(xPosition_pixels))
    nonNaN_mean = nanmean(xPosition_pixels);
    NaN_positions = isnan(xPosition_pixels);
    NaN_counter = sum(NaN_positions);
    xPosition_pixels(isnan(xPosition_pixels)) = 0;
    xPosition_pixels = ((1 - NaN_positions) .* xPosition_pixels) + (NaN_positions * nonNaN_mean);
end
if max(isnan(yPosition_pixels))
    nonNaN_mean = nanmean(yPosition_pixels);
    NaN_positions = isnan(yPosition_pixels);
    yPosition_pixels(isnan(yPosition_pixels)) = 0;
    yPosition_pixels = ((1 - NaN_positions) .* yPosition_pixels) + (NaN_positions * nonNaN_mean);
end


% include option to fit linear function to x and y traces to account for
% drift
if flag_correct_drift
    % imfindcircles
    x_polyfit_params = polyfit([1:length(xPosition_pixels)].',squeeze(xPosition_pixels),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(xPosition_pixels)],xPosition_pixels,'b'...
        ,[1:length(xPosition_pixels)],polyval(x_polyfit_params,[1:length(xPosition_pixels)]),'r');
    title('x-axis trace (before drift correction) (imfindcircles');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    xPosition_pixels = xPosition_pixels - polyval(x_polyfit_params,[1:length(xPosition_pixels)]).';
    subplot(2,1,2);plot([1:length(xPosition_pixels)],xPosition_pixels,'b');
    title('x-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    y_polyfit_params = polyfit([1:length(yPosition_pixels)].',squeeze(yPosition_pixels),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(yPosition_pixels)],yPosition_pixels,'g'...
        ,[1:length(yPosition_pixels)],polyval(y_polyfit_params,[1:length(yPosition_pixels)]),'r');
    title('y-axis trace (before drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
    yPosition_pixels = yPosition_pixels - polyval(y_polyfit_params,[1:length(yPosition_pixels)]).';
    subplot(2,1,2);plot([1:length(yPosition_pixels)],yPosition_pixels,'g');
    title('y-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
    
    % centre-of-mass
    x_polyfit_params_COM = polyfit([1:length(xPosition_pixels_COM)].',squeeze(xPosition_pixels_COM),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(xPosition_pixels_COM)],xPosition_pixels_COM,'b'...
        ,[1:length(xPosition_pixels_COM)],polyval(x_polyfit_params_COM,[1:length(xPosition_pixels_COM)]),'r');
    title('x-axis trace (before drift correction) (centre-of-mass)');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    xPosition_pixels_COM = xPosition_pixels_COM - polyval(x_polyfit_params_COM,[1:length(xPosition_pixels_COM)]).';
    subplot(2,1,2);plot([1:length(xPosition_pixels_COM)],xPosition_pixels_COM,'b');
    title('x-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    y_polyfit_params_COM = polyfit([1:length(yPosition_pixels_COM)].',squeeze(yPosition_pixels_COM),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(yPosition_pixels_COM)],yPosition_pixels_COM,'g'...
        ,[1:length(yPosition_pixels_COM)],polyval(y_polyfit_params_COM,[1:length(yPosition_pixels_COM)]),'r');
    title('y-axis trace (before drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
    yPosition_pixels_COM = yPosition_pixels_COM - polyval(y_polyfit_params_COM,[1:length(yPosition_pixels_COM)]).';
    subplot(2,1,2);plot([1:length(yPosition_pixels_COM)],yPosition_pixels_COM,'g');
    title('y-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
    
    % filtered centre-of-mass
    x_polyfit_params_fCOM = polyfit([1:length(xPosition_pixels_fCOM)].',squeeze(xPosition_pixels_fCOM),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(xPosition_pixels_fCOM)],xPosition_pixels_fCOM,'b'...
        ,[1:length(xPosition_pixels_fCOM)],polyval(x_polyfit_params_fCOM,[1:length(xPosition_pixels_fCOM)]),'r');
    title('x-axis trace (before drift correction) (filtered centre-of-mass)');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    xPosition_pixels_fCOM = xPosition_pixels_fCOM - polyval(x_polyfit_params_fCOM,[1:length(xPosition_pixels_fCOM)]).';
    subplot(2,1,2);plot([1:length(xPosition_pixels_fCOM)],xPosition_pixels_fCOM,'b');
    title('x-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('x-axis [pixels]');
    y_polyfit_params_fCOM = polyfit([1:length(yPosition_pixels_fCOM)].',squeeze(yPosition_pixels_fCOM),drift_polynomial_order);
    figure();subplot(2,1,1);plot([1:length(yPosition_pixels_fCOM)],yPosition_pixels_fCOM,'g'...
        ,[1:length(yPosition_pixels_fCOM)],polyval(y_polyfit_params_fCOM,[1:length(yPosition_pixels_fCOM)]),'r');
    title('y-axis trace (before drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
    yPosition_pixels_fCOM = yPosition_pixels_fCOM - polyval(y_polyfit_params_fCOM,[1:length(yPosition_pixels_fCOM)]).';
    subplot(2,1,2);plot([1:length(yPosition_pixels_fCOM)],yPosition_pixels_fCOM,'g');
    title('y-axis trace (after drift correction)');
    xlabel('frame #');
    ylabel('y-axis [pixels]');
end
    
% mean-centre data - imfindcircles
x_mean = mean(xPosition_pixels);
y_mean = mean(yPosition_pixels);
xPosition_pixels_meanSubtracted = xPosition_pixels - x_mean;
yPosition_pixels_meanSubtracted = yPosition_pixels - y_mean;

% mean-centre data - CoM
x_mean_COM = mean(xPosition_pixels_COM);
y_mean_COM = mean(yPosition_pixels_COM);
xPosition_pixels_meanSubtracted_COM = xPosition_pixels_COM - x_mean_COM;
yPosition_pixels_meanSubtracted_COM = yPosition_pixels_COM - y_mean_COM;

% mean-centre data - filtered CoM
x_mean_fCOM = mean(xPosition_pixels_fCOM);
y_mean_fCOM = mean(yPosition_pixels_fCOM);
xPosition_pixels_meanSubtracted_fCOM = xPosition_pixels_fCOM - x_mean_fCOM;
yPosition_pixels_meanSubtracted_fCOM = yPosition_pixels_fCOM - y_mean_fCOM;

% convert pixels to real-space (metres) - imfindcircles
xPosition = xPosition_pixels_meanSubtracted * camera_pixel_size(2) / imageScalingFactor / magnification;
yPosition = yPosition_pixels_meanSubtracted * camera_pixel_size(1) / imageScalingFactor / magnification;

% convert pixels to real-space (metres) - CoM
xPosition_COM = xPosition_pixels_meanSubtracted_COM * camera_pixel_size(2) / imageScalingFactor / magnification;
yPosition_COM = yPosition_pixels_meanSubtracted_COM * camera_pixel_size(1) / imageScalingFactor / magnification;

% convert pixels to real-space (metres) - filtered CoM
xPosition_fCOM = xPosition_pixels_meanSubtracted_fCOM * camera_pixel_size(2) / imageScalingFactor / magnification;
yPosition_fCOM = yPosition_pixels_meanSubtracted_fCOM * camera_pixel_size(1) / imageScalingFactor / magnification;

    % output plots
    % traces
    % imfindcircles
        figure();
        subplot(2,1,1);plot(xPosition * 1e6,'b');
        title('x-axis trace (imfindcircles)');
        xlabel('frame #');
        ylabel('x-axis [um]');
        subplot(2,1,2);plot(yPosition * 1e6,'g');
        title('y-axis trace (imfindcircles)');
        xlabel('frame #');
        ylabel('y-axis [um]');
    
    % centre-of-mass
        figure();
        subplot(2,1,1);plot(xPosition_COM * 1e6,'b');
        title('x-axis trace (centre-of-mass)');
        xlabel('frame #');
        ylabel('x-axis [um]');
        subplot(2,1,2);plot(yPosition_COM * 1e6,'g');
        title('y-axis trace (centre-of-mass)');
        xlabel('frame #');
        ylabel('y-axis [um]');
    
    % filtered centre-of-mass
        figure();
        subplot(2,1,1);plot(xPosition_fCOM * 1e6,'b');
        title('x-axis trace (filtered centre-of-mass)');
        xlabel('frame #');
        ylabel('x-axis [um]');
        subplot(2,1,2);plot(yPosition_fCOM * 1e6,'g');
        title('y-axis trace (filtered centre-of-mass)');
        xlabel('frame #');
        ylabel('y-axis [um]');
        
    %histograms
    %imfindcircles
        figure();
        subplot(1,2,1);hist(xPosition * 1e6,50);
        title('x-axis histogram (imfindcircles)');
        xlabel('x-axis [um]');
        ylabel('frequency');
        subplot(1,2,2);hist(yPosition * 1e6,50);
        title('y-axis histogram (imfindcircles)');
        xlabel('y-axis [um]');
        ylabel('frequency');
        hist_handle = findobj(gca,'Type','patch');
        set(hist_handle,'FaceColor','g');
        
    %centre-ofmass
        figure();
        subplot(1,2,1);hist(xPosition_COM * 1e6,50);
        title('x-axis histogram (centre-of-mass)');
        xlabel('x-axis [um]');
        ylabel('frequency');
        subplot(1,2,2);hist(yPosition_COM * 1e6,50);
        title('y-axis histogram (centre-of-mass)');
        xlabel('y-axis [um]');
        ylabel('frequency');
        hist_handle = findobj(gca,'Type','patch');
        set(hist_handle,'FaceColor','g');
        
    %filtered centre-ofmass
        figure();
        subplot(1,2,1);hist(xPosition_fCOM * 1e6,50);
        title('x-axis histogram (filtered centre-of-mass)');
        xlabel('x-axis [um]');
        ylabel('frequency');
        subplot(1,2,2);hist(yPosition_fCOM * 1e6,50);
        title('y-axis histogram (filtered centre-of-mass)');
        xlabel('y-axis [um]');
        ylabel('frequency');
        hist_handle = findobj(gca,'Type','patch');
        set(hist_handle,'FaceColor','g');
        
    
% variance - imfindcircles
var_x = var(xPosition);
var_y = var(yPosition);

% trap stiffness - imfindcircles
k_x = k_b * T / (var_x);
k_y = k_b * T / (var_y);

% variance - CoM
var_x_COM = var(xPosition_COM);
var_y_COM = var(yPosition_COM);

% trap stiffness - CoM
k_x_COM = k_b * T / (var_x_COM);
k_y_COM = k_b * T / (var_y_COM);

% variance - filtered CoM
var_x_fCOM = var(xPosition_fCOM);
var_y_fCOM = var(yPosition_fCOM);

% trap stiffness - filtered CoM
k_x_fCOM = k_b * T / (var_x_fCOM);
k_y_fCOM = k_b * T / (var_y_fCOM);


time_elapsed = toc(start_timer);

% outputs
fprintf('Trapping video analysed in %d s. \n',round(time_elapsed));
fprintf('Video parameters: \n');
fprintf('Bead diamter (um) = %d. \n',bead_diameter * 1e6);
fprintf('Trap power (mW) = %d. \n',trap_power * 1e3);
fprintf('Number of frames = %d. \n',noFrames);
fprintf('Original image size (pixels) = [%d, %d]. \n',[imHeight,imWidth]);
fprintf('Size of camera pixels (um) = [%d,%d]. \n',[camera_pixel_size(1),camera_pixel_size(2)]);
fprintf('Microscope magnification = %d. \n',magnification);
fprintf('Pixel upsampling rate = %d. \n',imageScalingFactor);
if flag_correct_drift
    fprintf('Drift correction used.\n %d th-order polynomial fitted. \n',drift_polynomial_order);
else
    fprintf('No drift correction used. \n');
end
fprintf('Number of frames where no circle could be found using "imfindcircles": %d. \n\n',NaN_counter);

% fprintf('Determining bead size based on particle tracking and estimated magnification: \n');
% fprintf('Mean paricle diamter (um) = %d. \n',mean(radius_pixels) * camera_pixel_size(2) / imageScalingFactor / magnification *1e6 * 2);
% fprintf('Determining magnification based on particle tracking and estimated bead size: \n');
% fprintf('Magnification = %d. \n \n',mean(radius_pixels) * camera_pixel_size(2) / imageScalingFactor / bead_diameter * 2);
% 
fprintf('Using "imfindcircles" method: \n');
fprintf('Trap stiffness along x-axis is %d pN/um. \n',k_x * 1e12 / 1e6);
fprintf('Trap stiffness along y-axis is %d pN/um. \n \n',k_y * 1e12 / 1e6);

fprintf('Using CoM method: \n');
fprintf('Trap stiffness along x-axis is %d pN/um. \n',k_x_COM * 1e12 / 1e6);
fprintf('Trap stiffness along y-axis is %d pN/um. \n \n',k_y_COM * 1e12 / 1e6);

fprintf('Using filtered CoM method: \n');
fprintf('Trap stiffness along x-axis is %d pN/um. \n',k_x_fCOM * 1e12 / 1e6);
fprintf('Trap stiffness along y-axis is %d pN/um. \n \n',k_y_fCOM * 1e12 / 1e6);
fprintf('\n \n');

