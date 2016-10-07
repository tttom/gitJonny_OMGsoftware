%%% HEADER %%%
%%% Simon to populate
%%%

%%% Simulates the (xz) light-sheet imaging process on a user defined image. 
%%% All imaging parameters are taken from processed data ".mat" files.
%%% 
%%% Inputs: 

function simulate_light_sheet_imaging(image_target_filename,imaging_parameters_filename,output_folder,z_interpolation_factor,invert_image)

    % Default inputs.
    if nargin < 5
        invert_image = 1;
    end
    
    if nargin < 4
%         z_interpolation_factor = 0;
%         z_interpolation_factor = 1;
        z_interpolation_factor = 2;
    end
    
    if nargin < 3
        output_folder = ...
            'G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\M2 Email Banner';
    end
    
    if nargin < 2
        imaging_parameters_filename = ...
            'G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\__Acquired Data\recording0_lambda532nm_alpha0_beta100.mat';
    end
    
    if nargin < 1
%         image_target_filename = ...
%             'C:\Users\Jonathan Nylk\Dropbox\03-foundation-colour.bmp';
        image_target_filename = ...
            'G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\M2 Email Banner\m2_email_banner.jpg';
    end
    
    % Additional user specificed variables.
    convolution_offset = 0;    % pixels.
    noise_magnitude = 0.01;
    
    % Load imaging parameters.
    load(imaging_parameters_filename,'lightSheetPsf','yRange','zRange','setupConfig');
    
    % Calculate key parameters to determine image scaling.
    numerical_aperture_of_illuminating_lightsheet = ...
        setupConfig.excitation.objective.numericalAperture.*setupConfig.excitation.fractionOfNumericalApertureUsed;
    wavelength_of_lightsheet = setupConfig.excitation.wavelength;
    lightsheet_axial_resolution = wavelength_of_lightsheet / 2 / numerical_aperture_of_illuminating_lightsheet * 1e6;
    x_spacing = (yRange(2) - yRange(1)) * 1e6;
    z_spacing = (zRange(2) - zRange(1)) * 1e6;
    recommended_interpolation_factor = ceil(lightsheet_axial_resolution/z_spacing);
    x_width_pixels = size(yRange,2);
    x_width_um = x_width_pixels * x_spacing;
    x_start = yRange(1) * 1e6;
    x_end = yRange(end) * 1e6; %#ok<COLND>
    z_width_pixels = size(zRange,2);
    z_width_um = size(zRange,2) * z_spacing;
    z_start = zRange(1) * 1e6;
    z_end = zRange(end) * 1e6; %#ok<COLND>
    
    % Display key parameters to determine image scaling.
    fprintf('Light-sheet calculated over %g um [%d pixels] (x-axis), between %g : %g um, \n',x_width_um,x_width_pixels,x_start,x_end)
    fprintf('and %g um [%d pixels] (z-axis), between %g : %g um. \n',z_width_um,z_width_pixels,z_start,z_end)
    fprintf('Numerical aperture of light-sheet is %g, giving axial resolution of %g um. \n',numerical_aperture_of_illuminating_lightsheet,lightsheet_axial_resolution)
    fprintf('Z-plane spacing is %g um. \n',z_spacing)
    fprintf('Recommend interpolating the z-axis of image by %d times. \n',recommended_interpolation_factor)
    fprintf('Current interpolation set to %d times. \n',z_interpolation_factor)
    
    % Load simulation image target.
    base_image_target = single(imread(image_target_filename));
    for m = 1:3
        interp_image_target(:,:,m) = interp2(base_image_target(:,:,m),z_interpolation_factor);
    end
    
    % Normalise image target.
    image_target = interp_image_target ./ max(interp_image_target(:));
    
    % Check image target size matches light-sheet PSF.
    % (Width).
    if size(image_target,2) < x_width_pixels
        % Zero pad image.
        full_image = zeros([size(image_target,1),x_width_pixels,size(image_target,3)],'single');
        full_image(:,1:size(image_target,2),:) = image_target;
        image_target = full_image;
        clear full_image;
    else
        % Crop image target.
        image_target = image_target(:,1:x_width_pixels,:);
    end
    
    % (Height).
    if size(image_target,1) < z_width_pixels
        % Zero pad image.
        full_image = zeros([z_width_pixels,size(image_target,2),size(image_target,3)],'single');
        full_image(1:size(image_target,1),:,:) = image_target;
        image_target = full_image;
        clear full_image;
    else
        % Crop image target.
        image_target = image_target(:,1:x_width_pixels,:);
    end
    
    
    image_target_DISPLAY = image_target;
    % Invert image?
    if invert_image
        image_target = 1 - image_target;
    end
    
    % Convolve image target with light-sheet PSF.
    convolved_image=zeros([size(conv(squeeze(image_target(:,1,1)),squeeze(lightSheetPsf(1,1,:))),1),size(image_target,2),size(image_target,3)],'single');
    for m = 1:size(image_target,3)
        for n = 1:size(image_target,2)
            convolved_image(:,n,m) = conv(squeeze(image_target(:,n,m)),squeeze(flipdim(lightSheetPsf(1,n,:),3)));
        end
    end
    
    % Crop convolved image.
    convolved_image = convolved_image(floor(z_width_pixels / 2) + convolution_offset:floor(z_width_pixels / 2) + convolution_offset - 1 + z_width_pixels,:,:);
    
    % Add Gaussian noise.
    image_noise = randn(size(convolved_image)) .* max(convolved_image(:)) * noise_magnitude;
    convolved_image = convolved_image + image_noise;
    
    % Re-normalisation of image.
    if min(convolved_image(:)) < 0
        convolved_image = convolved_image - min(convolved_image(:));    % Ensures all +ve values
    end
    
    convolved_image = convolved_image ./ max(convolved_image(:));
    
    % Invert image?
    if invert_image
        convolved_image_DISPLAY = 1 - convolved_image;
    else
        convolved_image_DISPLAY = convolved_image;
    end
    
    % Display results.
    figure();
    subplot(1,3,1); image(yRange * 1e6,zRange * 1e6,image_target_DISPLAY);axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Image target');
    subplot(1,3,2); imagesc(yRange * 1e6,zRange * 1e6,squeeze(1 - lightSheetPsf).');axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Light-sheet PSF');
    colormap hot;
    subplot(1,3,3); imagesc(yRange * 1e6,zRange * 1e6,convolved_image_DISPLAY);axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Light-sheet image');
    
    % Display data.
    figure();
    subplot(1,3,1); image(yRange * 1e6,zRange * 1e6,image_target);axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Image target - data');
    subplot(1,3,2); imagesc(yRange * 1e6,zRange * 1e6,squeeze(lightSheetPsf).');axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Light-sheet PSF');
    colormap jet;
    subplot(1,3,3); imagesc(yRange * 1e6,zRange * 1e6,convolved_image);axis image;
    xlabel('x-axis [um]'); ylabel('z-axis [um]');
    title('Light-sheet image - data');
    
    % Switch image dimensions.
    convolved_image = permute(convolved_image,[3 2 1]);   % Colour channels become "xz planes".
    
    % Save "recorded data"
    recordedImageStack = convolved_image;
    save(strcat(output_folder,'\simulated_image.mat'),'recordedImageStack');
    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%
    %%% To run after deconvolution code:
    
    %%% load in restoredDataCube, yRange, and zRange.
%     load('G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\__Acquired Data\recording0_lambda532nm_alpha0_beta100.mat','restoredDataCube','yRange','zRange');
%     
%     % Re-order datacube and normalise.
%     deconvolved_image = permute(restoredDataCube,[3 2 1]);
%     deconvolved_image = deconvolved_image ./ max(deconvolved_image(:));
%     deconvolved_image = deconvolved_image .* (deconvolved_image > 0);
% 
%     % Additional scaling for Gaussian only
%     deconvolved_image = deconvolved_image *2.4;
%     deconvolved_image = deconvolved_image .* (deconvolved_image < 1) + (deconvolved_image >1);
% 
%     deconvolved_image_DISPLAY = 1 - deconvolved_image;
%     
%     figure();imagesc(yRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);axis image;
%     xlabel('x-axis [um]'); ylabel('z-axis [um]');
%     title('Light-sheet reconstruction');
%    
%     deconvolved_image_DISPLAY = deconvolved_image;
%     
%     figure();imagesc(yRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);axis image;
%     xlabel('x-axis [um]'); ylabel('z-axis [um]');
%     title('Light-sheet reconstruction');
% 
%     load('G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Simulation_Examples\__Acquired Data\recording0_lambda532nm_alpha7_beta100.mat','restoredDataCube','yRange','zRange');
%     
%     % Re-order datacube and normalise.
%     deconvolved_image = permute(restoredDataCube,[3 2 1]);
%     deconvolved_image = deconvolved_image ./ max(deconvolved_image(:));
%     deconvolved_image = deconvolved_image .* (deconvolved_image > 0);
% 
%     deconvolved_image_DISPLAY = 1 - deconvolved_image;
%     
%     figure();imagesc(yRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);axis image;
%     xlabel('x-axis [um]'); ylabel('z-axis [um]');
%     title('Light-sheet reconstruction');
% 
%     deconvolved_image_DISPLAY = deconvolved_image;
%     
%     figure();imagesc(yRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);axis image;
%     xlabel('x-axis [um]'); ylabel('z-axis [um]');
%     title('Light-sheet reconstruction');


end