%%% name:           FastLightSheetProfileSimulation_attenuationCompensation
%%% author:         Jonathan Nylk
%%% date created:   19/07/2016
%%% description:    This function simulates attenuation compensated
%%%                 light-sheets and plots their relative pupil functions,
%%%                 light-sheet cross-sectional profiles, MTF as a
%%%                 function of propagation distance, and the maximum
%%%                 intensity as a function of propagation.
%%%
%%% updates (latest first):
%%%                 
%%%
%%%
%%% END %%%


function FastLightSheetProfileSimulation_attenuationCompensation(zRange,xRange,lambda,NA,alpha,compensation_type...
    ,sigma_exp,sigma_lin,sigma_peak,x_peak,C_abs,outputFolder,Flag_imageSimulation,deconvolution_lightSheet,amp_peak)

    %default inputs
    if nargin < 1
%         zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
        zRange = [-50:0.1:50] * 1e-6;  % transverse beam axis [metres]
    end
    if nargin < 2
%         xRange = [-100:1:100] * 1e-6; % propagation axis [metres]
        xRange = [-1:1:1] * 1e-6;
    end
    if nargin < 3
        lambda = 532e-9;    % [metres]
    end
    if nargin < 4
        NA = 0.42;
%         NA = 0.4;
    end
    if nargin < 5
        alpha = 7;
%         alpha = 0;
    end
    if nargin < 6
        compensation_type = 'Exponential';
%         compensation_type = 'Linear';
%         compensation_type = 'Peak';
    end
    if nargin < 7
        sigma_exp = 0.54;
%         sigma_exp = 0.27;
%         sigma_exp = 0;
%         sigma_exp = 0.46;
    end
    if nargin < 8
        sigma_lin = 0.54;
    end
    if nargin < 9
        sigma_peak = 0.1;
    end
    if nargin < 10
%         x_peak = 0 * 1e-6; % [metres]
        x_peak = 0.75;
    end
    if nargin < 11
        C_abs = 64.95 * 100;   % [metres^-1]
%         C_abs = 0;   % [metres^-1]
    end
    if nargin < 12
        outputFolder = 'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\CompensatedLightsheet\Theory and Simulations\2017-04-19_tweakedOTFcalc';
    end
    if nargin < 13
        Flag_imageSimulation = 1;   % performs image convolution and deconvolution
%         Flag_imageSimulation = 0;   % performs image convolution and deconvolution
    end
    if nargin < 14
%         deconvolution_lightSheet = 'Normal';      % Use calculated light-sheet with no attenuation included. 
        deconvolution_lightSheet = 'APriori';     % Use calculated light-sheet with attenuation included.
    end
    if nargin < 15
        amp_peak = 0.5;
    end
    
    % If imaging simulation is to be performed, override certain inputs to
    % comply with image size.
    if Flag_imageSimulation
        zRange = [-62:0.2:62.2] * 1e-6;  % transverse beam axis [metres]
        xRange = [-200:0.2:199.8] * 1e-6; % propagation axis [metres]
    end
    
    
    % other variable definitions
    ref_index = 1.33;
    illumination_mag = 40;
    illumination_tube_length = 0.2; % [metres]
    
    % define figure handles
    fig_pupil_amplitude = 101;
    fig_lightSheet_profile = 102;
    fig_intensity_maxima = 103;
    fig_MTF = 104;
    fig_MTF_5percent = 105;
    switch Flag_imageSimulation
        case 1
            fig_image_sims = 106;
            fig_image_sims_2 = 107;
        otherwise
            % Otherwise do nothing
    end
    
    % cylindrical simulation so no y-axis needed
    yRange = [0] * 1e-6;    % transverse beam axis (in plane of light-sheet) [metres]
%     yRange = zRange;

    % soft-aperture implemented, therefore extend hard aperture limit of NA
    NA_pupil = NA / 0.7;

    % define spatial frequency coords
    kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA * lambda;
    kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]
    kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) * 1 / 2 / xRange(end) / 2 / NA * lambda;
    kxRange = kxRange(ceil(length(kxRange) / 2):end);     % normalised [0:1]

    % set pupil limits
        VPupilSize = 1;       % z-axis
        UPupilSize = 0.01;    % y-axis
%         UPupilSize = 1;     % y-axis

    % set hard-aperture at pupil limit
    pupilAmplitudeMaskU = @(U,V) 1 * (U >= -UPupilSize & U <= UPupilSize);
    pupilAmplitudeMaskV = @(U,V) 1 * (V >= -VPupilSize & V <= VPupilSize);

    % set super-Gaussian apodization
    sigma_superGaussian = 8;    %super-Gaussian order (must be even!)
%     sigma_superGaussian = 100;    % TO REMOVE EFFECT OF SUPER-GUASSIAN APODIZATION
    superGaussian = @(U,V) exp(-(sqrt(2) * U).^ sigma_superGaussian) .* exp(-(sqrt(2) * V).^ sigma_superGaussian);
%     superGaussian = @(U,V) exp(-(sqrt(2) * sqrt(U.^2 + V.^2)).^ sigma_superGaussian); % CIRCULAR PUPIL
    
    pupilPhaseModulations = cell([1,2]);
    
    % compensated beam definition
    switch compensation_type
        case 'Exponential'
            pupilPhaseModulations{1} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * V.^3)...
                .* exp(sigma_exp / (NA / NA_pupil) * (V - 1)) ...
                .* superGaussian(U,V);    % exponential compenation factor (1+1D Airy)
%             pupilPhaseModulations{1} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * (U.^3 + V.^3))...
%                 .* exp(sigma_exp / (NA / NA_pupil) * (V - 1)) .* superGaussian(U,V);    % exponential compenation factor (2+1D Airy)
            fprintf('Pupil function: %s.\n Parameters:\n Alpha = %f.\n Sigma = %f.\n',func2str(pupilPhaseModulations{1}),alpha,sigma_exp);
        case 'Linear'
            pupilPhaseModulations{1} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * V.^3)...
                .* sigma_exp / (NA / NA_pupil) .* (V + 1) .* superGaussian(U,V);   % linear compensation factor
            fprintf('Pupil function: %s.\n Parameters:\n Alpha = %f.\n Sigma = %f.\n',func2str(pupilPhaseModulations{1}),alpha,sigma_lin);
        case 'Peak'
           pupilPhaseModulations{1} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * V.^3)...
                .* ((1 - amp_peak) + amp_peak .* exp(-((V - (x_peak * (NA / NA_pupil))) / sqrt(2) / (sigma_peak * (NA / NA_pupil))).^2)) .* superGaussian(U,V);    % no compensation, intensity peak
           fprintf('Pupil function: %s.\n Parameters:\n Alpha = %f.\n Sigma (width) = %f.\n Peak position = %f.\n Peak height = %f.\n',func2str(pupilPhaseModulations{1}),alpha,sigma_peak,x_peak,amp_peak);
        otherwise
           fprintf('Case: "%s" not recognised/implemented yet, using case: "Exponential" for compensation instead. \n',compensation_type)
           pupilPhaseModulations{1} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * V.^3)...
                .* exp(sigma_exp / (NA / NA_pupil) * (V - 1)) .* superGaussian(U,V);    % exponential compenation factor
           fprintf('Pupil function: %s.\n Parameters:\n Alpha = %f.\n Sigma = %f.\n',func2str(pupilPhaseModulations{1}),alpha,sigma_exp);
    end
    % standard beam definition
    pupilPhaseModulations{2} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * V.^3)...
        .* superGaussian(U,V);  % no compensation
%     pupilPhaseModulations{2} = @(U,V) exp(2i * pi * alpha / (NA / NA_pupil).^3 * (U.^3 + V.^3))...
%         .* superGaussian(U,V);  % no compensation    
    
    % allocate memory for PSFs
    PSF = zeros(length(zRange),length(xRange),length(pupilPhaseModulations));
    
    % evaluate compensated pupil function to determine normalisation constants
    pupilRangeU = [-UPupilSize:0.01:UPupilSize];
    pupilRangeV = [-VPupilSize:0.01:VPupilSize];
    if strcmp(compensation_type,'Peak')
        pupilMaxVal = 1;
        pupilNormVal = 1;
    else
        pupilPhaseModulation = pupilPhaseModulations{1};
        [maxVal,maxPos] = max(abs(pupilPhaseModulation(0,pupilRangeV)));
        pupilMaxVal = maxVal;
        pupilNormVal = abs(pupilPhaseModulation(0,-pupilRangeV(maxPos)));
    end
    
    for Idx = 1:length(pupilPhaseModulations)
        % scale pupil functions with normalisation constants
        pupilPhaseModulation = pupilPhaseModulations{Idx};
        pupilFunctor = @(U,V) pupilAmplitudeMaskU(U,V) .* pupilAmplitudeMaskV(U,V)...
            .* pupilPhaseModulation(U,V);
        if Idx == 1
            pupil_scaling_factor = 1 / pupilMaxVal;
        else
            pupil_scaling_factor = pupilNormVal / pupilMaxVal;
        end
        pupilFunctor = @(U,V) pupilFunctor(U,V) * pupil_scaling_factor;
        fprintf('Pupil scaling factor = %f.\n',pupil_scaling_factor);
        
        % sum amplitude over pupil
        total_pupil_amplitude = sum(abs(pupilFunctor(0,pupilRangeV)));
        fprintf('Total pupil amplitude = %f.\n',total_pupil_amplitude);
        total_pupil_squared_amplitude = sum(abs(pupilFunctor(0,pupilRangeV)).^2);
        fprintf('Total pupil squared amplitude = %f.\n',total_pupil_squared_amplitude);
        total_pupil_amplitude_squared = sum(abs(pupilFunctor(0,pupilRangeV))).^2;
        fprintf('Total pupil amplitude squared = %f.\n',total_pupil_amplitude_squared);
        
        % plot pupil amplitude
        figure(fig_pupil_amplitude);
        subplot(1,2,Idx);
        plot(pupilRangeV,abs(pupilFunctor(0,pupilRangeV)));
        xlim([-1 1]); ylim([0 1]); axis square;
        xlabel('u-axis'); ylabel('Amplitude [a.u.]');
        title(sprintf('Pupil amplitude. \n\n Total pupil amplitude = %f. \n\n Pupil scaling factor = %f.'...
            ,total_pupil_amplitude,pupil_scaling_factor));
        drawnow;shg;
        
        % plot 2D pupil amplitude and phase
        [pupilRangeU_2D,pupilRangeV_2D] = meshgrid(pupilRangeU,pupilRangeV);
        figure(51);
        subplot(1,2,Idx);
        imagesc(pupilRangeV,pupilRangeU,abs(pupilFunctor(pupilRangeU_2D,pupilRangeV_2D)).');
        axis image;
        xlim([-1 1]); ylim([-1 1]);
        xlabel('u-axis'); ylabel('v-axis');
        title('amplitude');
        figure(52);
        subplot(1,2,Idx);
        imagesc(pupilRangeV,pupilRangeU,angle(pupilFunctor(pupilRangeU_2D,pupilRangeV_2D)).');
        axis image;
        xlim([-1 1]); ylim([-1 1]);
        xlabel('u-axis'); ylabel('v-axis');
        title('phase');
        drawnow;shg;
        
        
        % calculate PSF
        psf_timer = tic;
        psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
            ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
            ,NA_pupil,ref_index,illumination_mag,illumination_tube_length,[]);
        psf_time_elapsed = toc(psf_timer);
        fprintf('Time taken to calculate light-sheet PSF: %f seconds.\n',psf_time_elapsed);
        
        % if a 2D pupil function, then 3D beam profile.
        % project along y-axis to determine light-sheet profile.
        if max(yRange~=0)
            psf = sum(psf,1);
        end
        
        % re-orient PSF
        psf = squeeze(permute(psf,[1,3,2])).';
        
        % take finite depth-of-field of detection lens into account
        [~,zCoords] = meshgrid(xRange,zRange);
        DOF_half_width = 16e-6;  % [metres]
        detection_DOF = exp(-1 * zCoords.^2 / 2 / DOF_half_width^2);
        psf = psf .* detection_DOF;
        
% % % % % % % %         % TOTFC
% % % % % % % %         % normalise
% % % % % % % %         psf = psf ./ max(psf(:));
        
        % take a copy of the PSF before applying absorption
        psf_before_absorption = psf;
        
        % simulate "perfect" absorption
        absorption_decay = repmat(exp(-C_abs * xRange),[length(zRange) 1]);
        psf = psf .* absorption_decay;
        
% % % % % % % %         % TOTFC
% % % % % % % %         % re-normalise
% % % % % % % %         psf = psf ./ max(psf(:));
        
        % write psf
        PSF(:,:,Idx) = psf;
        
        % plot attenuated light-sheet
        figure(fig_lightSheet_profile);
        subplot(2,1,Idx);
        imagesc(xRange * 1e6,zRange * 1e6,psf); axis image;
        title('Light-sheet profile');
        xlabel('x-axis [um]'); ylabel('z-axis [um]');
        drawnow;shg;

        % determine transverse intensity maxima
        [C_max,I_max]=max(psf,[],1);
        
        % plot position and intensity of maxima
        figure(fig_intensity_maxima);
        subplot(2,2,Idx)
        plot(xRange * 1e6,C_max); axis square;
        title('Light-sheet peak intensity value');
        xlabel('x-axis [um]'); ylabel('Intensity [a.u.]');
        xlim([xRange(1) xRange(end)] * 1e6);
        subplot(2,2,Idx + 2);
        plot(xRange * 1e6,zRange(I_max) * 1e6); axis square;
        title('Light-sheet peak intensity position');
        xlabel('x-axis [um]'); ylabel('z-axis [um]');
        xlim([xRange(1) xRange(end)] * 1e6);
        drawnow;shg;
        
        % determine MTF_z as a function of x-axis coordinate
        MTF_z = abs(fftshift(fft(psf,[],1),1));
        % normalisation
% % % % % %             % TOTFC
                        xRange0=find(xRange>=0,1); % approx. closest value to x=0
                        MTF_z = MTF_z ./ max(MTF_z(:,xRange0));
% % % % % %              MTF_z = MTF_z ./ max(MTF_z(:));
        MTF_z = MTF_z(ceil(size(MTF_z,1) / 2):end,:); % range, kz = [0:1]
        
% % % % % % % %         % TOTFC
% % % % % % % %         if Idx == 1
% % % % % % % %             MTF_z = MTF_z / pupilMaxVal / pupilNormVal;
% % % % % % % %             fprintf('MTF scaling factor = %f.\n',1 / pupilMaxVal / pupilNormVal);
% % % % % % % %         end
        
        % plot MTF
        figure(fig_MTF);
        subplot(2,1,Idx);
        imagesc(xRange * 1e6,kzRange,MTF_z);
        if Idx == 1
            title(sprintf('MTF(x,kz). \n\n MTF scaling factor = %f.\n',1 / pupilMaxVal / pupilNormVal));
        else
            title('MTF(x,kz)');
        end
        xlabel('x [um]');ylabel('kz');
        ylim([0 1]);
        drawnow;shg;
        
        % plot thresholded MTF
        figure(fig_MTF_5percent);
        subplot(2,1,Idx);
        imagesc(xRange * 1e6,kzRange,MTF_z >= 0.05);
        if Idx == 1
            title(sprintf('MTF(x,kz). \n\n MTF scaling factor = %f.\n',1 / pupilMaxVal / pupilNormVal));
        else
            title('MTF(x,kz)');
        end
        xlabel('x [um]');ylabel('kz');
        colormap gray;
        ylim([0 1]);
        drawnow;shg;
        
        
        % Simulate imaging of target
        switch Flag_imageSimulation
            case 1
                % Load image_target
                image_target_base =...
                    single(imread('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\CompensatedLightsheet\Theory and Simulations\03-foundation-colour.bmp'));
%                 image_target = imresize(image_target,2);
                image_target_base = repmat(image_target_base,[1 2 1]); 
                image_target = ones(size(repmat(image_target_base,[2 1 1]))) * max(image_target_base(:));
                image_target(floor(size(image_target_base,1) / 2):floor(size(image_target_base,1) / 2) -1 + size(image_target_base,1),:,:) = image_target_base;
                % Normalise
                if min(image_target(:)) < 0
                    image_target = image_target .* (image_target>0);    % Ensures all +ve values
                end
                image_target = image_target / max(image_target(:));

                % Image target to display
                image_target_DISPLAY = image_target;

                % invert image target
                image_target = 1 - image_target;


                % Convolve image target with light-sheet PSF
                convolved_image = zeros([size(image_target,1),size(image_target,2),size(image_target,3)],'single');
                for m = 1:size(image_target,3)
                    for n = 1:size(image_target,2)
                        convolved_image(:,n,m) = conv(squeeze(image_target(:,n,m)),squeeze(psf(:,n)),'same');
                    end
                end

                % Add Gaussian noise
                noise_magnitude = 0.01;
%                 noise_magnitude = 0; % no noise
                image_noise = randn(size(convolved_image)) .* max(convolved_image(:)) * noise_magnitude;
                convolved_image = convolved_image + image_noise;

                % Re-normalisation of image
                if min(convolved_image(:)) < 0
                    convolved_image = convolved_image - min(convolved_image(:));    % Ensures all +ve values
                end
                convolved_image = convolved_image / max(convolved_image(:));

                % Convolved image to display
                convolved_image_DISPLAY = 1 - convolved_image;

                % Deconvolve image with light-sheet PSF
                
                % Deconvolve using light-sheet profile without including absorption
                switch deconvolution_lightSheet
                    case 'Normal'
                        psf = psf_before_absorption;
                    otherwise
                end
                
                % Set attenuation of light-sheet for deconvolution purposes
                % incorrectly
                
                C_abs_deconv = C_abs * 1;
                psf = psf_before_absorption;
                absorption_decay_deconv = repmat(exp(-C_abs_deconv * xRange),[length(zRange) 1]);
                psf = psf .* absorption_decay_deconv;
% % % % %                 % TOTFC
% % % % %                 % re-normalise
% % % % %                 psf = psf ./ max(psf(:));
                
                
                % Zero-pad and centre PSF
                padded_psf = zeros(size(psf,1) + 2 * ceil(size(psf,1) / 2),size(psf,2));
                padded_psf(ceil(size(psf,1) / 2) + [1:size(psf,1)],:) = psf;
                
                % Model z-OTF as inverse function
                ZOtf = ([1:size(padded_psf,1)] - floor(size(padded_psf,1) / 2) - 1) / ((zRange(2)-zRange(1)) * size(padded_psf,1));
                excitationOpticalCutOffSpFreq = 2 * NA / lambda; % Independent of refractive index of medium
                excitationNoiseToSignalRatio = ZOtf / excitationOpticalCutOffSpFreq / 5;
                excitationNoiseToSignalRatio =excitationNoiseToSignalRatio.';
                
                % Calculate Deconvolution filter
                lightSheetOtf = fftshift(fft(ifftshift(padded_psf,1),[],1),1);
% % % % % %                 % TOTFC
                            xRange0=find(xRange>=0,1); % approx. closest value to x=0
                            lightSheetOtf = lightSheetOtf ./ max(abs(lightSheetOtf(:,xRange0)));
% % % % % % %                 lightSheetOtf = lightSheetOtf ./ max(abs(lightSheetOtf(:))); % Maintain the mean brightness
                
% % % % % %                 % TOTFC
% % % % % %                 if Idx == 1
% % % % % %                     lightSheetOtf = lightSheetOtf / pupilMaxVal / pupilNormVal;
% % % % % %                     fprintf('OTF scaling factor = %f.\n',1 / pupilMaxVal / pupilNormVal);
% % % % % %                 end
                
                lightSheetDeconvFilter = conj(lightSheetOtf) ./ (abs(lightSheetOtf).^2 ...
                    + repmat(excitationNoiseToSignalRatio.^2,[1 size(lightSheetOtf,2)]));
                
                % Pad and extend edges of image
                padded_convolved_image = convolved_image([1:end, end * ones(1,floor(end / 2))...
                        ,ones(1,floor((end + 1) / 2))],:,:);
                
                % Deconvolve
                imageSliceFft = fft(padded_convolved_image,[],1);
                deconvImage = ifft(imageSliceFft .* repmat(ifftshift(lightSheetDeconvFilter,1),[1 1 size(convolved_image,3)]),[],1,'symmetric');
                deconvolved_image = deconvImage(1:end / 2,:,:); % drop wrapped kernel padding;

                % Re-normlaise
                deconvolved_image = deconvolved_image ./ max(deconvolved_image(:));
                deconvolved_image = deconvolved_image .* (deconvolved_image > 0);
                
                % Re-normlaise within specified xRange
                xRange_min = -100e-6;   % [metres]
                xRange_max = 100e-6;    % [metres]
                [~,xLim_lower] = max(xRange >= xRange_min);
                [~,xLim_upper] = max(xRange >= xRange_max);
                image_min_level = min(min(min(deconvolved_image(:,xLim_lower:xLim_upper,:))));
                image_max_level = max(max(max(deconvolved_image(:,xLim_lower:xLim_upper,:))));
                deconvolved_image = deconvolved_image - image_min_level;
                deconvolved_image = deconvolved_image / image_max_level;
                deconvolved_image = deconvolved_image .* (deconvolved_image > 0);
                deconvolved_image = deconvolved_image .* (deconvolved_image < 1) + (deconvolved_image >= 1);

                % Deconvolved image to display
                deconvolved_image_DISPLAY = 1 - deconvolved_image;

                figure(fig_image_sims);
                subplot(3,2,Idx);
                imagesc(xRange * 1e6,zRange * 1e6,image_target_DISPLAY);
                axis image;
                xlabel('x-axis [um]');
                ylabel('z-axis [um]');
                title('image target');
                subplot(3,2,2+Idx);
                imagesc(xRange * 1e6,zRange * 1e6,convolved_image_DISPLAY);
                axis image;
                xlabel('x-axis [um]');
                ylabel('z-axis [um]');
                title('convolved image');
                subplot(3,2,4+Idx);
                imagesc(xRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);
                axis image;
                xlabel('x-axis [um]');
                ylabel('z-axis [um]');
                title('deconvolved image');
                drawnow;shg;
                
                figure(fig_image_sims_2);
                subplot(2,1,Idx);
                imagesc(xRange * 1e6,zRange * 1e6,deconvolved_image_DISPLAY);
                axis image;
                xlabel('x-axis [um]');
                ylabel('z-axis [um]');
                title('deconvolved image');
                drawnow;shg;
                
            otherwise
                % Otherwise do nothing
        end
        
    end
    
    % save figures
    switch compensation_type
        case 'Exponential'
            fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_ExponentialCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_exp);
        case 'Linear'
            fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_LinearCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_lin);
        case 'Peak'
 
           fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_Peak_sigma_%f_centre_%f_amp_%f',outputFolder,alpha,C_abs,sigma_peak,x_peak,amp_peak);
        otherwise
           fprintf('Case: "%s" not recognised/implemented yet, using case: "Exponential" for compensation instead. \n',compensation_type);
           fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_ExponentialCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_exp);
    end
    mkdir(fullOutputFolder);
    saveas(fig_pupil_amplitude,sprintf('%s\\pupilAmplitude.fig',fullOutputFolder));
    saveas(fig_pupil_amplitude,sprintf('%s\\pupilAmplitude.png',fullOutputFolder));
    saveas(fig_lightSheet_profile,sprintf('%s\\lightSheetProfile.fig',fullOutputFolder));
    saveas(fig_lightSheet_profile,sprintf('%s\\lightSheetProfile.png',fullOutputFolder));
    saveas(fig_intensity_maxima,sprintf('%s\\intensityMaxima.fig',fullOutputFolder));
    saveas(fig_intensity_maxima,sprintf('%s\\intensityMaxima.png',fullOutputFolder));
    saveas(fig_MTF,sprintf('%s\\MTF.fig',fullOutputFolder));
    saveas(fig_MTF,sprintf('%s\\MTF.png',fullOutputFolder));
    saveas(fig_MTF_5percent,sprintf('%s\\MTF_5.fig',fullOutputFolder));
    saveas(fig_MTF_5percent,sprintf('%s\\MTF_5.png',fullOutputFolder));
     switch Flag_imageSimulation
            case 1
                saveas(fig_image_sims,sprintf('%s\\imaging_simulation.fig',fullOutputFolder));
                saveas(fig_image_sims,sprintf('%s\\imaging_simulation.png',fullOutputFolder));
                saveas(fig_image_sims_2,sprintf('%s\\deconvolved_image.fig',fullOutputFolder));
                saveas(fig_image_sims_2,sprintf('%s\\deconvolved_image.png',fullOutputFolder));
         otherwise
                % Otherwise do nothing
     end
    
%     close all
        
end
