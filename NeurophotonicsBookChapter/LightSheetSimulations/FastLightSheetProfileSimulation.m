%%% name:           FastLightSheetProfileSimulation
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


function FastLightSheetProfileSimulation(zRange,xRange,lambda,NA,outputFolder,beamTypes)

    %default inputs
    if nargin < 1
        zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
%         zRange = [-50:0.1:50] * 1e-6;  % transverse beam axis [metres]
    end
    if nargin < 2
%         xRange = [-10:0.5:175] * 1e-6; % propagation axis [metres]
%         xRange = [-30:2:30] * 1e-6;
        xRange = [-100:100:3000] * 1e-6;
    end
    if nargin < 3
%         lambda = 488e-9;    % [metres]
        lambda = 333e-9;    % [metres]
    end
    if nargin < 4
%         NA = 0.4;
        NA = 0.23;
    end
    if nargin < 5
%         outputFolder = 'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\NeurophotonicsBookChapter\LightSheetSimulations\2018-02-02';
        outputFolder = 'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\NeurophotonicsBookChapter\LightSheetSimulations\3PhotonExcitationSims';
    end
    if nargin < 6
%         beamTypes = {'Gaussian','Gaussian (confocal)','Bessel','Bessel (confocal)','Bessel (2PE)','Airy'};
%         beamTypes = {'Bessel (2PE)'};
%         beamTypes = {'Gaussian','Gaussian (2PE)','Gaussian (3PE)','Bessel','Bessel (2PE)','Bessel (3PE)'};
        beamTypes = {'Bessel (3PE)'};
    end
    
    
    % other variable definitions
    ref_index = 1.33;
    illumination_mag = 40;
    illumination_tube_length = 0.2; % [metres]
    
    yRange = zRange;
    % cylindrical simulation so no y-axis needed
%     yRange = [0] * 1e-6;    % transverse beam axis (in plane of light-sheet) [metres]


    % define spatial frequency coords
    kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA * lambda;
    kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]
    kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) * 1 / 2 / xRange(end) / 2 / NA * lambda;
    kxRange = kxRange(ceil(length(kxRange) / 2):end);     % normalised [0:1]

    for bIdx = 1:length(beamTypes)
        beamType = beamTypes{bIdx};

        switch beamType
            case 'Gaussian'
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = 1; %[metres] (no confocal detection)
                intensity_order = 1; %1PE
            case 'Gaussian (2PE)'
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = 1; %[metres] (no confocal detection)
                intensity_order = 2; %2PE
            case 'Gaussian (3PE)'
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = 1; %[metres] (no confocal detection)
                intensity_order = 3; %3PE
            case 'Gaussian (confocal)'
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = lambda / 2 / NA; %[metres] (half-width of Gaussian beam)
                intensity_order = 1; %1PE
            case 'Bessel'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            case 'Bessel (confocal)'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 2.405 * lambda / 2 / pi / NA / (1 - beta / 2); %[metres] (half-widht of Bessel core)
                intensity_order = 1; %1PE
            case 'Bessel (2PE)'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 2; %2PE
            case 'Bessel (3PE)'
                beta = 0.01;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
            case 'Airy'
                alpha = 7;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* exp(2* pi * 1i * alpha .* (U.^3 + V.^3));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            otherwise
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
        end

        [psf,~,intensity_orders] = calcVectorialPsf(yRange,zRange,xRange,lambda .* intensity_order...
            ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
            ,NA,ref_index,illumination_mag,illumination_tube_length,1);
        
        switch beamType
            case 'Gaussian (2PE)'
                psf = intensity_orders(:,:,:,intensity_order - 1);
            case 'Gaussian (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            otherwise
                clear intensity_orders
        end

        psf = psf / max(psf(:)); % normalise


            % re-orient PSF
            psf = squeeze(permute(psf,[1,3,2])).';

            % take finite depth-of-field of detection lens into account
            [~,zCoords] = meshgrid(xRange,zRange);
            DOF_half_width = 2 * lambda * ref_index / NA / NA;  % [metres]
            detection_DOF = exp(-1 * zCoords.^2 / 2 / DOF_half_width^2);
            psf = psf .* detection_DOF;
            % apply confocal slit
            confocal_slit = exp(-1 * zCoords.^2 / 2 / confocal_slit_width^2);
            psf = psf .* confocal_slit;

            psf = psf ./ max(psf(:));


            % plot light-sheet
            figure();
            subplot(3,1,1);
            imagesc(xRange * 1e6,zRange * 1e6,psf); axis image;
            title(beamType);
            ylim([-15 15]);
            xlabel('x-axis [um]'); ylabel('z-axis [um]');
            colormap(hot(2^8));
            drawnow;shg;

            % determine MTF_z as a function of x-axis coordinate
            MTF_z = abs(fftshift(fft(psf,[],1),1));
            % normalisation
            MTF_z = MTF_z ./ max(MTF_z(:));
            MTF_z = MTF_z(ceil(size(MTF_z,1) / 2):end,:); % range, kz = [0:1]


            % plot MTF
            subplot(3,1,2);
            imagesc(xRange * 1e6,kzRange,MTF_z);
            title('MTF(x,f_z)');
            xlabel('x [um]');ylabel('f_z');
            ylim([0 1]);
            drawnow;shg;

            % plot thresholded MTF
            subplot(3,1,3);
            imagesc(xRange * 1e6,kzRange,MTF_z >= 0.05);
            title('MTF(x,f_z)');
            xlabel('x [um]');ylabel('f_z');
            ylim([0 1]);
            drawnow;shg;


    end
    
%     % save figures
%     switch compensation_type
%         case 'Exponential'
%             fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_ExponentialCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_exp);
%         case 'Linear'
%             fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_LinearCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_lin);
%         case 'Peak'
%  
%            fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_Peak_sigma_%f_centre_%f_amp_%f',outputFolder,alpha,C_abs,sigma_peak,x_peak,amp_peak);
%         otherwise
%            fprintf('Case: "%s" not recognised/implemented yet, using case: "Exponential" for compensation instead. \n',compensation_type);
%            fullOutputFolder = sprintf('%s\\Alpha_%f_Abs_%f_ExponentialCompensation_sigma_%f',outputFolder,alpha,C_abs,sigma_exp);
%     end
%     mkdir(fullOutputFolder);
%     saveas(fig_pupil_amplitude,sprintf('%s\\pupilAmplitude.fig',fullOutputFolder));
%     saveas(fig_pupil_amplitude,sprintf('%s\\pupilAmplitude.png',fullOutputFolder));
%     saveas(fig_lightSheet_profile,sprintf('%s\\lightSheetProfile.fig',fullOutputFolder));
%     saveas(fig_lightSheet_profile,sprintf('%s\\lightSheetProfile.png',fullOutputFolder));
%     saveas(fig_intensity_maxima,sprintf('%s\\intensityMaxima.fig',fullOutputFolder));
%     saveas(fig_intensity_maxima,sprintf('%s\\intensityMaxima.png',fullOutputFolder));
%     saveas(fig_MTF,sprintf('%s\\MTF.fig',fullOutputFolder));
%     saveas(fig_MTF,sprintf('%s\\MTF.png',fullOutputFolder));
%     saveas(fig_MTF_5percent,sprintf('%s\\MTF_5.fig',fullOutputFolder));
%     saveas(fig_MTF_5percent,sprintf('%s\\MTF_5.png',fullOutputFolder));
%      switch Flag_imageSimulation
%             case 1
%                 saveas(fig_image_sims,sprintf('%s\\imaging_simulation.fig',fullOutputFolder));
%                 saveas(fig_image_sims,sprintf('%s\\imaging_simulation.png',fullOutputFolder));
%                 saveas(fig_image_sims_2,sprintf('%s\\deconvolved_image.fig',fullOutputFolder));
%                 saveas(fig_image_sims_2,sprintf('%s\\deconvolved_image.png',fullOutputFolder));
%          otherwise
%                 % Otherwise do nothing
%      end
%     
% %     close all
        
end
