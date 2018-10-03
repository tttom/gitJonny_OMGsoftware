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


function FastLightSheetProfileSimulation_3PE(zRange,xRange,lambda,NA,outputFolder,beamTypes)


    %default inputs
    if nargin < 1
%         zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
%         zRange = [-50:0.1:50] * 1e-6;  % transverse beam axis [metres]
        zRange = [-25:0.05:25] * 1e-6;  % transverse beam axis [metres]
    end
    if nargin < 2
%         xRange = [-10:0.5:175] * 1e-6; % propagation axis [metres]
%         xRange = [0:5:50] * 1e-6;
%         xRange = [-1:1:15] * 1e-6;
%         xRange = [-100:100:3000] * 1e-6;
        xRange = [-10:0.5:50] * 1e-6;
    end
    if nargin < 3
        lambda = 488e-9;    % [metres]
%         lambda = 333e-9;    % [metres]
    end
    if nargin < 4
        NA = 0.4;
%         NA = 0.23;
%         NA = 0.20;
    end
    if nargin < 5
%         outputFolder = 'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\NeurophotonicsBookChapter\LightSheetSimulations\2018-02-02';
        outputFolder = ...
            'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\3PELSM\beamAnimation';
    end
    if nargin < 6
%         beamTypes = {'Gaussian','Gaussian (confocal)','Bessel','Bessel (confocal)','Bessel (2PE)','Airy'};
%         beamTypes = {'Bessel (2PE)'};
%         beamTypes = {'Gaussian','Gaussian (2PE)','Gaussian (3PE)'};
%         beamTypes = {'Gaussian','Gaussian (2PE)','Gaussian (3PE)'...
%             ,'Bessel10','Bessel10 (2PE)','Bessel10 (3PE)'...
%             ,'Bessel5','Bessel5 (2PE)','Bessel5 (3PE)'...
%             ,'Bessel2','Bessel2 (2PE)','Bessel2 (3PE)'...
%             ,'Bessel1','Bessel1 (2PE)','Bessel1 (3PE)'...
%             };
%         beamTypes = {'Bessel5 (3PE)'};
%         beamTypes = {'Gaussian (3PE)'};
        beamTypes = {'Airy'};
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
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
                
            case 'Bessel10'
                beta = 0.1;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            case 'Bessel10 (2PE)'
                beta = 0.1;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 2; %2PE
            case 'Bessel10 (3PE)'
                beta = 0.1;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
            case 'Bessel5'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            case 'Bessel5 (2PE)'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 2; %2PE
            case 'Bessel5 (3PE)'
                beta = 0.05;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
            case 'Bessel2'
                beta = 0.02;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            case 'Bessel2 (2PE)'
                beta = 0.02;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 2; %2PE
            case 'Bessel2 (3PE)'
                beta = 0.02;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
            case 'Bessel1'
                beta = 0.01;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            case 'Bessel1 (2PE)'
                beta = 0.01;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 2; %2PE
            case 'Bessel1 (3PE)'
                beta = 0.01;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 3; %3PE
                
                
            case 'Airy'
                alpha = 5;
%                 alpha = 7;
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* exp(2* pi * 1i * alpha .* (U.^3 + V.^3));
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
            otherwise
                pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
                confocal_slit_width = 1; %[metres]  (no confocal detection)
                intensity_order = 1; %1PE
        end

%         [psf,~,intensity_orders] = calcVectorialPsf(yRange,zRange,xRange,lambda .* intensity_order...
%             ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
%             ,NA,ref_index,illumination_mag,illumination_tube_length,1);
        [psf,~,~] = calcVectorialPsf(yRange,zRange,xRange,lambda .* intensity_order...
            ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
            ,NA,ref_index,illumination_mag,illumination_tube_length,[],intensity_order);
        
        switch beamType
            case 'Gaussian (2PE)'
                psf = intensity_orders(:,:,:,intensity_order - 1);
            case 'Gaussian (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
                
            case 'Bessel10 (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel10 (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel5 (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel5 (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel2 (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel2 (3PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel1 (2PE)'
                psf = intensity_orders(:,:,:,intensity_order  - 1);
            case 'Bessel1 (3PE)'
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
            Fig_beamProfile = figure(1);
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

            
            
            thresholdedMTF = MTF_z >= 0.05;
            resolutionAcrossFOV = zeros([1,size(thresholdedMTF,2)]);
            for xIdx = 1:size(thresholdedMTF,2)
                firstZero = find(thresholdedMTF(:,xIdx) == 0,1,'first');
                thresholdedMTF(firstZero:end,xIdx) = 0;
                resolutionAcrossFOV(xIdx) = kzRange(max(firstZero - 1,1));
            end
            
            % plot thresholded MTF
%             subplot(3,1,3);
%             imagesc(xRange * 1e6,kzRange,thresholdedMTF);
%             title('MTF(x,f_z)');
%             xlabel('x [um]');ylabel('f_z');
%             ylim([0 1]);
%             drawnow;shg;
            subplot(3,1,3);
            contour(xRange * 1e6,kzRange,thresholdedMTF,1);
            title('MTF(x,f_z)');
            xlabel('x [um]');ylabel('f_z');
            ylim([0 1]);
            drawnow;shg;
            
            %find x = 0
            x0_pix = find(xRange == 0,1,'first');
            
            %%% resolution (MTF)
            bestResolution = max(resolutionAcrossFOV);
            
            %%% FOV by resolution (MTF)
            halfFOVpixel1 = find(resolutionAcrossFOV(x0_pix:end) <= bestResolution / sqrt(2),1,'first');
            halfFOV1 = xRange(halfFOVpixel1); % MTF drops to 1 / sqrt(2)
            halfFOVpixel2 = find(resolutionAcrossFOV(x0_pix:end) <= bestResolution / 2,1,'first');
            halfFOV2 = xRange(halfFOVpixel2); % MTF drops to 1 / 2
            
            %%% FOV by intensity
            intensityProfile = max(psf,[],1);
            halfFOVpixel3 = find(intensityProfile(x0_pix:end) <= 1 / sqrt(2),1,'first');
            halfFOV3 = xRange(halfFOVpixel3 + (x0_pix - 1)); % Peak intensity drops to 1 / 2
            halfFOVpixel4 = find(intensityProfile(x0_pix:end) <= 1 / 2,1,'first');
            halfFOV4 = xRange(halfFOVpixel4  + (x0_pix - 1)); % peak intensity drops to 1 / 4
            
            %%% resolution from PSF FWHM
            psfFWHM = zeros(size(intensityProfile));
            for xIdx = 1:length(xRange)
                start_coord = find(psf(:,xIdx) >= (intensityProfile(xIdx) / 2),1,'first');
                end_coord = find(psf(:,xIdx) >= (intensityProfile(xIdx) / 2),1,'last');
                psfFWHM(xIdx) = (zRange(end_coord) - zRange(start_coord)) * 1e6;
            end
            
            Fig_intensityProfile = figure(2);
            subplot(2,1,1);plot(xRange * 1e6,intensityProfile);
            title('Longitudinal intensity profile');
            xlabel('x-axis [um]');
            ylabel('Intensity [a.u.]');
            subplot(2,1,2);plot(xRange * 1e6,psfFWHM);
            title('FWHM of psf');
            xlabel('x-axis[um]');
            ylabel('z-axis FWHM [um]');
            ylim([0 10]);
            drawnow;shg;
            
            % save figures
            mkdir(outputFolder);
            figure(Fig_beamProfile)
            saveas(gcf,strcat(outputFolder,'\',beamType,'_beamProfile.fig'));
            figure(Fig_intensityProfile)
            saveas(gcf,strcat(outputFolder,'\',beamType,'_intensityProfile.fig'));
            save(strcat(outputFolder,'\',beamType,'.mat')...
                ,'bestResolution','halfFOV1','halfFOV2','halfFOV3','halfFOV4','psfFWHM','intensityProfile');

    end
        
end


% 
% 
%         MTF_z = get(get(gca,'Children'),'CData');
%         kzRange = get(get(gca,'Children'),'YData');
%         xRange = get(get(gca,'Children'),'XData');
%         %find x = 0
%         x0_pix = find(xRange == 0,1,'first');
%         thresholdedMTF = MTF_z >= 0.05;
%         resolutionAcrossFOV = zeros([1,size(thresholdedMTF,2)]);
%         for xIdx = 1:size(thresholdedMTF,2)
%             firstZero = find(thresholdedMTF(:,xIdx) == 0,1,'first');
%             thresholdedMTF(firstZero:end,xIdx) = 0;
%             resolutionAcrossFOV(xIdx) = kzRange(max(firstZero - 1,1));
%         end
%         %%% resolution
%         bestResolution = max(resolutionAcrossFOV);
%         %%% FOV by resolution
%         halfFOVpixel1 = find(resolutionAcrossFOV >= bestResolution / sqrt(2),1,'last');
%         halfFOV1 = xRange(halfFOVpixel1);
%         halfFOVpixel2 = find(resolutionAcrossFOV >= bestResolution / 2,1,'last');
%         halfFOV2 = xRange(halfFOVpixel2);
%         
%         
%         
%         
%         psf = get(get(gca,'Children'),'CData');
%         zRange = get(get(gca,'Children'),'YData');
%         xRange = get(get(gca,'Children'),'XData');
%         %%% FOV by intensity
%         intensityProfile = max(psf,[],1);
%         halfFOVpixel3 = find(intensityProfile >= 1 / 2,1,'last');
%         halfFOV3 = xRange(halfFOVpixel3);
%         halfFOVpixel4 = find(intensityProfile >= 1 / 4,1,'last');
%         halfFOV4 = xRange(halfFOVpixel4);
%         
%         
%         
%         saveas(gcf,'Bessel1 (3PE).fig');
%         save('Bessel1 (3PE).mat','bestResolution','halfFOV1','halfFOV2','halfFOV3','halfFOV4');
%         