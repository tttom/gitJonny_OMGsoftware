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


function [psf,psfDOF,MTF_z] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction)

    %default inputs
    if nargin < 1
        zRange = [-25:0.25:25] * 1e-6;  % transverse beam axis [metres]
    end
    if nargin < 2
        xRange = [-100:1:100] * 1e-6;
    end
    if nargin < 3
        lambda = 532e-9;    % [metres]
    end
    if nargin < 4
        NA_ill = 0.4;
    end
    if nargin < 5
        NA_det = 0.4;
    end
    if nargin < 6
        outputFolder = [];
    end
    if nargin < 7
        inputPupilFunction = @(U,V) exp(2* pi * 1i * 7 .* (V.^3));
    end
    
    % other variable definitions
    ref_index = 1.33;
    illumination_mag = 40;
    illumination_tube_length = 0.2; % [metres]
    intensity_order = 1; %1-photon excitation
    confocal_slit_width = 1; %[metres]  (no confocal detection)
    
    yRange = zRange;
    % cylindrical simulation so no y-axis needed
%     yRange = [0] * 1e-6;    % transverse beam axis (in plane of light-sheet) [metres]


    % define spatial frequency coords
    kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA_ill * lambda;
    kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]
%     kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) * 1 / 2 / xRange(end) / 2 / NA * lambda;
%     kxRange = kxRange(ceil(length(kxRange) / 2):end);     % normalised [0:1]  (not actually used)

    % sampling check
    if max(kzRange) <= 1
        disp('Transverse dimension (zRange) step size is too small to accurately determine MTF at this NA, reduce step size.');
    end

    pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* inputPupilFunction(U,V);


    [psf,~,intensity_order_psfs] = calcVectorialPsf(yRange,zRange,xRange,lambda .* intensity_order...
        ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
        ,NA_ill,ref_index,illumination_mag,illumination_tube_length,1);

    if intensity_order ~= 1
        psf = intensity_order_psfs;
    end

    psf = psf / max(psf(:)); % normalise


    % re-orient PSF
    psf = squeeze(permute(psf,[1,3,2])).';
    
    %psf_before applying DOF
    original_psf = psf;
    
    % take finite depth-of-field of detection lens into account
    [~,zCoords] = meshgrid(xRange,zRange);
    DOF_half_width = 2 * lambda * ref_index / NA_det / NA_det;  % [metres]
    detection_DOF = exp(-1 * zCoords.^2 / 2 / DOF_half_width^2);
    psf = psf .* detection_DOF;
    % apply confocal slit
    confocal_slit = exp(-1 * zCoords.^2 / 2 / confocal_slit_width^2);
    psf = psf .* confocal_slit;

    psf = psf ./ max(psf(:));

    % determine MTF_z as a function of x-axis coordinate
    MTF_z = abs(fftshift(fft(psf,[],1),1));
    % normalisation
    MTF_z = MTF_z ./ max(MTF_z(:));
    MTF_z = MTF_z(ceil(size(MTF_z,1) / 2):end,:); % range, kz = [0:1]
            
    
    %rename outputs
    psfDOF = psf;
    psf = original_psf;
    
            [pupilCoordU,pupilCoordV] = meshgrid([-1.1:0.02:1.1],[-1.1:0.02:1.1]);
            
            figure();
            subplot(5,2,1);
            imagesc(pupilCoordU(1,:),pupilCoordV(:,1),abs(pupilFunctor(pupilCoordU,pupilCoordV))); axis image;
            xlim([-1.1 1.1]);
            ylim([-1.1 1.1]);
            xlabel('u'); ylabel('v');
            colormap(hot(2^8));
            drawnow;shg;
            
            subplot(5,2,2);
            imagesc(pupilCoordU(1,:),pupilCoordV(:,1),angle(pupilFunctor(pupilCoordU,pupilCoordV))); axis image;
            xlim([-1.1 1.1]);
            ylim([-1.1 1.1]);
            xlabel('u'); ylabel('v');
            colormap(hot(2^8));
            drawnow;shg;
            
            % plot light-sheet
            subplot(5,2,[3 4]);
            imagesc(xRange * 1e6,zRange * 1e6,psf); axis image;
            xlabel('x-axis [um]'); ylabel('z-axis [um]');
            colormap(hot(2^8));
            drawnow;shg;

            % plot light-sheet
            subplot(5,2,[5 6]);
            imagesc(xRange * 1e6,zRange * 1e6,psfDOF); axis image;
            xlabel('x-axis [um]'); ylabel('z-axis [um]');
            colormap(hot(2^8));
            drawnow;shg;
            
            % plot MTF
            subplot(5,2,[7 8]);
            imagesc(xRange * 1e6,kzRange,MTF_z);
            title('MTF(x,f_z)');
            xlabel('x [um]');ylabel('f_z');
            ylim([0 1]);
            drawnow;shg;

            % plot thresholded MTF
            subplot(5,2,[9 10]);
            imagesc(xRange * 1e6,kzRange,MTF_z >= 0.05);
            title('MTF(x,f_z)');
            xlabel('x [um]');ylabel('f_z');
            ylim([0 1]);
            drawnow;shg;
end
