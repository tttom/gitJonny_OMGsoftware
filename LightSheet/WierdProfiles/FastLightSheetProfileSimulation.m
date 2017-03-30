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


%%%


function FastLightSheetProfileSimulation(zRange,xRange,lambda,NA)

    %default inputs
    if nargin < 1
%         zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
        zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
    end
    if nargin < 2
%         xRange = [-100:1:100] * 1e-6; % propagation axis [metres]
        xRange = [-250:1:250] * 1e-6;
    end
    if nargin < 3
        lambda = 532e-9;    % [metres]
    end
    if nargin < 4
        NA = 0.42;
    end
    
    % other variable definitions
    ref_index = 1.33;
    illumination_mag = 40;
    illumination_tube_length = 0.2; % [metres]
    
    % cylindrical simulation so no y-axis needed
    yRange = [0] * 1e-6;    % transverse beam axis (in plane of light-sheet) [metres]
%     yRange = zRange;

    % define spatial frequency coords
    kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA * lambda;
    kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]
    kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) * 1 / 2 / xRange(end) / 2 / NA * lambda;
    kxRange = kxRange(ceil(length(kxRange) / 2):end);     % normalised [0:1]

    % set pupil limits
        VPupilSize = 1;       % z-axis
        UPupilSize = 0.01;    % y-axis

    % set hard-aperture at pupil limit
    pupilAmplitudeMaskU = @(U,V) 1 * (U >= -UPupilSize & U <= UPupilSize);
    pupilAmplitudeMaskV = @(U,V) 1 * (V >= -VPupilSize & V <= VPupilSize);
    
    alpha = 3;
    defoc_param = @(x) ref_index / lambda * (1 - sqrt(1 - (NA / ref_index)^2)).* x;
    tilt_param = @(z) illumination_mag * NA * (8.8e-3 / 1.1) / 1.6 / illumination_tube_length / lambda .* z;
    
%     pupilPhaseModulation = @(U,V) exp(2i * pi * alpha.* (V.^3));
    
%     pupilPhaseModulation = @(U,V) exp(2i * pi * alpha.* (V.^3)) .* exp(2i * pi * defoc_param(-50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5.2e-6) .* V) / sqrt(2) ...
%         + exp(2i * pi * -alpha .* (V.^3)) .* exp(2i * pi * defoc_param(+50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5.2e-6) .* V) .* exp(1i * pi) / sqrt(2);

%     pupilPhaseModulationH = @(U,V) exp(2i * pi * alpha.* (V.^3)) .* exp(2i * pi * defoc_param(-50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5.2e-6) .* V);
%     pupilPhaseModulationV = @(U,V) exp(2i * pi * -alpha.* (V.^3)) .* exp(2i * pi * defoc_param(+50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5.2e-6) .* V);

    pupilPhaseModulation = @(U,V) exp(2i * pi * alpha.* (V.^3)) .* exp(2i * pi * defoc_param(-50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5e-6) .* V) / 2 ...
        + exp(2i * pi * -alpha .* (V.^3)) .* exp(2i * pi * defoc_param(+50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5e-6) .* V) / 2 ...
        + exp(2i * pi * -alpha .* (V.^3)) .* exp(2i * pi * defoc_param(-150e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5e-6) .* V) / 2 ...
        + exp(2i * pi * alpha .* (V.^3)) .* exp(2i * pi * defoc_param(+150e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5e-6) .* V) / 2 ;

%     pupilPhaseModulationH = @(U,V) exp(2i * pi * alpha.* (V.^3)) .* exp(2i * pi * defoc_param(-50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5.2e-6) .* V) / sqrt(2) ...
%         + exp(2i * pi * alpha.* (V.^3)) .* exp(2i * pi * defoc_param(150e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(5.2e-6) .* V) / sqrt(2);
%     pupilPhaseModulationV = @(U,V) exp(2i * pi * -alpha.* (V.^3)) .* exp(2i * pi * defoc_param(+50e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5.2e-6) .* V) / sqrt(2) ...
%         + exp(2i * pi * -alpha.* (V.^3)) .* exp(2i * pi * defoc_param(-150e-6) .* (V.^2)) .* exp(2i * pi * tilt_param(-5.2e-6) .* V) / sqrt(2);


%         pupilFunctor = @(U,V) pupilAmplitudeMaskU(U,V) .* pupilAmplitudeMaskV(U,V)...
%             .* pupilPhaseModulation(U,V);

        pupilFunctorH = @(U,V) pupilAmplitudeMaskU(U,V) .* pupilAmplitudeMaskV(U,V)...
            .* pupilPhaseModulationH(U,V);
        pupilFunctorV = @(U,V) pupilAmplitudeMaskU(U,V) .* pupilAmplitudeMaskV(U,V)...
            .* pupilPhaseModulationV(U,V);
        
%         % plot pupil amplitude
%         pupilRangeV = [-1:0.01:1];
%         figure();
%         subplot(1,2,1);
%         plot(pupilRangeV,abs(pupilFunctor(0,pupilRangeV)));
%         xlim([-1 1]); axis square;
%         xlabel('u-axis'); ylabel('Amplitude [a.u.]');
%         subplot(1,2,2);
%         plot(pupilRangeV,angle(pupilFunctor(0,pupilRangeV)));
%         xlim([-1 1]); axis square;
%         xlabel('u-axis'); ylabel('Phase [rad]');
%         drawnow;shg;
        
        % calculate PSF
        psf_timer = tic;
%         psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
%             ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
%             ,NA,ref_index,illumination_mag,illumination_tube_length,[]);
        psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
            ,@(U,V) pupilFunctorH(U,V) / sqrt(2),@(U,V) pupilFunctorV(U,V) / sqrt(2)...
            ,NA,ref_index,illumination_mag,illumination_tube_length,[]);
        psf_time_elapsed = toc(psf_timer);
        fprintf('Time taken to calculate light-sheet PSF: %f seconds.\n',psf_time_elapsed);
        
        % if a 2D pupil function, then 3D beam profile.
        % project along y-axis to determine light-sheet profile.
        if max(yRange~=0)
            psf = sum(psf,1);
        end
        
        % re-orient PSF
        psf = squeeze(permute(psf,[1,3,2])).';
        
        % normalise
        psf = psf ./ max(psf(:));
        
        % plot light-sheet
        figure();
        imagesc(xRange * 1e6,zRange * 1e6,psf); axis image;
        title('Light-sheet profile');
        xlabel('x-axis [um]'); ylabel('z-axis [um]');
        drawnow;shg;
        
        % determine transverse intensity maxima
        [C_max,I_max]=max(psf,[],1);
        
        % plot position and intensity of maxima
        figure();
        subplot(2,1,1)
        plot(xRange * 1e6,C_max);
        title('Light-sheet peak intensity value');
        xlabel('x-axis [um]'); ylabel('Intensity [a.u.]');
        xlim([xRange(1) xRange(end)] * 1e6);
        subplot(2,1,2);
        plot(xRange * 1e6,zRange(I_max) * 1e6);
        title('Light-sheet peak intensity position');
        xlabel('x-axis [um]'); ylabel('z-axis [um]');
        xlim([xRange(1) xRange(end)] * 1e6);
        drawnow;shg;
        
        % determine MTF_z as a function of x-axis coordinate
        MTF_z = abs(fftshift(fft(psf,[],1),1));
        % normalisation
        MTF_z = MTF_z ./ max(MTF_z(:));
        MTF_z = MTF_z(ceil(size(MTF_z,1) / 2):end,:); % range, kz = [0:1]
        
        % plot MTF
        figure();
        subplot(2,1,1);
        imagesc(xRange * 1e6,kzRange,MTF_z);
        title('MTF(x,kz)');
        xlabel('x [um]');ylabel('kz');
        ylim([0 1]);
        subplot(2,1,2);
        imagesc(xRange * 1e6,kzRange,MTF_z >= 0.05);
        title('MTF(x,kz)');
        xlabel('x [um]');ylabel('kz');
        ylim([0 1]);
        drawnow;shg;
        
end
