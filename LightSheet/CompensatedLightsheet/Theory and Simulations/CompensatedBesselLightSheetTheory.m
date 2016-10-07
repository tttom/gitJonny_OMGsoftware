%%% name:           CompensatedBesselLightSheetTheory
%%% author:         Jonathan Nylk
%%% date created:   10/08/2016
%%% description:    This function produces equations and plots based on the
%%%                 theoretical parameters for a compensated Bessel
%%%                 light-sheet.
%%%
%%% updates (latest first):
%%%     05/09/2016: Changed k_zRange and zRange to 1D vectors to make code
%%%                 more memory efficient. This reduces the speed of the
%%%                 function slightly.
%%%
%%%
%%% END %%%


lambda = 532e-9;    % [metres]
NA = 0.42;
objective_tube_lens_length = 0.2;   % [metres]
magnification = 40;
ref_index = 1.33;
% beta = 0.025;    % Bessel parameter
beta = 0.25;    % Bessel parameter

objective_half_cone_angle = asin(NA / ref_index);   % [radians]
f = objective_tube_lens_length / magnification; % [metres]
pupilScaleFactor = f * tan(objective_half_cone_angle);  % [metres]
k = 2 * pi / lambda;    % [metres^-1]
R = pupilScaleFactor * (1 - beta / 2);  % [metres]
r_annulus = pupilScaleFactor * beta;    % [metres]

% need to find bug in this section. Fourier coordinates scale with the
% spacing of the pupil range. If everything correct this shouldn't
% happen...?
norm_pupilHalfRange = [0:0.0001:1];
pupilHalfRange = norm_pupilHalfRange * pupilScaleFactor;    % [metres]
k_pupilHalfRange = norm_pupilHalfRange * 2 * NA / lambda;   % [metres^-1]
k_r = k_pupilHalfRange; % [metres^-1]
k_z = sqrt(k^2 - k_r.^2);   % [metres^-1]
k_r0 = k_r(round(R / pupilHalfRange(length(pupilHalfRange)) * length(k_r)));    % [metres^-1]
k_z0 = sqrt(k^2 - k_r0^2);  % [metres^-1]

z_max = R / atan (r_annulus / f);   % [metres]
z = [-10:0.01:10] * z_max;   % [metres]

%[k_zRange,zRange] = meshgrid(k_z,z);
k_zRange = k_z;
clear k_z;
zRange = z;
clear z;

% on-axis intensity function
% A = (abs(zRange) <= z_max); % flat intensity profile for |z| <= z_max
sigma = 10;
A = (abs(zRange) <= z_max) .* exp (sigma .* zRange); % linearly increasing profile for |z| <= z_max
% A = (0.5 + exp(0.1 * zRange)) .* (abs(zRange) <= z_max); % expential increase for |z| <= z_max

S = zeros(size(k_zRange));
U = zeros(size(zRange));

% S = 1 / 2 / pi ./ k_z .* sum(A .* exp(1i * k_z0 .* zRange) .* exp (-1i .* k_zRange .* zRange),1);
% figure;plot(k_r / k,real(S));
% 
% U = sum(k_zRange .* repmat(S,[length(z),1]) .*exp(1i .* k_zRange .* zRange),2);
% figure;plot(z,abs(U).^2);

fprintf('Calculating "S".\n');
for kIdx = 1:length(k_zRange)
    S(kIdx) = 1 / 2 / pi / k_zRange(kIdx) * sum(A .* exp (1i * k_z0 .* zRange) .* exp (-1i * k_zRange(kIdx) .* zRange),2);
    if rem(kIdx,round(length(k_zRange)/10)) == 0
        fprintf('Calculating "S", loop %d of %d.\n',[kIdx,length(k_zRange)]);
    end
end

fprintf('Calculating "U".\n');
for zIdx = 1:length(zRange)
    U(zIdx) = sum(k_zRange .* S .* exp (1i .* k_zRange * zRange(zIdx)),2);
     if rem(zIdx,round(length(zRange)/10)) == 0
        fprintf('Calculating "U", loop %d of %d.\n',[zIdx,length(zRange)]);
    end
end

pupilFilter = zeros(size(S));
% pupilFilter = exp(-((k_r - k_r0) / sqrt(2) / (k_r0 / 4)) .^ 2);
pupilFilter = k_r / (k_r0 * pupilScaleFactor / R) >= 1 - beta;
S_Filtered = S .* pupilFilter;

U_Filtered = zeros(size(zRange));

fprintf('Calculating filtered "U".\n');
for zIdx = 1:length(zRange)
    U_Filtered(zIdx) = sum(k_zRange .* S_Filtered .* exp (1i .* k_zRange * zRange(zIdx)),2);
     if rem(zIdx,round(length(zRange)/10)) == 0
        fprintf('Calculating filtered "U", loop %d of %d.\n',[zIdx,length(zRange)]);
    end
end

figure;
subplot(5,1,1);plot(zRange / z_max,abs(A).^2);
title('|A|^2');
xlabel('z / z_m_a_x');
ylabel('Intensity [a.u.]');
xlim([-2 2]);
subplot(5,1,2);plot(k_r / (k_r0 * pupilScaleFactor / R),abs(S));
title('|S|');
ylabel('Amplitude [a.u.]');
xlabel('k_r [0:1]');
xlim([0 1]);
subplot(5,1,3);plot(zRange / z_max,abs(U).^2);
title('|U|^2');
xlabel('z / z_m_a_x');
ylabel('Intensity [a.u.]');
xlim([-2 2]);
subplot(5,1,4);plot(k_r / (k_r0 * pupilScaleFactor / R),abs(S_Filtered));
title('|S|');
xlim([0 1]);
xlabel('k_r [0:1]');
ylabel('Amplitude [a.u.]');
subplot(5,1,5);plot(zRange / z_max,abs(U_Filtered).^2);
title('|U|^2');
xlabel('z / z_m_a_x');
ylabel('Intensity [a.u.]');
xlim([-2 2]);
drawnow; shg;




%                     zRange = [-20:0.01:20] * 1e-6;  % transverse beam axis [metres]
%                     xRange = [-150:1:150] * 1e-6; % propagation axis [metres]
%                     lambda = 532e-9;    % [metres]
%                     NA = 0.42;
%                     % other variable definitions
%                     ref_index = 1.33;
%                     illumination_mag = 40;
%                     illumination_tube_length = 0.2; % [metres]
%                     % set pupil limits
%                     PupilSize = 1;
%                     
%                     % cylindrical simulation so no y-axis needed
%                     yRange = [0] * 1e-6;    % transverse beam axis (in plane of light-sheet) [metres]
%                     
%                     % set hard-aperture at pupil limit
%                     beta = 0.05;   % Bessel parameter
%                     pupilAmplitudeMask = @(U,V) 1 * (sqrt(U.^2 + V.^2) <= PupilSize & (sqrt(U.^2 + V.^2) >= PupilSize - beta));
%             
%                     % set super-Gaussian apodization
%                     sigma_superGaussian = 0;    %super-Gaussian order (must be even!)
%                     superGaussian = @(U,V) exp(-((sqrt(U.^2 + V.^2) - (PupilSize - beta /2)) / sqrt(2) / (beta / 2)).^ sigma_superGaussian);
%             %         superGaussian = @(U,V) exp(-(U+V).^ sigma_superGaussian);
%             
%                     % set compensation mask
%                     sigma = 0;
%                     pupilCompensation = @(U,V) exp(sigma * (sqrt(U.^2 + V.^2) - 1));    % exponential
%             %         pupilCompensation = @(U,V)(sqrt(U.^2 + V.^2).^sigma);    % polynomial
%             
% %                     params = [1.3413, 0.9500, 0.0001, 0.1427, 0.0009, 0.0190, 0.0082]; % sigma = 0;
%             %         params = [1.6364, 0.9500, 0.0001, 0.2039, 0.0008, 0.0334, 0.0054]; % sigma = 10;
%                     
% %                     pupilCompensation = @(U,V) min(max(max(params(1) .* exp(-((sqrt(U.^2 + V.^2) - params(2)) / sqrt(2) / params(3)).^2),params(4) .* exp(-((sqrt(U.^2 + V.^2) - params(2)) / sqrt(2) / params(5)).^2)),params(6) .* exp(-((sqrt(U.^2 + V.^2) - params(2)) / sqrt(2) / params(7)).^2)),1);        
%                     
%                     pupilFunctor = @(U,V) pupilAmplitudeMask(U,V) .* superGaussian(U,V) .* pupilCompensation(U,V);
%                     
%                     pupilRange = [-1:0.001:1];
%                     [U,V] = meshgrid(pupilRange,pupilRange);
%                     
%                     % normalise pupil
%                     pupilNorm = 1 / max(max(pupilFunctor(pupilRange,pupilRange)));
%                     pupilFunctor = @(U,V) pupilFunctor(U,V) * pupilNorm;
%                     
%                     figure;
%                     subplot(1,2,1);
%                     imagesc(pupilRange,pupilRange,abs(pupilFunctor(U,V)));
%                     xlim([-1 1]); ylim([-1 1]); axis square;
%                     xlabel('u-axis'); ylabel('v-axis');
%                     subplot(1,2,2);
%                     plot(pupilRange,abs(pupilFunctor(0,pupilRange)));
%                     xlim([-1 1]); ylim([0 1]); axis square;
%                     xlabel('u-axis'); ylabel('Amplitude [a.u.]');
%                     drawnow;shg;
%                     
%                     % calculate PSF
%                     psf_timer = tic;
%                     psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
%                         ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
%                         ,NA,ref_index,illumination_mag,illumination_tube_length,[]);
%                     psf_time_elapsed = toc(psf_timer);
%                     fprintf('Time taken to calculate light-sheet PSF: %f seconds.\n',psf_time_elapsed);
%                     
%                     % re-orient PSF
%                     psf = squeeze(permute(psf,[1,3,2])).';
%                     
%                     % normalise
%                     psf = psf ./ max(psf(:));
%                     
%                     % plot attenuated light-sheet
%                     figure;
%                     imagesc(xRange * 1e6,zRange * 1e6,psf); axis image;
%                     title('Light-sheet profile');
%                     xlabel('x-axis [um]'); ylabel('z-axis [um]');
%                     drawnow;shg;
%                     
%                     % plot attenuated light-sheet
%                     figure;
%                     plot(xRange * 1e6,squeeze(psf(ceil(size(psf,1) /2),:)));
%                     title('Light-sheet profile');
%                     xlabel('x-axis [um]'); ylabel('Intensity [a.u]');
%                     drawnow;shg;
%                     
%                     
%                     
%                     % calculate PSF
%                     psf_timer = tic;
%                     psf = calcVectorialPsf(zRange,zRange,0,lambda...
%                         ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
%                         ,NA,ref_index,illumination_mag,illumination_tube_length,[]);
%                     psf_time_elapsed = toc(psf_timer);
%                     fprintf('Time taken to calculate light-sheet PSF: %f seconds.\n',psf_time_elapsed);
%                     
%                     % normalise
%                     psf = psf ./ max(psf(:));
%                     
%                     % plot attenuated light-sheet
%                     figure;
%                     imagesc(zRange * 1e6,zRange * 1e6,squeeze(psf)); axis image;
%                     title('Light-sheet profile');
%                     xlabel('y-axis [um]'); ylabel('z-axis [um]');
%                     xlim([-5 5]); ylim([-5 5]);
%                     drawnow;shg;