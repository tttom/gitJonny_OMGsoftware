%%% name:           CompensatedBesselLightSheetTheory
%%% author:         Jonathan Nylk
%%% date created:   10/08/2016
%%% description:    This function produces equations and plots based on the
%%%                 theoretical parameters for a compensated Bessel
%%%                 light-sheet.
%%%
%%% updates (latest first):
%%%     28/11/2016: Formatted code into function format to allow easy
%%%                 looping with a range of parameters. 
%%%     05/09/2016: Changed k_zRange and zRange to 1D vectors to make code
%%%                 more memory efficient. This reduces the speed of the
%%%                 function slightly.
%%%
%%%
%%% END %%%


function [cartesianPupilFunction,cartesianPupilFunction_amplitude,cartesianPupilFunction_phase,X_cart,Y_cart] = CompensatedBesselLightSheetTheory(beta,sigma)

    if nargin < 1
        beta = 0.05;
    end
    
    if nargin < 2
        sigma = 0.25;
    end

    lambda = 532e-9;    % [metres]
    NA = 0.42;
    objective_tube_lens_length = 0.2;   % [metres]
    magnification = 40;
    ref_index = 1.33;
    % beta = 0.25;    % Bessel parameter
%     beta = 0.05;    % Bessel parameter

    objective_half_cone_angle = asin(NA / ref_index);   % [radians]
    f = objective_tube_lens_length / magnification; % [metres]
    pupilScaleFactor = f * tan(objective_half_cone_angle);  % [metres]
    k = 2 * pi / lambda;    % [metres^-1]
    % R = pupilScaleFactor * (1 - beta / 2);  % [metres]
    % r_annulus = pupilScaleFactor * beta;    % [metres]
    R = pupilScaleFactor; % [metres]
    r_annulus = pupilScaleFactor * (1 - beta / 2); % [metres]


    norm_pupilHalfRange = [0:1/500:1] * sqrt(2);
    % k_pupilHalfRange = norm_pupilHalfRange * 2 * NA / lambda;   % [metres^-1]
    k_pupilHalfRange = norm_pupilHalfRange * 2 * NA / lambda;   % [metres^-1]
    k_r = k_pupilHalfRange; % [metres^-1]
    k_z = sqrt(k^2 - k_r.^2);   % [metres^-1]
    k_r0 = 2 * NA / lambda * (1 - beta / 2);    % [metres^-1]
    % k_r0 = k_r(round(R / pupilHalfRange(length(pupilHalfRange)) * length(k_r)));    % [metres^-1]
    k_z0 = sqrt(k^2 - k_r0^2);  % [metres^-1]

    % z_max = R / atan (r_annulus / f);   % [metres] (McGloin 2005)
    % z_max = R / tan (objective_half_cone_angle);   % [metres] (Durnin 1987)
    % z_max = lambda / 2 / ref_index / beta / (1 - sqrt(1 - (NA / ref_index)^2)); % [metres] (Vettenburg 2014)
    z_max = 500e-6;
    z_max = z_max * 0.05 / beta; % "0.05 / beta" factor to allow propagation length scaling as it has been fudged.
    z = [-10:0.01:10] * z_max;   % [metres]

    %[k_zRange,zRange] = meshgrid(k_z,z);
    k_zRange = k_z;
    clear k_z;
    zRange = z;
    clear z;

    % on-axis intensity function
    % A = (abs(zRange) <= z_max); % flat intensity profile for |z| <= z_max
%     sigma = 2;
    A = exp(sigma .* zRange / (z_max)) .* (abs(zRange) <= z_max); % expential increase for |z| <= z_max

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
%         if rem(kIdx,round(length(k_zRange)/10)) == 0
%             fprintf('Calculating "S", loop %d of %d.\n',[kIdx,length(k_zRange)]);
%         end
    end

    fprintf('Calculating "U".\n');
    for zIdx = 1:length(zRange)
        U(zIdx) = sum(k_zRange .* S .* exp (1i .* k_zRange * zRange(zIdx)),2);
%          if rem(zIdx,round(length(zRange)/10)) == 0
%             fprintf('Calculating "U", loop %d of %d.\n',[zIdx,length(zRange)]);
%         end
    end

    pupilFilter = zeros(size(S));
    pupilFilter = exp(-((k_r - k_r0) / sqrt(2) / (k_r0 / 4)) .^ 2);
    S_Filtered = S .* pupilFilter;
    pupilFilterStandard = (k_r / k_r0 * (1 - beta / 2) >= 1 - beta) .* (k_r / k_r0 * (1 - beta / 2) <= 1);
    pupilFilterStandard = pupilFilterStandard .* max(S_Filtered);
    % S_Filtered = pupilFilter;

    U_Filtered = zeros(size(zRange));

    fprintf('Calculating filtered "U".\n');
    for zIdx = 1:length(zRange)
        U_Filtered(zIdx) = sum(k_zRange .* S_Filtered .* exp (1i .* k_zRange * zRange(zIdx)),2);
%          if rem(zIdx,round(length(zRange)/10)) == 0
%             fprintf('Calculating filtered "U", loop %d of %d.\n',[zIdx,length(zRange)]);
%         end
    end

%     figure;
%     subplot(5,1,1);plot(zRange / z_max,abs(A).^2);
%     title('|A|^2');
%     xlabel('z / z_m_a_x');
%     ylabel('Intensity [a.u.]');
%     xlim([-2 2]);
%     subplot(5,1,2);plot(k_r / (2 * NA / lambda),real(S));
%     title('|S|');
%     ylabel('Amplitude [a.u.]');
%     xlabel('k_r [0:1]');
%     xlim([0 1.4]);
%     subplot(5,1,3);plot(zRange / z_max,abs(U).^2);
%     title('|U|^2');
%     xlabel('z / z_m_a_x');
%     ylabel('Intensity [a.u.]');
%     xlim([-2 2]);
%     subplot(5,1,4);plot(k_r / (2 * NA / lambda),real(S_Filtered),k_r / (2 * NA / lambda),real(pupilFilterStandard));
%     title('|S|');
%     xlim([0 1.4]);
%     xlabel('k_r [0:1]');
%     ylabel('Amplitude [a.u.]');
%     subplot(5,1,5);plot(zRange / z_max,abs(U_Filtered).^2);
%     title('|U|^2');
%     xlabel('z / z_m_a_x');
%     ylabel('Intensity [a.u.]');
%     xlim([-2 2]);
%     drawnow; shg;


    % convert to a 2d pupil function array
    radial_Coord = k_r / (2 * NA / lambda);
    angular_Coord = ([1:length(radial_Coord)] - 1) * 2 * pi / (length(radial_Coord) - 1);
    polarPupilFunction = repmat(S_Filtered,[size(angular_Coord,2),1]); % compensated Bessel beam
%     polarPupilFunction = repmat(abs(pupilFilterStandard),[size(angular_Coord,2),1]); % normal Bessel beam
    polarPupilFunction = reshape(polarPupilFunction,[],1);
    [rad_Coords,ang_Coords] = meshgrid(radial_Coord,angular_Coord);
    [x_polar,y_polar] = pol2cart(ang_Coords,rad_Coords);
    x_polar = reshape(x_polar,[],1);
    y_polar = reshape(y_polar,[],1);
    F = TriScatteredInterp(x_polar,y_polar,polarPupilFunction);
    [X_cart,Y_cart] = meshgrid(([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(x_polar(:)) - min(x_polar(:))) + min(x_polar(:))...
        ,([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(y_polar(:)) - min(y_polar(:))) + min(y_polar(:)));
    cartesianPupilFunction = F(X_cart,Y_cart);
    cartesianPupilFunction(isnan(cartesianPupilFunction)) = 0;
    
    cartesianPupilFunction_amplitude = abs(cartesianPupilFunction);
    cartesianPupilFunction_phase = angle(cartesianPupilFunction);
    
%     figure;
%     subplot(2,1,1);
%     imagesc(([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(x_polar(:)) - min(x_polar(:))) + min(x_polar(:))...
%         ,([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(y_polar(:)) - min(y_polar(:))) + min(y_polar(:))...
%         ,cartesianPupilFunction_amplitude);
%     axis image;
%     title('amplitude');
%     subplot(2,1,2);
%     imagesc(([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(x_polar(:)) - min(x_polar(:))) + min(x_polar(:))...
%         ,([1:length(radial_Coord)] - 1) / (length(radial_Coord) - 1) .* (max(y_polar(:)) - min(y_polar(:))) + min(y_polar(:))...
%         ,cartesianPupilFunction_phase);
%     axis image;
%     title('phase');

end
