%%% name:           fitCompensatedBesselPupil2Function
%%% author:         Jonathan Nylk
%%% date created:   08/09/2016
%%% description:    This funciton fits a lorentzian function to user
%%%                 specificed sections of the numerically determined
%%%                 attenuation compnesated Bessel beam pupil function (S)
%%%                 from the script "CompensatedBesselLightSheetTheory.m".
%%%
%%% updates (latest first):
%%%
%%%
%%% END %%%

% pupilCoords = radial pupil coordinate
% Idx = index positions of pupilCoords to fit to
% dataVector = data to fit function to

Functor = @(params,Idx,pupilCoords) params(1) / pi / params(3) ./ (1 + ((pupilCoords(Idx) - params(2)) / params(3)).^2); % Lorentzian
% params(1) = peak amplitude
% params(2) = x-coordinate of peak amplitude
% params(3) = gamma of distribution (width)

% Functor = @(params,Idx,pupilCoords) params(1) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(3)).^2); % Gaussian
% params(1) = peak amplitude
% params(2) = x-coordinate of peak amplitude
% params(3) = width

% Functor = @(params,Idx,pupilCoords) params(1) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(3)).^2) + params(4) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(5)).^2); % 2 Gaussian
% params(1) = peak amplitude of peak1
% params(2) = x-coordinate of both peaks
% params(3) = width of peak 1
% params(4) = peak amplitude of peak2
% params(5) = width of peak2

% Functor = @(params,Idx,pupilCoords) params(1) .* abs(sinc((pupilCoords(Idx) - params(2)) / params(3))); %sinc
% params(1) = peak amplitude
% params(2) = x-coordinate of peak amplitude
% params(3) = width

% Functor = @(params,Idx,pupilCoords) params(1) .* abs(sinc((pupilCoords(Idx) - params(2)) / params(3))) .* params(4) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(5)).^2); %sinc + Gaussian
% params(1) = peak amplitude (sinc)
% params(2) = x-coordinate of peak amplitude
% params(3) = width (sinc)
% params(4) = peak amplitude (Gaussian)
% params(5) = width (Gaussian)

% Functor = @(params,Idx,pupilCoords) min(max(params(1) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(3)).^2),params(4) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(5)).^2)),1); % 2 Gaussian
% params(1) = peak amplitude of peak1
% params(2) = x-coordinate of both peaks
% params(3) = width of peak 1
% params(4) = peak amplitude of peak2
% params(5) = width of peak2

Functor = @(params,Idx,pupilCoords) min(max(max(params(1) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(3)).^2),params(4) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(5)).^2)),params(6) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(7)).^2)),1); % 3 Gaussian
% params(1) = peak amplitude of peak1
% params(2) = x-coordinate of both peaks
% params(3) = width of peak 1
% params(4) = peak amplitude of peak2
% params(5) = width of peak2
% params(4) = peak amplitude of peak3
% params(5) = width of peak3

% Functor = @(params,Idx,pupilCoords) min(params(1) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(3)).^2) + params(4) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(5)).^2) + params(6) .* exp(-((pupilCoords(Idx) - params(2)) / sqrt(2) / params(7)).^2),1); % 3 Gaussian
% params(1) = peak amplitude of peak1
% params(2) = x-coordinate of both peaks
% params(3) = width of peak 1
% params(4) = peak amplitude of peak2
% params(5) = width of peak2
% params(4) = peak amplitude of peak3
% params(5) = width of peak3

FunctorMinimisation = @(params,Idx,pupilCoords,dataVector) sum((Functor(params,Idx,pupilCoords) - (abs(dataVector(Idx)) / max(abs(dataVector(Idx))))).^2);


% determine start and end index positions
dk_r = (k_r(2)-k_r(1)) / (k_r0 * pupilScaleFactor / R);
[~,Idx_start] = max(k_r / (k_r0 * pupilScaleFactor / R) >= 1 - beta);
[~,Idx_end] = max(k_r / (k_r0 * pupilScaleFactor / R) > 1 - dk_r);
Idx = [Idx_start:Idx_end];

% initial_Params = [1 / 300 ,k_r(Idx_start + round(length(Idx) / 2)) / (k_r0 * pupilScaleFactor / R), 10 * dk_r];
% initial_Params = [1, k_r(Idx_start + round(length(Idx) / 2)) / (k_r0 * pupilScaleFactor / R), dk_r];
% initial_Params = [1, k_r(Idx_start + round(length(Idx) / 2)) / (k_r0 * pupilScaleFactor / R) ,dk_r ,1 / 5, 50 * dk_r];
initial_Params = [1, k_r(Idx_start + round(length(Idx) / 2)) / (k_r0 * pupilScaleFactor / R) ,dk_r ,1 / 5, 10 * dk_r, 1 / 10, 50 * dk_r];

[params] = fminsearch(@(params) FunctorMinimisation(params,Idx,k_r / (k_r0 * pupilScaleFactor / R),S),initial_Params)

residual = sum(sqrt((abs(S(Idx)) / max(abs(S(Idx))) - Functor(params,Idx,k_r / (k_r0 * pupilScaleFactor / R))).^2))

figure();
plot(k_r(Idx) / (k_r0 * pupilScaleFactor / R),abs(S(Idx)) / max(abs(S(Idx))),k_r(Idx) / (k_r0 * pupilScaleFactor / R),Functor(params,Idx,k_r / (k_r0 * pupilScaleFactor / R)));
title('fit');

figure();
plot(k_r(Idx) / (k_r0 * pupilScaleFactor / R),abs(S(Idx)) / max(abs(S(Idx))) - Functor(params,Idx,k_r / (k_r0 * pupilScaleFactor / R)));
title('residual');
