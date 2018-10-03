zRange = [-25:0.25:25] * 1e-6;  % transverse beam axis [metres]
xRange = [-100:5:100] * 1e-6;
lambda = 532e-9;    % [metres]
NA_ill = 0.4;
NA_det = 0.4;
outputFolder = [];
%inputPupilFunction
    coeffRangeMin = -10;
    coeffRangeMax = 10;
    phasePolynomialCoeffsAsymm = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* rand([10,1]);
    phasePolynomialCoeffsAsymm(end-1:end) = [0 0]; % set defocus, lateral shift, and constant phase offset to zero
    phasePolynomialCoeffsSymm = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* rand([10,1]);
    phasePolynomialCoeffsSymm(end-2) = 0; % set defocus, lateral shift, and constant phase offset to zero
    phasePolynomialCoeffsSymm(end) = 0; % set defocus, lateral shift, and constant phase offset to zero
    phaseFunctionAsymm = anonymousPolynomialAsymmetric('V',phasePolynomialCoeffsAsymm); %asymmetric versions of polynomial function only
    phaseFunctionSymm = anonymousPolynomialSymmetric('V',phasePolynomialCoeffsSymm); %symmetric versions of polynomial function only
    
%     phaseFunction = anonymousPolynomial('V',phasePolynomialCoeffsSymm); %true polynomial function
%     phaseFunction = @(V) phaseFunctionAsymm(V); %asymmetric versions of polynomial function only
%     phaseFunction = @(V) phaseFunctionSymm(V); %symmetric versions of polynomial function only
    phaseFunction = @(V) phaseFunctionSymm(V) + phaseFunctionAsymm(V); %general case

    inputPupilFunction = @(U,V) exp(2 * pi * 1i .* phaseFunction(V));

tic
[psf,psfDOF,MTF_z] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction);
toc

%%% psf = light-sheet psf
%%% psfDOF = light-sheet psf after depth-of-focus exclusion
%%% MTF_z = modulation transfer function in z-direction along propagation


%%% COST FUNCTIONS
kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA_ill * lambda;
kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]

% ideal MTF
ideal_MTF = zeros([length(kzRange),length(xRange)]);
ideal_MTF_range = [-50,50] * 1e-6;
ideal_MTF(1:find(kzRange <= 1,1,'last')...
    ,find(xRange >= ideal_MTF_range(1),1,'first'):find(xRange <= ideal_MTF_range(2),1,'last')) = 1;

% binary thresholded (simulated) MTF
MTF_threshold = 0.05;
thresholded_MTF_z = MTF_z >= MTF_threshold;
    % remove MTF "ripples"
    for x_idx = 1:length(xRange)
        test_vector = thresholded_MTF_z(:,x_idx);
        cuttoff_idx = find(test_vector == 0,1,'first');
        thresholded_MTF_z(cuttoff_idx:end,x_idx) = 0;
    end
    
    
% maximise continuous MTF
cf1 = min(min((ideal_MTF - MTF_z).^2));

% maximise binary MTF
cf2 = min(min((ideal_MTF - thresholded_MTF_z).^2));

% maximise light efficiency (localise light-sheet within DOF)
cf3 = min(min((psfDOF - psf).^2));

k_res = 1; %weighting factor for resolution
k_eff = 1; %weighting factor for efficiency

% maximise continuous MTF + light efficiency
cf4 = min(min(k_res .* (ideal_MTF - MTF_z).^2 + k_eff .* (psfDOF - psf).^2));

% maximise binary MTF + light efficiency
cf5 = min(min(k_res .* (ideal_MTF - thresholded_MTF_z).^2 + k_eff .* (psfDOF - psf).^2));