function [cf] = LightSheetOptimisation_polynomial(inputcoeffs,make_fig)
%make_fig = 1;

%default inputs
    if nargin < 2
        make_fig = 0;
    end
    if nargin < 1
        inputcoeffs=[0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
    end

zRange = [-25:0.25:25] * 1e-6;  % transverse beam axis [metres]
xRange = [-100:5:100] * 1e-6;
lambda = 532e-9;    % [metres]
NA_ill = 0.4; %might be optimised, increased NA_ill should increase axial resolution (STAY BELOW 0.8)
NA_det = 0.4; %might be optimised, increased NA_det gives better collection efficiency (STAY BELOW 0.8)
outputFolder = [];
%inputPupilFunction
    coeffRangeMin = -10;
    coeffRangeMax = 10;
    
    
        %%% for coefficinets to be zero: set them as:
        %%% -1 * coeffRangeMin / (coeffRangeMax - coeffRangeMin)
    
        for coeffIdx = 1:length(inputcoeffs)
            if inputcoeffs(coeffIdx,1) == 0
                inputcoeffs(coeffIdx,1) = -1 * coeffRangeMin / (coeffRangeMax - coeffRangeMin);
            end
        end
        
    
    phasePolynomialCoeffsAsymm = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* inputcoeffs(1:size(inputcoeffs,1)/2);
    phasePolynomialCoeffsAsymm(end-1:end) = [0 0]; % set lateral shift and constant phase offset to zero
    phasePolynomialCoeffsSymm = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* inputcoeffs(1+size(inputcoeffs,1)/2:size(inputcoeffs,1));
    phasePolynomialCoeffsSymm(end-2) = 0; % set defocus to zero
    phasePolynomialCoeffsSymm(end) = 0; % set constant phase offset to zero
    phaseFunctionAsymm = anonymousPolynomialAsymmetric('V',phasePolynomialCoeffsAsymm); %asymmetric versions of polynomial function only
    phaseFunctionSymm = anonymousPolynomialSymmetric('V',phasePolynomialCoeffsSymm); %symmetric versions of polynomial function only
    
%     phaseFunction = anonymousPolynomial('V',phasePolynomialCoeffsSymm); %true polynomial function
%     phaseFunction = @(V) phaseFunctionAsymm(V); %asymmetric versions of polynomial function only
%     phaseFunction = @(V) phaseFunctionSymm(V); %symmetric versions of polynomial function only
    phaseFunction = @(V) phaseFunctionSymm(V) + phaseFunctionAsymm(V); %general case

    inputPupilFunction = @(U,V) exp(2 * pi * 1i .* phaseFunction(V));

[psf,psfDOF,MTF_z] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction,make_fig);

%%% psf = light-sheet psf
%%% psfDOF = light-sheet psf after depth-of-focus exclusion
%%% MTF_z = modulation transfer function in z-direction along propagation


%%% COST FUNCTIONS
kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA_ill * lambda;
kzRange = kzRange(ceil(length(kzRange) / 2):end);   % normalised [0:1]

% ideal MTF
ideal_MTF = zeros([length(kzRange),length(xRange)]);
% ideal_MTF_range = [-50,50] * 1e-6;   %Might change this later
ideal_MTF(1:find(kzRange <= 1,1,'last'),:)=1;
    %,find(xRange >= ideal_MTF_range(1),1,'first'):find(xRange <= ideal_MTF_range(2),1,'last')) = 1;

% binary thresholded (simulated) MTF
MTF_threshold = 0.05;
thresholded_MTF_z = MTF_z >= MTF_threshold;
    % remove MTF "ripples"
    for x_idx = 1:length(xRange)
        test_vector = thresholded_MTF_z(:,x_idx);
        cuttoff_idx = find(test_vector == 0,1,'first');
        thresholded_MTF_z(cuttoff_idx:end,x_idx) = 0;
    end

% intensity profile
int_profile = max(psfDOF,[],1);
[~,IntRange] = meshgrid(xRange,[0:0.01:1]);
thresholded_intensity_profile = zeros(size(IntRange));
ideal_intensity_profile = ones(size(IntRange));
for x_idx = 1:length(xRange)
    thresholded_intensity_profile(:,x_idx) = IntRange(:,x_idx) <= int_profile(x_idx);
end
    
%     load('airy.mat')
%     % plot thresholded MTF
%     figure;
%             imagesc(xRange * 1e6,kzRange,thresholded_MTF_z-0.5*airy(:,:));
%             title('filtered MTF(x,f_z)');
%             xlabel('x [um]');ylabel('f_z');
%             ylim([0 1]);
%             drawnow;shg;
    
    
% maximise continuous MTF
cf(1) = mean(mean((ideal_MTF - MTF_z).^2));

% maximise binary MTF
cf(2) = mean(mean((ideal_MTF - thresholded_MTF_z).^2));

% maximise continuous (but thresholded) MTF
cf(3) = mean(mean((ideal_MTF - (thresholded_MTF_z .* MTF_z)).^2));

% maximise light efficiency (localise light-sheet within DOF)
cf(4) = mean(mean((psfDOF - psf).^2));

% maximise intensity envelope
cf(5) = mean(mean(ideal_intensity_profile - thresholded_intensity_profile).^2);


k_res = .5; %weighting factor for resolution
k_eff = .5; %weighting factor for efficiency
k_int = .5; %weighting factor for int

% maximise continuous MTF + light efficiency
jcf(1) = k_res*cf(1) + k_eff*cf(3);

% maximise binary MTF + light efficiency
jcf(2) = k_res*cf(2) + k_eff*cf(3);
%<A*|B>