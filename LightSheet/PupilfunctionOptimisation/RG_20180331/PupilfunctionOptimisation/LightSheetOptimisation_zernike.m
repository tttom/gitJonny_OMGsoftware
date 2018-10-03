function [cf] = LightSheetOptimisation_zernike(inputcoeffs,make_fig)
%make_fig = 1;

%default inputs
    if nargin < 2
        make_fig = 1;
    end
    if nargin < 1
        inputcoeffs=rand([1,21]);
    end

zRange = [-25:0.25:25] * 1e-6;  % transverse beam axis [metres]
xRange = [-100:5:100] * 1e-6;
lambda = 532e-9;    % [metres]
NA_ill = 0.4; %might be optimised, increased NA_ill should increase axial resolution (STAY BELOW 0.8)
NA_det = 0.4; %might be optimised, increased NA_det gives better collection efficiency (STAY BELOW 0.8)
outputFolder = [];


    %inputPupilFunction
    coeffRangeMin = -1;
    coeffRangeMax = 1;
    
    
        %%% for coefficinets to be zero: set them as:
        %%% -1 * coeffRangeMin / (coeffRangeMax - coeffRangeMin)
    
        for coeffIdx = 1:length(inputcoeffs)
            if inputcoeffs(1,coeffIdx) == 0
                inputcoeffs(1,coeffIdx) = -1 * coeffRangeMin / (coeffRangeMax - coeffRangeMin);
            end
        end
        

    %%% Probably want to put this outside of any optimisation loop to save it recalculating the zernike modes every time (BEGIN)
    %sets and generates the number of zernike polynomials to use
    n_max = 5;
    [ns,ms,~] = calculateNumberZernikeTerms(n_max);
    anonZernikeTerms = generateAnonymousZernikeTerms(ns,ms,'U','V');
    %%% Probably want to put this outside of any optimisation loop to save it recalculating the zernike modes every time (END)
    
    inputZernikeCoeffs = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* inputcoeffs;
    fullZernikePolynomial = constructAnonymousZernikePolynomial(anonZernikeTerms,inputZernikeCoeffs);
    
    inputPupilFunction = @(U,V) exp(2 * pi * 1i .* fullZernikePolynomial(U,V));

% calculate aberration-free PSF    
[psf,psfDOF,MTF_z] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction,make_fig);


Z_4_0 = anonymousZernike('U','V',4,0);
% calculate PSF with positive primary spherical aberration
pupilFunctor_Z_4_0_p = @(U,V) exp(2 * pi * 1i .* Z_4_0(U,V) * 1);
inputPupilFunction_Z_4_0_p = @(U,V) inputPupilFunction(U,V) .* pupilFunctor_Z_4_0_p(U,V);
[psf_Z_4_0_p,psfDOF_Z_4_0_p,MTF_z_Z_4_0_p] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction_Z_4_0_p,make_fig);
% calculate PSF with negative primary spherical aberration
pupilFunctor_Z_4_0_n = @(U,V) exp(2 * pi * 1i .* Z_4_0(U,V) * -1);
inputPupilFunction_Z_4_0_n = @(U,V) inputPupilFunction(U,V) .* pupilFunctor_Z_4_0_n(U,V);
[psf_Z_4_0_n,psfDOF_Z_4_0_n,MTF_z_Z_4_0_n] = LightSheetSimulation_PupilfunctionPropagation(zRange,xRange,lambda,NA_ill,NA_det,outputFolder,inputPupilFunction_Z_4_0_n,make_fig);

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
cf(5) = mean(mean((ideal_intensity_profile - thresholded_intensity_profile).^2));

% minimise effect of aberrations
cf(6) = mean(mean((psf - psf_Z_4_0_p).^2));
cf(7) = mean(mean((psf - psf_Z_4_0_n).^2));


k_res = .5; %weighting factor for resolution
k_eff = .5; %weighting factor for efficiency
k_int = .5; %weighting factor for int

% maximise continuous MTF + light efficiency
jcf(1) = k_res*cf(1) + k_eff*cf(3);

% maximise binary MTF + light efficiency
jcf(2) = k_res*cf(2) + k_eff*cf(3);
%<A*|B>