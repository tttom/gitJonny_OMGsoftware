% Attenuation-compensation Bessel Fourier Synthesis test

z = [-40:20:40] * 1e-6;

% Generate Pupil mask
[cartesianPupilFunction,~,~,~,~] = CompensatedBesselLightSheetTheory();

blankMask = zeros(size(cartesianPupilFunction));

fullPsfs = zeros([501,501,length(z),size(cartesianPupilFunction,1)]);


% run beam propagation on full pupil
    [~,~,transverse_Range,longitudinal_Range,fullPsf_normal]...
        = besselPsfSimulation_ZYang(cartesianPupilFunction,z);


for sIdx = 1:size(cartesianPupilFunction,1)
    
    sIdx
    
    
    % take slice out of pupil mask
    pupilMask = blankMask;
    pupilMask(sIdx,:) = cartesianPupilFunction(sIdx,:);
    
    
    % run beam propagation on slice
    [~,~,transverse_Range,longitudinal_Range,fullPsf]...
        = besselPsfSimulation_ZYang(pupilMask,z);
    
    fullPsfs(:,:,:,sIdx) = fullPsf(2501-250:2501+250,2501-250:2501+250,:);
end

clear fullPsf

%sum all slices to give average sheet
fullPsf_FourierSynthesis = sum(fullPsfs,4,'omitnan');