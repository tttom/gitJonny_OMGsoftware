zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
xRange = [-10:10:2000] * 1e-6; % propagation axis [metres]
lambda = 333e-9;    % [metres]
% NA = 0.20;
% outputFolder = ...
%             'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\3PELSM\Results_20180419\NA020_lowRes';
% beamTypes = {'Gaussian','Gaussian (2PE)','Gaussian (3PE)'};
beamTypes = {'Bessel10','Bessel10 (2PE)','Bessel10 (3PE)'...
            ,'Bessel5','Bessel5 (2PE)','Bessel5 (3PE)'...
            ,'Bessel2','Bessel2 (2PE)','Bessel2 (3PE)'...
            ,'Bessel1','Bessel1 (2PE)','Bessel1 (3PE)'...
            };

NAs = [0.16,0.17,0.18];

for naIdx = 1:length(NAs)
    NA = NAs(naIdx);
    outputFolder = strcat('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\3PELSM\Results_20180425\NA'...
        ,num2str(NA));
    FastLightSheetProfileSimulation_3PE(zRange,xRange,lambda,NA,outputFolder,beamTypes)
end