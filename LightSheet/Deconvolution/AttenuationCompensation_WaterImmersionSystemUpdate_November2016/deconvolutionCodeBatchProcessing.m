%set up list of folders and corresponding center offsets
inputParams(4).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs27pt5cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 11_50_38.298';
inputParams(4).centerOffset=[0 -5.6]*1e-6;
inputParams(4).sampleAttenuation=27.5*100;

inputParams(3).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs27pt5cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 11_57_01.018';
inputParams(3).centerOffset=[0 -4.6]*1e-6;
inputParams(3).sampleAttenuation=27.5*100;

inputParams(1).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 12_33_32.528';
inputParams(1).centerOffset=[0 -7.3]*1e-6;
inputParams(1).sampleAttenuation=55*100;

inputParams(2).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 12_36_42.195';
inputParams(2).centerOffset=[0 -3.5]*1e-6;
inputParams(2).sampleAttenuation=55*100;



% run processWaterImmersionLightSheetVideos_OfflineCodeOnly on various
% folders with different center offsets
for n=1:length(inputParams)
    folderName=inputParams(n).folderName;
    centerOffset=inputParams(n).centerOffset;
%     sampleAttenuation=inputParams(n).sampleAttenuation;
    sampleAttenuation=0;
    processWaterImmersionLightSheetVideos_OfflineCodeOnly({folderName},true,false,centerOffset,[-0.0594 -0.0044],[0.0008344 -0.0002917],false,sampleAttenuation)
end



% processWaterImmersionLightSheetVideos_OfflineCodeOnly(folderNames,reprocessData,deleteAvi,centerOffset,scanShear,perspectiveScaling,showFigures,sampleAttenuation)