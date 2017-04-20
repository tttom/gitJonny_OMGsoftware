%set up list of folders and corresponding center offsets
% inputParams(4).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs27pt5cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 11_50_38.298';
% inputParams(4).centerOffset=[0 -5.6]*1e-6;
% inputParams(4).sampleAttenuation=27.5*100;
% 
% inputParams(3).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs27pt5cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 11_57_01.018';
% inputParams(3).centerOffset=[0 -4.6]*1e-6;
% inputParams(3).sampleAttenuation=27.5*100;
% 
% inputParams(1).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 12_33_32.528';
% inputParams(1).centerOffset=[0 -7.3]*1e-6;
% inputParams(1).sampleAttenuation=55*100;
% 
% inputParams(2).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_noAttenuationDeconvolution\2017-01-25 12_36_42.195';
% inputParams(2).centerOffset=[0 -3.5]*1e-6;
% inputParams(2).sampleAttenuation=55*100;

% inputParams(1).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_02\2017-04-14 15_18_02.825';
% inputParams(1).centerOffset=[0 -1.5]*1e-6;
% inputParams(1).sampleAttenuation=50*100;
% 
% inputParams(2).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_04\2017-04-14 15_41_54.028';
% inputParams(2).centerOffset=[0 -1.5]*1e-6;
% inputParams(2).sampleAttenuation=85*100;
% 
% inputParams(3).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_05\2017-04-14 15_46_28.226';
% inputParams(3).centerOffset=[0 1.5]*1e-6;
% inputParams(3).sampleAttenuation=85*100;
% 
% inputParams(4).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_06\2017-04-14 15_51_53.940';
% inputParams(4).centerOffset=[0 0]*1e-6;
% inputParams(4).sampleAttenuation=85*100;
% 
% inputParams(5).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_11\2017-04-14 16_17_27.871';
% inputParams(5).centerOffset=[0 -4]*1e-6;
% inputParams(5).sampleAttenuation=25*100;
% 
% inputParams(6).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_12\2017-04-14 16_23_26.727';
% inputParams(6).centerOffset=[0 -4]*1e-6;
% inputParams(6).sampleAttenuation=30*100;
% 
% inputParams(7).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_13\2017-04-14 16_32_05.751';
% inputParams(7).centerOffset=[0 -35]*1e-6;
% inputParams(7).sampleAttenuation=75*100;
% 
% inputParams(8).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_15\2017-04-14 16_41_33.327';
% inputParams(8).centerOffset=[0 -3]*1e-6;
% inputParams(8).sampleAttenuation=75*100;
% 
% inputParams(9).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_16\2017-04-14 16_51_22.343';
% inputParams(9).centerOffset=[0 -5]*1e-6;
% inputParams(9).sampleAttenuation=50*100;
% 
% inputParams(10).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_17\2017-04-14 16_56_24.888';
% inputParams(10).centerOffset=[0 -2.5]*1e-6;
% inputParams(10).sampleAttenuation=80*100;
% 
% inputParams(11).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_18\2017-04-14 17_00_15.377';
% inputParams(11).centerOffset=[0 -5]*1e-6;
% inputParams(11).sampleAttenuation=80*100;

% inputParams(1).folderName='H:\NEW_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_13\2017-04-14 16_32_05.751';
% inputParams(1).centerOffset=[0 -35]*1e-6;
% inputParams(1).sampleAttenuation=75*100;

inputParams(1).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_modifiedDeconvolutionProcedure - 2017-04-20\2017-01-25 12_33_32.528';
inputParams(1).centerOffset=[0 -7.3]*1e-6;
inputParams(1).sampleAttenuation=55*100;
inputParams(1).scanShear = [-0.0594,-0.0044];
inputParams(1).perspectiveScaling = [0.0008344,-0.0002917];

inputParams(2).folderName='F:\NEW_SYSTEM_RESULTS\2017-01-25_AttnComp\abs55cm-1\25msInt_modifiedDeconvolutionProcedure - 2017-04-20\2017-01-25 12_36_42.195';
inputParams(2).centerOffset=[0 -3.5]*1e-6;
inputParams(2).sampleAttenuation=55*100;
inputParams(2).scanShear = [-0.0594,-0.0044];
inputParams(2).perspectiveScaling = [0.0008344,-0.0002917];

% run processWaterImmersionLightSheetVideos_OfflineCodeOnly on various
% folders with different center offsets
for n=1:length(inputParams)
    folderName=inputParams(n).folderName;
    centerOffset=inputParams(n).centerOffset;
    sampleAttenuation=inputParams(n).sampleAttenuation;
    scanShear = inputParams(n).scanShear;
    perspectiveScaling = inputParams(n).perspectiveScaling;
%     sampleAttenuation=0;
%     processWaterImmersionLightSheetVideos_OfflineCodeOnly({folderName},true,false,centerOffset,[0.0002 0.0302],[0.5 1.2]*1e-3,false,sampleAttenuation)
    processWaterImmersionLightSheetVideos_OfflineCodeOnly({folderName},true,false,centerOffset,scanShear,perspectiveScaling,false,sampleAttenuation)
    
%     convertLSMatfileData2ImageStack(folderName);
    
end


% processWaterImmersionLightSheetVideos_OfflineCodeOnly(folderNames,reprocessData,deleteAvi,centerOffset,scanShear,perspectiveScaling,showFigures,sampleAttenuation)