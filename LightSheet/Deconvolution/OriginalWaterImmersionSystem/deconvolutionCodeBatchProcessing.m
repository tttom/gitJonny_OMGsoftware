%set up list of folders and corresponding center offsets
inputParams(1).folderName='G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Bead Calibration\400nmSteps_200umRange\2015-08-10 15_06_32.814'; inputParams(1).centerOffset=[0 0];
inputParams(2).folderName='G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Cleared_Tissue\02_DetSide_surface\2015-08-13 12_13_10.668'; inputParams(2).centerOffset=[0 1.2e-5]; 
inputParams(3).folderName='G:\Stored Files\M2_DeconvolutionExampleFiles_CONFIDENTIAL\Non-cleared_Tissue\Middle_09_surface\2015-08-11 16_46_43.518'; inputParams(3).centerOffset=[0 1e-5];





% run processWaterImmersionLightSheetVideos_OfflineCodeOnly on various
% folders with different center offsets
for n=1:length(inputParams)
    folderName=inputParams(n).folderName;
    centerOffset=inputParams(n).centerOffset;
    processWaterImmersionLightSheetVideos_OfflineCodeOnly({folderName},false,false,centerOffset,[-0.0594 -0.0044],[0.0008344 -0.0002917])
end