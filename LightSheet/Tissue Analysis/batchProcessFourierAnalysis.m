%batchProcessFourierAnalysis

% %% Datasets of neurons expressing mCherry (2048x1024 pixel camera frame)
%
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\DetSide_01_surface\2015-08-11 14_43_29.666'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\DetSide_07_surface\2015-08-11 16_18_50.353'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\IllSide_05_surface\2015-08-11 15_41_16.642'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle_03_surface\2015-08-11 15_12_12.348'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle_09_surface\2015-08-11 16_46_43.518'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Det_02_surface\2015-08-11 14_57_12.904'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Det_08_surface\2015-08-11 16_32_21.885'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Ill_04_surface\2015-08-11 15_26_05.195'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Ill_10_surface\2015-08-11 17_00_45.692'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_01_surface\2015-08-11 10_50_09.802'...
%     ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_01_-100um\2015-08-11 11_08_58.482'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_04_surface\2015-08-11 13_10_22.513'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_03_surface\2015-08-11 12_25_40.095'...
%     ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_03_-50um\2015-08-11 12_39_39.235'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_04_surface\2015-08-11 12_55_59.077'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_05_surface_longExposure\2015-08-11 13_26_32.087'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_surface\2015-08-11 11_27_23.705'...
%     ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_-50um\2015-08-11 11_56_41.187'...
%     ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_-100um\2015-08-11 11_40_19.870'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\01_Middle_surface\2015-08-13 12_02_42.848'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_surface\2015-08-13 12_13_10.668'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_-100um\2015-08-13 12_22_55.902'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_-200um\2015-08-13 12_36_58.950'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\03_Middle_surface\2015-08-13 14_17_39.052'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\03_Middle_-100um\2015-08-13 14_30_34.294'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\04_IllSide_surface\2015-08-13 14_44_31.941'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\05_IllSide_surface\2015-08-13 14_58_23.297'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\05_IllSide_-100um\2015-08-13 15_12_26.620'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\06_Middle-Ill_surface\2015-08-13 15_29_23.649'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\07_Middle-Det_surface\2015-08-13 15_42_54.606'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\07_Middle-Det_-100um\2015-08-13 15_55_41.448'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_surface\2015-08-13 16_09_37.889'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-100um\2015-08-13 16_22_51.417'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-200um\2015-08-13 16_36_27.459'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-300um\2015-08-13 16_49_14.294'});
% FourierAnalysis_LSMTissueData({'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_surface\2015-08-13 17_04_02.042'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-100um\2015-08-13 17_17_04.168'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-200um\2015-08-13 17_30_12.623'...
%     ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-300um\2015-08-13 17_43_38.152'});


% %% Datasets of tissue injected with fluorescent beads (2048x1024 pixel camera frame)
%
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\center\2015-10-19 12_57_39.967'});
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\2015-10-19_Javier_beadbrain_contd\center2\2015-10-19 15_33_22.510'});
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\detectionaberrations\2015-10-19 13_13_22.215'});
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\2015-10-19_Javier_beadbrain_contd\detectionab2\2015-10-19 15_00_22.047'});
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\excitationaberration\2015-10-19 12_43_12.138'});
% FourierAnalysis_LSMTissueData({'G:\LSM Data\2015-10-19_Javier_beadbrainsample1\2015-10-19_Javier_beadbrain_contd\excitationab2\2015-10-19 14_44_01.608'});

%% Datasets of tissue injected with fluorescent beads (2048x1024 pixel camera frame)
    % PERFORM fliplr(restoredDataCube) IN ForuierAnalysis_LSMTissueData.m
    % ON THESE DATASETS BECAUSE THE CAMERA IMAGE HAS FLIPPED, OTHERWISE
    % ROTATION WILL BE INCORRECT!!!
    % ALSO NEED TO CHANGE LINE 51 to 
    % fileName=strcat(folderName,'/recording0_lambda532nm_alpha-7_beta100.mat');
    % BECAUSE THE ALPHA VALUE HAS CHANGED SIGN!!!!

FourierAnalysis_LSMTissueData({'F:\LSM Files\2016-02-17_Javier-beads-cleared\constdet_top\2016-02-17 14_39_37.490'...
    ,'F:\LSM Files\2016-02-17_Javier-beads-cleared\constdet_2\2016-02-17 14_41_48.556'...
    ,'F:\LSM Files\2016-02-17_Javier-beads-cleared\constdet_3\2016-02-17 14_45_17.552'...
    });
FourierAnalysis_LSMTissueData({'F:\LSM Files\2016-02-17_Javier-beads-cleared\constexc_top\2016-02-17 15_28_08.952'...
    ,'F:\LSM Files\2016-02-17_Javier-beads-cleared\constexc_2\2016-02-17 15_30_42.538'...
    ,'F:\LSM Files\2016-02-17_Javier-beads-cleared\constexc_3\2016-02-17 15_35_40.624'...
    });

