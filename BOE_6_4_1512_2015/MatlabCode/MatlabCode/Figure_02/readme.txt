Files: Figure_02\Analysis_JC_Data_140522_3um
       Figure_02\Analysis_JC_Data_140522_4um
       Figure_02\Avi2Mat_140522_3um
       Figure_02\Avi2Mat_140522_4um
       Figure_02\Tracking_JC_Data_140522
       Figure_02\functions\dftregistration_v1
       Figure_02\functions\expfit
       Figure_02\functions\filenames
       Figure_02\functions\function_noiseRemoval
       Figure_02\functions\transformFourier
Author: Martin V.G. Kristensen
Last Updated: 12/05/2015

These programs are part of a body of work that is published and described 
in Nylk, J. et al, "Development of a graded index microlens based fiber
optical trap and its characterization using principal component
analysis", Biomedical Optics Express 6(4) 1512-1519 (2015) doi:
10.1364/BOE.6.001512.

These programs load video files from trap stiffness measurements, run
them through a principal component analysis (PCA) filtering algrorithm
and then perform particle tracking and power spectrum desnity analysis to
determine the trap siffness of the optical trap. Video files for 3um and
4um diameter beads were used. To run the code for particles of Xum, the
procedure is as follows:

- Run Avi2Mat_140522_Xum to convert the video files to matfiles and to
perform PCA filtering of the data.
- Run Analysis_JC_Data_140522_Xum to perfrom tracking and power spectrum
density analysis.

Other files are subfunctions called by these programs.