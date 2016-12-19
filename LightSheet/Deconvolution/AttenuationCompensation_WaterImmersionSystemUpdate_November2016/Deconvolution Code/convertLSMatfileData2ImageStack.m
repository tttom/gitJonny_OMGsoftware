function convertLSMatfileData2ImageStack(folderNames,dataTypes)
%Convert selected data cubes into stacks of .png images (XY cross sections)
%for display in, e.g., Image J.  Input folders contain one or more
%deconvolved data cubes (the .mat file output of the processWaterImmersion
%algorithm) or subfolders containing them.

%Select data sets to convert.  Options are 'lightSheetPsf', 'recorded', and
%'deconvolved'
    if nargin<2
%        dataTypes={'recorded','deconvolved'};
        dataTypes={'lightSheetPsf','recorded','deconvolved'};
    end
    
%Error catching in the case that there is no input folder
    if nargin<1
        Flag_noFolders=true;
    else
        Flag_noFolders=false;
    end
    if exist('folderNames','var')
        if ischar(folderNames)
            folderNames={folderNames};
        end
    end
    if exist('dataTypes','var')
        if ischar(dataTypes)
            dataTypes={dataTypes};
        end
    end  

    if Flag_noFolders
        disp('No input folders')
 
%.png writing code
    else
        for folderNo=1:length(folderNames)
             folderName=folderNames{folderNo};
             fileDir=dir(folderName);
             if length(fileDir)>2 %=2 means empty folder
                 for fileNo=3:length(fileDir)
                      fileName=fileDir(fileNo).name;
                      fullFileName=strcat(folderName,'\',fileName);
                      
                      %Regression to search nested folders for .mat files
                      if fileDir(fileNo).isdir==1
                          convertLSMatfileData2ImageStack(fullFileName,dataTypes);
                      
                      %Check current folder for .mat files
                      else
                          if strcmp(fileName(end-3:end),'.mat')
                             
                              %check matfile contains correct variables
                              Flag_wrongMatfile=0;
                              try
                                  load(fullFileName,'xRange');
                                  clear xRange
                              catch err
                                  Flag_wrongMatfile=1;
                              end
                              
                              %if correct matfile variables, load datacubes
                              %and define folder(s) to save image stacks
                              if ~Flag_wrongMatfile
                                  disp(strcat('Now converting file: ',fullFileName))
                                  
                                  %Handle data types one at a time
                                  for dataNo=1:length(dataTypes)
                                       if strcmp(dataTypes{dataNo},'lightSheetPsf')
                                           disp('Converting lightSheetPsf to image')
                                           load(fullFileName,'lightSheetPsf')
                                           dataCube=squeeze(lightSheetPsf)';
                                           clear lightSheetPsf
                                           outputFolder=strcat(fullFileName(1:end-4),'\lightSheetPsf');
                                       end
                                       if strcmp(dataTypes{dataNo},'deconvolved')
                                           disp('Converting deconvolved data to image stack')
                                           load(fullFileName,'restoredDataCube')
                                           dataCube=restoredDataCube;
                                           clear restoredDataCube
                                           outputFolder=strcat(fullFileName(1:end-4),'\deconvolvedImages');
                                       end
                                       if strcmp(dataTypes{dataNo},'recorded')
                                           disp('Converting recorded data to image stack')
                                           load(fullFileName,'recordedImageStack')
                                           dataCube=recordedImageStack;
                                           clear recordedImageStack
                                           outputFolder=strcat(fullFileName(1:end-4),'\recordedImages');
                                       end
                                       dataMax=max(dataCube(:));
                                       mkdir(outputFolder);
                                       
                                       %Write images of each designated
                                       %data cube's XY slices to
                                       %appropriate folder
                                       for imageNo=1:size(dataCube,3)
                                            outputFullFileName=strcat(outputFolder,'\',dataTypes{dataNo},'Image_',num2str(10000+imageNo),'.png');
                                            outputImage=squeeze(dataCube(:,:,imageNo));
                                            outputImage=outputImage.*(outputImage>0);
                                            outputImage=outputImage./dataMax;
                                            imwrite(outputImage,outputFullFileName,'png','bitdepth',16);
%                                             imwrite(outputImage,outputFullFileName,'png','bitdepth',8);
                                       end
                                       clear dataCube
                                       disp('Image stack written to file')
                                  end
                              end
                          end
                      end
                 end
             end
        end
    end
end