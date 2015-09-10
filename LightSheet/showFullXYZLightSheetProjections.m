%%% name:           showFullXYZLightSheetProjections
%%% author:         Jonathan Nylk
%%% date created:   14/08/2015
%%% description:    Loads deconvoled light-sheet volumetric data and shows
%%%                 maximum intensity projections along all 3 Cartesian
%%%                 axes. Can also optionally output projections where a
%%%                 number of pixels around the edge of the dataset has
%%%                 been cropped in case of artefacts at the edge.
%%%
%%% updates (latest first):
%%%                 09/09/2015:
%%%                 - Given option to crop a different number of
%%%                 pixels in x,y axes and z axis.
%%%                 - Stopped the cropping method from using the "permute"
%%%                 function as this was using too much memory.
%%%
%%% END %%%

function showFullXYZLightSheetProjections(folderNames,crop,cropPixels_xy,cropPixels_z)

    if nargin<1
        folderNames={'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Bi-axonalNeuron\2015-08-11 17_14_03.104'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\DetSide_01_surface\2015-08-11 14_43_29.666'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\DetSide_06_surface_waterLost\2015-08-11 16_03_23.281'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\DetSide_07_surface\2015-08-11 16_18_50.353'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\IllSide_05_surface\2015-08-11 15_41_16.642'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle_03_surface\2015-08-11 15_12_12.348'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle_09_surface\2015-08-11 16_46_43.518'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Det_02_surface\2015-08-11 14_57_12.904'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Det_08_surface\2015-08-11 16_32_21.885'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Ill_04_surface\2015-08-11 15_26_05.195'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_01\Middle-Ill_10_surface\2015-08-11 17_00_45.692'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_01_-100um\2015-08-11 11_08_58.482'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_01_surface\2015-08-11 10_50_09.802'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\DetSide_04_surface\2015-08-11 13_10_22.513'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_03_-50um\2015-08-11 12_39_39.235'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_03_surface\2015-08-11 12_25_40.095'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_04_surface\2015-08-11 12_55_59.077'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\IllSide_05_surface_longExposure\2015-08-11 13_26_32.087'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_-50um\2015-08-11 11_56_41.187'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_-100um\2015-08-11 11_40_19.870'...
            ,'F:\Stored Files\2015-08-11_Javier\UnclearedTissue_02\Middle_02_surface\2015-08-11 11_27_23.705'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\01_Middle_surface\2015-08-13 12_02_42.848'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_-100um\2015-08-13 12_22_55.902'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_-200um\2015-08-13 12_36_58.950'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\02_DetSide_surface\2015-08-13 12_13_10.668'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\03_Middle_-100um\2015-08-13 14_30_34.294'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\03_Middle_surface\2015-08-13 14_17_39.052'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\04_IllSide_surface\2015-08-13 14_44_31.941'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\05_IllSide_-100um\2015-08-13 15_12_26.620'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\05_IllSide_surface\2015-08-13 14_58_23.297'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\06_Middle-Ill_surface\2015-08-13 15_29_23.649'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\07_Middle-Det_-100um\2015-08-13 15_55_41.448'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\07_Middle-Det_surface\2015-08-13 15_42_54.606'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-100um\2015-08-13 16_22_51.417'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-200um\2015-08-13 16_36_27.459'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_-300um\2015-08-13 16_49_14.294'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\08_DetSide_surface\2015-08-13 16_09_37.889'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-100um\2015-08-13 17_17_04.168'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-200um\2015-08-13 17_30_12.623'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_-300um\2015-08-13 17_43_38.152'...
            ,'F:\Stored Files\2015-08-13_Javier\47TDE-53PBS\09_DetSide_surface\2015-08-13 17_04_02.042'...
            };
    end
    if nargin<2
        crop=true;
    end
    if nargin<3
        cropPixels_xy=40;
    end
    if nargin<4
        cropPixels_z=10;
    end
    
    if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        matFileList=dir(strcat(folderName,'/*.mat'));
        for (fileName={matFileList.name})
            fileName=fileName{1};
            filePathAndName=strcat(folderName,'/',fileName);
            alphaStartPos=strfind(fileName,'alpha')+5;
            alphaEndPos=strfind(fileName,'_beta')-1;
            alphaVal=fileName(alphaStartPos:alphaEndPos);
            storedVariables = whos('-file',filePathAndName);
            if (ismember('restoredDataCube', {storedVariables.name}))
                if str2num(alphaVal)==0
                    load(filePathAndName,'recordedImageStack','xRange','yRange','zRange');
                    dataCube=recordedImageStack;
                    clear recordedImageStack
                else
                    load(filePathAndName,'restoredDataCube','xRange','yRange','zRange');
                    dataCube=restoredDataCube;
                    clear restoredDataCube
                end
                
            xRange=xRange*1e6;
            yRange=yRange*1e6;
            zRange=zRange*1e6;
            xStep=xRange(2)-xRange(1);
            yStep=yRange(2)-yRange(1); %should be equal to x-step
            zStep=zRange(2)-zRange(1);
            
            if crop
                %set first and last cropPixel planes to zero and permute to
                %artificially crop the edge of datacube
%                 tic
%                 for n=1:3
%                     dataCube(:,:,1:cropPixels_xy)=0;
%                     dataCube(:,:,end-cropPixels_xy:end)=0;
%                     dataCube=permute(dataCube,[2,3,1]);
%                 end
%                 toc
                tic
                dataCube(1:cropPixels_xy,:,:)=0;
                dataCube(end-cropPixels_xy:end,:,:)=0;
                dataCube(:,1:cropPixels_xy,:)=0;
                dataCube(:,end-cropPixels_xy:end,:)=0;
                dataCube(:,:,1:cropPixels_z)=0;
                dataCube(:,:,end-cropPixels_z:end)=0;
                toc
            end
            
            %get projection images
            xyProj=squeeze(max(dataCube,[],3));
            xzProj=squeeze(max(dataCube,[],1)).';
            zyProj=squeeze(max(dataCube,[],2));
            clear dataCube %memory efficient
            
            %interpolate z-axis to give square pixels
            zSampling=zStep/xStep;
            newZLength=floor(length(zRange)*zSampling);
            zRange_int=interp1([1:length(zRange)],zRange,[1:newZLength]/zSampling,'linear');
            xzProj_int=interp2(zRange.',yRange,xzProj.',zRange_int.',yRange,'linear');
            xzProj_int=xzProj_int.';
            zyProj_int=interp2(xRange.',zRange,zyProj.',xRange.',zRange_int,'linear');
            zyProj_int=zyProj_int.';
            
            %normalise data range ([0 1])
            xyProj=xyProj./max(xyProj(:));
            xzProj_int=xzProj_int./max(xzProj_int(:));
            zyProj_int=zyProj_int./max(zyProj_int(:));
            
            %make large image of all projections
            allProj=ones(size(xyProj,1)+size(xzProj_int,1)+2,size(xyProj,2)+size(zyProj_int,2)+2);
            allProj(1:size(xyProj,1),1:size(xyProj,2))=xyProj;
            allProj(1:size(xyProj,1),size(xyProj,2)+3:size(allProj,2))=zyProj_int;
            allProj(size(xyProj,1)+3:size(allProj,1),1:size(xyProj,2))=xzProj_int;
            %scalebar in empty section
            scaleBarLength=floor(50/xStep); %50um
            scaleBarHeight=floor(5/yStep);
            allProj(size(xyProj,1)+100:size(xyProj,1)+100+scaleBarHeight...
                ,size(xyProj,2)+200:size(xyProj,2)+200+scaleBarLength)=0;
            
%             % display image
%             figure();
%             imagesc(allProj);axis image;
%             colormap gray;
%             title(sprintf('File: %s, alpha=%s',fileName,alphaVal));
%             drawnow;shg;
            
            % save image
            imwrite(allProj,strcat(folderName,'/',fileName(1:end-4),'Projections.png'));
            end
        end
        
    end
end