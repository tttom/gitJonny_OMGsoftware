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
        folderNames={'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\centre1\2015-11-27 15_52_26.595'};
%         folderNames={'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\centre1\2015-11-27 15_52_26.595'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantdetection1\2015-11-27 15_19_29.157'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantdetection2\2015-11-27 15_27_13.131'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantdetection3\2015-11-27 15_31_05.960'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantdetection4\2015-11-27 15_34_56.405'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantillumination1\2015-11-27 15_41_30.519'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantillumination2\2015-11-27 15_45_02.608'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-11-27_Javier_unclearedtissue2\constantillumination3\2015-11-27 15_48_29.364'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\03v1_constdetpath\2015-12-03 15_29_37.251'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\03v2_constdetpath\2015-12-03 15_32_33.922'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\03v3_constdetpath\2015-12-03 15_35_27.759'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\04v1_constdetpath\2015-12-03 15_39_27.254'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\04v2_constdetpath\2015-12-03 15_42_31.384'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\04v3_constdetpath\2015-12-03 15_45_21.334'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\05v1_constexcpath\2015-12-03 15_49_49.354'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\05v2_constexcpath\2015-12-03 15_52_55.943'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\05v3_constexcpath\2015-12-03 15_55_44.299'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\06v1_constexcpath\2015-12-03 15_58_58.736'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\06v2_constexcpath\2015-12-03 16_02_14.476'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\06v3_constexcpath\2015-12-03 16_05_02.092'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\07v1_constdetpath\2015-12-03 16_08_31.823'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\07v2_constdetpath\2015-12-03 16_12_27.940'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\07v3_constdetpath\2015-12-03 16_15_10.338'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\08v1_constexcpath\2015-12-03 16_19_15.583'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\08v2_constexcpath\2015-12-03 16_22_05.792'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\08v3_constexcpath\2015-12-03 16_25_27.867'...
%             ,'D:\PhD Data\LightSheet\ClearedTissueStudy\2015-12-3_Tello_ClearedBrainTissue\08v3-2_constexcpath\2015-12-03 16_32_06.091'...
%             };
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
                %tic
                dataCube(1:cropPixels_xy,:,:)=0;
                dataCube(end-cropPixels_xy:end,:,:)=0;
                dataCube(:,1:cropPixels_xy,:)=0;
                dataCube(:,end-cropPixels_xy:end,:)=0;
                dataCube(:,:,1:cropPixels_z)=0;
                dataCube(:,:,end-cropPixels_z:end)=0;
                %toc
            end
            
            %get projection images
            xyProj=squeeze(max(dataCube,[],3));
            xzProj=squeeze(max(dataCube,[],1)).';
            zyProj=squeeze(max(dataCube,[],2));
            clear dataCube %memory efficient
            
            %interpolate z-axis to give square pixels
            zSampling=zStep/xStep;
            newZLength=floor(length(zRange)*zSampling);
            zRange_int=interp1([1:length(zRange)],zRange,[1:newZLength]/zSampling,'linear','extrap');
            xzProj_int=interp2(zRange.',yRange,xzProj.',zRange_int.',yRange,'linear',0);
            xzProj_int=xzProj_int.';
            zyProj_int=interp2(xRange.',zRange,zyProj.',xRange.',zRange_int,'linear',0);
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