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
%%%
%%%
%%% END %%%

function showFullXYZLightSheetProjections(folderNames,crop,cropPixels)

    if nargin<1
        folderNames='';
    end
    if nargin<2
        crop=true;
    end
    if narging<3
        cropPixels=20;
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
            alphaVal=num2str(fileName(alphaStartPos:alphaEndPos));
            storedVariables = whos('-file',filePathAndName);
            if (ismember('restoredDataCube', {storedVariables.name}))
                if alphaVal==0
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
                xRange=xRange(cropPixels:end-cropPixels);
                yRange=yRange(cropPixels:end-cropPixels);
                zRange=zRange(cropPixels:end-cropPixels);
                dataCube2=dataCube(cropPixels:end-cropPixels,cropPixels:end-cropPixels,cropPixels:end-cropPixels);
                clear dataCube
                dataCube=dataCube2;
                clear dataCube2 %memory efficient?
            end
            
            %get projection images
            xyProj=squeeze(max(restoredDataCube,[],3));
            xzProj=squeeze(max(restoredDataCube,[],1)).';
            zyProj=squeeze(max(restoredDataCube,[],2));
            clear datacube %memory efficient
            
            %interpolate z-axis to give square pixels
            zSampling=zStep/xStep;
            newZLength=floor(length(zRange)*zSampling);
            zRange_int=interp1([1:length(zRange)],zRange,[1:newZLength]/zSampling,'linear');
            xz_Proj_int=interp2(zRange.',yRange,xz_Proj,zRange_int.',yRange,'linear');
            zy_Proj_int=interp2(xRange.',zRange,zy_Proj,xRange.',zRange_int,'linear');
            
            %make large image of all projections
            allProj=ones(size(xyProj,1)+size(xzProj_int,1)+2,size(xyProj,2)+size(zyProj_int,2)+2);
            allProj(1:size(xyProj,1),1:size(xyProj,2))=xyProj;
            allProj(1:size(xyProj,1),size(xyProj,2)+2:end)=zyProj_int;
            allProj(size(xyProj,1)+2:end,1:size(xyProj,2))=xzProj_int;
            %scalebar in empty section
            scaleBarLength=floor(50/xStep);
            scaleBarHeight=floor(5/yStep);
            allProj(size(xyProj,1)+100:size(xyProj,1)+50+scaleBarHeight...
                ,size(xyProj,2)+200:size(xyProj,2)+200+scaleBarLength)=0;
            
            
            figure();
            imagesc(allProj);axis image;
            end
        end
        
    end
end