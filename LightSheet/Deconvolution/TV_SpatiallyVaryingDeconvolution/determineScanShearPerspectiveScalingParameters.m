function [dataStruct]=determineScanShearPerspectiveScalingParamters(recordedImageStack,xRange,yRange,zRange)
    % From a 3d datacube "recordedImageStack", allows determination of
    % scanShear and perspectiveScaling paramters by plotting the angle of
    % Airy PSF tails relative to the z-axis as a function of x and y.
    % Use on datasets of fluorescent beads only.
    
    if nargin<1
%        % "Mock-up" some data
%        xRange=single(1:512);
%        yRange=single(1:1024);
%        zRange=single(1:250);
%        xRange=xRange*1e-6;
%        yRange=yRange*1e-6;
%        zRange=zRange*1e-6;
%        recordedImageStack=randn([length(xRange),length(yRange),length(zRange)],'single');
%        recordedImageStack=recordedImageStack.*(recordedImageStack>0);
%        recordedImageStack=recordedImageStack./max(recordedImageStack(:));
%        recordedImageStack=recordedImageStack.*(recordedImageStack>0.7);
       
       % OR use real data
%        xRange=single(1:1024);
%        yRange=single(1:2048);
%        zRange=single(1:500);
%        fileName='E:\Stored Files\2015-04-21_TelloJavier\01_BeadCalibration\2015-04-21 14_37_15.684\recording0_lambda532nm_alpha7_beta100.mat';
%         fileName='F:\Stored Files\2015-02-18_Calibration\600nm_beads\2015-02-18 16_18_30.459\Airy\recording0_lambda532nm_alpha7_beta100.mat';
        fileName='F:\Stored Files\2015-05-29_BeadTest\Aperture_06\2015-05-29 12_38_50.809\recording0_lambda532nm_alpha7_beta100.mat';
        load(fileName,'recordedImageStack','xRange','yRange','zRange');
%        recordedImageStack=single(recordedImageStack);
%        recordedImageStack=recordedImageStack(1:512,1:1024,1:250); %save time and memory
    end
    
    xRange=xRange*1e6;
    yRange=yRange*1e6;
    zRange=zRange*1e6;
    
    xSpacing=xRange(2)-xRange(1);
    ySpacing=yRange(2)-yRange(1);
    zSpacing=zRange(2)-zRange(1);
    
    
    % Generate initial grid for localised searching for beads.
    % Let x and y refer to matlab coordinate system for clarity.
    noDataPoints=4000;
    xStartingCoords=rand([1,noDataPoints],'single');
    yStartingCoords=rand([1,noDataPoints],'single');
    zStartingCoords=rand([1,noDataPoints],'single');
    xStartingCoords=round(1+(length(xRange)-1)*xStartingCoords);
    yStartingCoords=round(1+(length(yRange)-1)*yStartingCoords);
    zStartingCoords=round(1+(length(zRange)-1)*zStartingCoords);
    
    
    % Generate a data structure for storing data and populate some fields;
    dataStruct=[];
    for n=1:noDataPoints
       dataStruct(n).beadNumber=n;
       dataStruct(n).startingXCoord=xStartingCoords(n);
       dataStruct(n).startingYCoord=yStartingCoords(n);
       dataStruct(n).startingZCoord=zStartingCoords(n);
       dataStruct(n).searchRangeX=[max(1,dataStruct(n).startingXCoord-50),min(length(xRange),dataStruct(n).startingXCoord+50)];
       dataStruct(n).searchRangeY=[max(1,dataStruct(n).startingYCoord-50),min(length(yRange),dataStruct(n).startingYCoord+50)];
       dataStruct(n).searchRangeZ=[max(1,dataStruct(n).startingZCoord-50),min(length(zRange),dataStruct(n).startingZCoord+50)];
       dataStruct(n).Flag_noBead=0;
       dataStruct(n).maxXCoord=[];
       dataStruct(n).maxYCoord=[];
       dataStruct(n).maxZCoord=[];
       dataStruct(n).analysisRangeX=[];
       dataStruct(n).analysisRangeY=[];
       dataStruct(n).analysisRangeZ=[];
       dataStruct(n).Flag_duplicatedBead=0;
       dataStruct(n).psfXZTilt=[];
       dataStruct(n).psfYZTilt=[];
    end
    clear xStartingCoords yStartingCoords zStartingCoords
    
    
    for n=1:noDataPoints
        croppedDataCube=recordedImageStack(dataStruct(n).searchRangeX(1):dataStruct(n).searchRangeX(2),dataStruct(n).searchRangeY(1):dataStruct(n).searchRangeY(2),dataStruct(n).searchRangeZ(1):dataStruct(n).searchRangeZ(2));
        meanVal=mean(croppedDataCube(:));
        
        [maxVal,maxPos]=max(croppedDataCube(:));
%             % Plots for testing
%             figure(1);
%             subplot(3,1,1);imagesc(dataStruct(n).searchRangeY,dataStruct(n).searchRangeZ,squeeze(max(croppedDataCube,[],1))');axis image;
%             xlabel('y [pixels]');ylabel('z [pixels]');title(strcat('yz-proj - dataPoint #:',num2str(n)));
%             subplot(3,1,2);imagesc(dataStruct(n).searchRangeX,dataStruct(n).searchRangeZ,squeeze(max(croppedDataCube,[],2))');axis image;
%             xlabel('x [pixels]');ylabel('z [pixels]');title(strcat('xz-proj - dataPoint #:',num2str(n)));
%             subplot(3,1,3);imagesc(dataStruct(n).searchRangeX,dataStruct(n).searchRangeY,squeeze(max(croppedDataCube,[],3))');axis image;
%             xlabel('x [pixels]');ylabel('y [pixels]');title(strcat('xy-proj - dataPoint #:',num2str(n)));
%             %
        if maxVal==0
            dataStruct(n).Flag_noBead=1;
%             % Plots for testing
%             figure(2);
%             subplot(1,1,1);
%             %
        elseif maxVal<2*meanVal;
            dataStruct(n).Flag_noBead=1;
%             % Plots for testing
%             figure(2);
%             subplot(1,1,1);
%             %
        else
            dataStruct(n).maxZCoord=floor(maxPos/size(croppedDataCube,2)/size(croppedDataCube,1))+dataStruct(n).searchRangeZ(1);
            if rem(maxPos,size(croppedDataCube,2)*size(croppedDataCube,1))==0
                maxPos=size(croppedDataCube,2)*size(croppedDataCube,1);
                dataStruct(n).maxZCoord=dataStruct(n).maxZCoord-1;
            else
                maxPos=rem(maxPos,size(croppedDataCube,2)*size(croppedDataCube,1));
            end
            dataStruct(n).maxYCoord=floor(maxPos/size(croppedDataCube,1))+dataStruct(n).searchRangeY(1);
            if rem(maxPos,size(croppedDataCube,1))==0
                maxPos=size(croppedDataCube,1);
                dataStruct(n).maxYCoord=dataStruct(n).maxYCoord-1;
            else
                maxPos=rem(maxPos,size(croppedDataCube,1));
            end
            dataStruct(n).maxXCoord=maxPos-1+dataStruct(n).searchRangeX(1);
%             [dataStruct(n).maxXCoord,dataStruct(n).maxYCoord,dataStruct(n).maxZCoord]
            
            dataStruct(n).analysisRangeX=[max(1,dataStruct(n).maxXCoord-10),min(length(xRange),dataStruct(n).maxXCoord+10)];
            dataStruct(n).analysisRangeY=[max(1,dataStruct(n).maxYCoord-10),min(length(yRange),dataStruct(n).maxYCoord+10)];
            dataStruct(n).analysisRangeZ=[max(1,dataStruct(n).maxZCoord-20),min(length(zRange),dataStruct(n).maxZCoord+5)];
%             croppedDataCube=recordedImageStack(dataStruct(n).analysisRangeX(1):dataStruct(n).analysisRangeX(2),dataStruct(n).analysisRangeY(1):dataStruct(n).analysisRangeY(2),dataStruct(n).analysisRangeZ(1):dataStruct(n).analysisRangeZ(2));
%             % Plots for testing
%             figure(2);
%             subplot(3,1,1);imagesc(dataStruct(n).analysisRangeY,dataStruct(n).analysisRangeZ,squeeze(max(croppedDataCube,[],1))');axis image;
%             xlabel('y [pixels]');ylabel('z [pixels]');title(strcat('yz-proj - dataPoint #:',num2str(n)));
%             subplot(3,1,2);imagesc(dataStruct(n).analysisRangeX,dataStruct(n).analysisRangeZ,squeeze(max(croppedDataCube,[],2))');axis image;
%             xlabel('x [pixels]');ylabel('z [pixels]');title(strcat('xz-proj - dataPoint #:',num2str(n)));
%             subplot(3,1,3);imagesc(dataStruct(n).analysisRangeX,dataStruct(n).analysisRangeY,squeeze(max(croppedDataCube,[],3))');axis image;
%             xlabel('x [pixels]');ylabel('y [pixels]');title(strcat('xy-proj - dataPoint #:',num2str(n)));
%             %
        end
        [n]
    end
    
    for n=1:noDataPoints
        if dataStruct(n).Flag_noBead~=1;
            % Plot CoM and Max value in each z-plane to determine PSF tilt
            
            croppedDataCube=recordedImageStack(dataStruct(n).analysisRangeX(1):dataStruct(n).analysisRangeX(2),dataStruct(n).analysisRangeY(1):dataStruct(n).analysisRangeY(2),dataStruct(n).analysisRangeZ(1):dataStruct(n).analysisRangeZ(2));
            % Plots for testing
            figure(2);
            subplot(3,1,1);imagesc(squeeze(max(croppedDataCube,[],1))');axis image;
            xlabel('y [pixels]');ylabel('z [pixels]');title(strcat('yz-proj - dataPoint #:',num2str(n)));
            subplot(3,1,2);imagesc(squeeze(max(croppedDataCube,[],2))');axis image;
            xlabel('x [pixels]');ylabel('z [pixels]');title(strcat('xz-proj - dataPoint #:',num2str(n)));
            subplot(3,1,3);imagesc(squeeze(max(croppedDataCube,[],3))');axis image;
            xlabel('x [pixels]');ylabel('y [pixels]');title(strcat('xy-proj - dataPoint #:',num2str(n)));
            %
            
            zVector=[1:size(croppedDataCube,3)];
            for IdZ=1:size(croppedDataCube,3)
                imagePlane=croppedDataCube(:,:,IdZ);
                meanVal=mean(imagePlane(:));
                [maxVal,maxPos]=max(imagePlane(:));
                if maxVal<2*meanVal
                    maxYPos(IdZ)=NaN;
                    maxXPos(IdZ)=NaN;
                    zVector(IdZ)=NaN;
                else
                    maxYPos(IdZ)=floor(maxPos/size(imagePlane,1))+1;
                    if rem(maxPos,size(imagePlane,1))==0
                        maxPos=size(croppedDataCube,1);
                        maxYPos(IdZ)=maxYPos(IdZ)-1;
                    else
                        maxPos=rem(maxPos,size(imagePlane,1));
                    end
                    maxXPos(IdZ)=maxPos;             
                end
%                 [maxXPos(IdZ),maxYPos(IdZ)]
                figure(3);
                subplot(2,1,1);imagesc(imagePlane');axis image;
                xlabel('x [pixels]');ylabel('y [pixels]');title(strcat('Image - dataPoint #:',num2str(n),', zPlane=',num2str(IdZ)));
                guessImage=zeros(size(imagePlane'));
                if isnan(maxXPos(IdZ))==0
                    guessImage(maxXPos(IdZ),maxYPos(IdZ))=1;
                end
                subplot(2,1,2);imagesc(guessImage');axis image;
                xlabel('x [pixels]');ylabel('y [pixels]');title(strcat('Max Pixel - dataPoint #:',num2str(n),', zPlane=',num2str(IdZ)));
            end
            
            for IdZ=1:size(croppedDataCube,3)
                if isnan(maxXPos(IdZ))
                    maxXPos(IdZ)=0;
                    maxYPos(IdZ)=0;
                    zVector(IdZ)=0;
                end
            end
            maxXPos=nonzeros(maxXPos);
            maxYPos=nonzeros(maxYPos);
            zVector=nonzeros(zVector);
            
            pX=polyfit(zVector*zSpacing,maxXPos*xSpacing,1);
            pY=polyfit(zVector*zSpacing,maxYPos*ySpacing,1);
            fitX=polyval(pX,zVector*zSpacing);
            fitY=polyval(pY,zVector*zSpacing);
            figure(4);
            subplot(2,1,1);plot(maxXPos*xSpacing,zVector*zSpacing,fitX,zVector*zSpacing);
            xlabel('x [um]');ylabel('z [um]');title('x-z gradient');
            subplot(2,1,2);plot(maxYPos*ySpacing,zVector*zSpacing,fitY,zVector*zSpacing);
            xlabel('y [um]');ylabel('z [um]');title('y-z gradient');
            
            dataStruct(n).psfXZTilt=pX(1);
            dataStruct(n).psfYZTilt=pY(1);
        
        elseif dataStruct(n).Flag_noBead==1;
            dataStruct(n).psfXZTilt=NaN;
            dataStruct(n).psfYZTilt=NaN;
        end
        clear zVector maxXPos maxYPos
    end

    for m=1:length(dataStruct)
        if dataStruct(m).Flag_noBead==0
            xPosition(m)=dataStruct(m).maxXCoord;
            yPosition(m)=dataStruct(m).maxYCoord;
            xGradient(m)=dataStruct(m).psfXZTilt;
            yGradient(m)=dataStruct(m).psfYZTilt;
        else
            xPosition(m)=NaN;
            yPosition(m)=NaN;
            xGradient(m)=NaN;
            yGradient(m)=NaN;
        end
    end

    xPosition=xPosition(isnan(xPosition)~=1);
    yPosition=yPosition(isnan(yPosition)~=1);
    xGradient=xGradient(isnan(xGradient)~=1);
    yGradient=yGradient(isnan(yGradient)~=1);
    
    pXGrad=polyfit((xPosition-floor(length(xRange)/2))*xSpacing,xGradient,1);
    pYGrad=polyfit((yPosition-floor(length(yRange)/2))*ySpacing,yGradient,1);
    XGradfit=polyval(pXGrad,(xPosition-floor(length(xRange)/2))*xSpacing);
    YGradfit=polyval(pYGrad,(yPosition-floor(length(yRange)/2))*ySpacing);
    
    scanShear=[pYGrad(2) pXGrad(2)]*-1
    scaling=[pYGrad(1) pXGrad(1)]*-1
    
    figure(5);scatter((xPosition-floor(length(xRange)/2))*xSpacing,xGradient);
    xlabel('x [um]');ylabel('tilt angle [rad]');
    hold on;
    plot((xPosition-floor(length(xRange)/2))*xSpacing,XGradfit,'r');
    hold off;
    figure(6);scatter((yPosition-floor(length(yRange)/2))*ySpacing,yGradient);
    xlabel('y [um]');ylabel('tilt angle [rad]');
    hold on;
    plot((yPosition-floor(length(yRange)/2))*ySpacing,YGradfit,'r');
    hold off;
    
end