%%% Function:           determineTrapBeamProfile
%%% Author:             Jonathan Nylk
%%% Created:            21/07/2015
%%% Description:        This function determines x-, y-, and z-axis beam
%%%                     profile through the maximum pixels in a z-stack of
%%%                     images acquired of an optical trap.
%%%                     The z-stack data can either be input as a 3d array
%%%                     or a details for a file structure leading to the
%%%                     image stack.
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function []=determineTrapBeamProfile(dataMatrix,imageFolder,baseFileName,fileType,magnification,zStepSize)

    %default values
    if nargin<1
        dataMatrix=[];
    end
    if nargin<2
        imageFolder='X:\SIM&Trapping\TrapBeamProfile\TrapProfileImages';
    end
    if nargin<3
        baseFileName='TrapBeamProfile_800nmGRINlens_croppedz';
    end
    if nargin<4
        fileType='.tif';
    end
    if nargin<5
        magnification=0.11; % um/pix
    end
    if nargin<6
        zStepSize=0.8; %um
    end


     %check for an input dataset
    if isempty(dataMatrix)
        %dataMatrix is empty, load data from file
        %generateImageFileList
        fullImageFileList=dir(strcat(imageFolder,'/',baseFileName,'*',fileType));
        testImage=single(imread(strcat(imageFolder,'/',fullImageFileList(1).name)));
        imHeight=size(testImage,1);
        imWidth=size(testImage,2);
        clear testImage
        %allocate memory for dataMatrix
        dataMatrix=zeros(imHeight,imWidth,length(fullImageFileList),'single');
        %load images
        for fileIdx=1:length(fullImageFileList)
            dataMatrix(:,:,fileIdx)=single(imread(strcat(imageFolder,'/',fullImageFileList(fileIdx).name)));
        end
    else %dataMatrix is not empty, use dataMatrix
        dataMatrix=single(dataMatrix);
    end
    
    %get maximum pixel value in ROI as true centre
    [~,maxPos]=max(dataMatrix(:));
    maxCoords=zeros(1,3);
    [maxCoords(1),maxCoords(2),maxCoords(3)]=ind2sub(size(dataMatrix),maxPos);
    %determine maximum array size centred on true centre
    croppedYHalfLength=min(maxCoords(1)-1,size(dataMatrix,1)-maxCoords(1)-1);
    croppedXHalfLength=min(maxCoords(2)-1,size(dataMatrix,2)-maxCoords(2)-1);
    croppedZHalfLength=min(maxCoords(3)-1,size(dataMatrix,3)-maxCoords(3)-1);
    %crop ROI around true centre if required
    if 2*croppedYHalfLength+1~=size(dataMatrix,1) && 2*croppedXHalfLength+1~=size(dataMatrix,2) && 2*croppedZHalfLength+1~=size(dataMatrix,3)
        croppedData=dataMatrix(maxCoords(1)-croppedYHalfLength:maxCoords(1)+croppedYHalfLength...
            ,maxCoords(2)-croppedXHalfLength:maxCoords(2)+croppedXHalfLength...
            ,maxCoords(3)-croppedZHalfLength:maxCoords(3)+croppedZHalfLength);
    else
        croppedData=dataMatrix;
    end
    clear dataMatrix
    %create real and pixel based zero-centred coordinate system
    yRange_pix=(1:size(croppedData,1))-floor(size(croppedData,1)/2);
    yRange=yRange_pix*magnification;
    xRange_pix=(1:size(croppedData,2))-floor(size(croppedData,2)/2);
    xRange=xRange_pix*magnification;
    zRange_pix=(1:size(croppedData,3))-floor(size(croppedData,3)/2);
    zRange=zRange_pix*zStepSize;
    
    %determine tilt of beam (in pixels) using centre-of-mass in each plane
    xCentre=zeros(1,length(zRange_pix));
    yCentre=zeros(1,length(zRange_pix));
    for zIdx=1:length(zRange_pix)
        [xCentre(zIdx),yCentre(zIdx)]=determine2DCentreOfMass(squeeze(croppedData(:,:,zIdx)),xRange_pix,yRange_pix);
    end
    
    %fit line to CoM traces
    fitRange=1:length(zRange_pix);
    polyX=polyfit(zRange_pix(fitRange),xCentre(fitRange),1);
    polyY=polyfit(zRange_pix(fitRange),yCentre(fitRange),1);
    fitX=polyval(polyX,zRange_pix);
    fitY=polyval(polyY,zRange_pix);
            
    %rotate datacube to counteract tilt
    angleX=atan(polyX(1));
    angleY=atan(polyY(1));
%     angleX=45/360*2*pi;
%     angleY=0/360*2*pi;
    workingArray=singleAxisRotate3DArray(croppedData,1,angleX);
    workingArray=singleAxisRotate3DArray(workingArray,2,angleY);
    %determine new "unit vectors" for x,y, and z axes
    rotXMag=magnification*cos(angleX)+zStepSize*sin(angleX);
    rotYMag=magnification*cos(angleY)+zStepSize*sin(angleY);
    rotZStepSize=-1*magnification*sin(angleX)+zStepSize*cos(angleX);
    rotZStepSize=-1*magnification*sin(angleY)+rotZStepSize*cos(angleY);
    
    %rescale rotated coordinates
    yRange=yRange_pix*rotXMag;
    xRange=xRange_pix*rotYMag;
    zRange=zRange_pix*rotZStepSize;

    %plot beam profiles through centre
    figure()
    subplot(1,2,1);
    imagesc(xRange,zRange,squeeze(workingArray(floor(length(yRange)/2),:,:)).');axis image;
    xlabel('x [um]');ylabel('z [um]');
    title('x-z cross-section');
    subplot(1,2,2);
    imagesc(yRange,zRange,squeeze(workingArray(:,floor(length(xRange)/2),:)).');axis image;
    xlabel('y [um]');ylabel('z [um]');
    title('y-z cross-section');
    
    %plot line profiles through centre
    figure()
    subplot(1,2,1);
    plot(xRange,workingArray(floor(length(yRange)/2),:,floor(length(zRange)/2))...
        ,yRange,workingArray(:,floor(length(xRange)/2),floor(length(zRange)/2)));
    xlabel('y [um]');ylabel('Intensity [a.u.]');
    title('x-, and y-axis line profiles');
    xlim([-5 5]);
    subplot(1,2,2);
    plot(zRange,squeeze(workingArray(floor(length(yRange)/2),floor(length(xRange)/2),:)));
    xlabel('z [um]');ylabel('Intensity [a.u.]');
    title('z-axis line profile');
    xlim([zRange(1) zRange(end)]);
    
    %add Gaussian fitting code to determine widths of peaks
    %output x-, y-, and z-axis widths
    
end