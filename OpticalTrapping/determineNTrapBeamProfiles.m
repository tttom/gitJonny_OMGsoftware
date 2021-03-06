%%% Function:           determineNTrapBeamProfiles
%%% Author:             Jonathan Nylk
%%% Created:            20/07/2015
%%% Description:        This function determines x-, y-, and z-axis beam
%%%                     profiles through the maximum pixels in a z-stack of
%%%                     images acquired of a distribution of N optical
%%%                     traps.
%%%                     The approximate x-y plane centre coordinates of
%%%                     each trap must be manually determined by clicking
%%%                     on an image which is a z-axis maximum intensity
%%%                     projection of the data (which will make it
%%%                     difficult for this code to handle traps which are
%%%                     very close in x and y but offset in z).
%%%                     This code is a generalisation of
%%%                     "determineTrapBeamProfile.m" which this function
%%%                     calls individaully for each trap.
%%%                     The z-stack data can either be input as a 3d array
%%%                     or a details for a file structure leading to the
%%%                     image stack.
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function []=determineNTrapBeamProfiles(dataMatrix,imageFolder,baseFileName,fileType,magnification,zStepSize,sizeROI,noTraps)

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
    if nargin<7
        sizeROI=100; %half-length of ROI box
    end
    if nargin<8
        noTraps=1;
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
    
    %create real and pixel based zero-centred coordinate system
    yRange_pix=(1:size(dataMatrix,1))-floor(size(dataMatrix,1)/2);
    yRange=yRange_pix*magnification;
    xRange_pix=(1:size(dataMatrix,2))-floor(size(dataMatrix,2)/2);
    xRange=xRange_pix*magnification;
    zRange_pix=(1:size(dataMatrix,3))-floor(size(dataMatrix,3)/2);
    zRange=zRange_pix*zStepSize;
    
    %generate user input beam centres for each trap
    for trapIdx=1:noTraps
        figure(99);
        imagesc(xRange_pix,yRange_pix,squeeze(max(dataMatrix,[],3)));
        axis image;
        xlabel('x [pixels]');ylabel('y [pixels]');
        title(strcat('Click on the centre of beam (',num2str(trapIdx),') of (',num2str(noTraps),')'));
%         drawnow;shg;
        [xClick,yClick]=ginput(1);
        yClick=yClick+floor(size(dataMatrix,1)/2);
        xClick=xClick+floor(size(dataMatrix,2)/2);
        close 99
    end
    
    %for each trap, crop dataMatrix around each user input centre and run
    %"determineTrapBeamProfile.m"
    for trapIdx=1:noTraps
        %crop ROI for current trap
        croppedData=dataMatrix(round(yClick(trapIdx))-sizeROI:round(yClick(trapIdx))+sizeROI...
            ,round(xClick(trapIdx))-sizeROI:round(xClick(trapIdx))+sizeROI,:);
        %get maximum pixel value in ROI as true centre
        [~,maxPos]=max(croppedData(:));
        maxCoords=zeros(1,3);
        [maxCoords(1),maxCoords(2),maxCoords(3)]=ind2sub(size(croppedData),maxPos);
        %determine maximum z-length centred on true centre
        croppedZHalfLength=min(maxCoords(3)-1,length(zRange)-maxCoords(3)-1);
        %re-crop ROI around true centre
        clear croppedData
        croppedData=dataMatrix((round(yClick(trapIdx))-sizeROI:round(yClick(trapIdx))+sizeROI)-sizeROI-1+maxCoords(1)...
            ,(round(xClick(trapIdx))-sizeROI:round(xClick(trapIdx))+sizeROI)-sizeROI-1+maxCoords(2)...
            ,maxCoords(3)-croppedZHalfLength:maxCoords(3)+croppedZHalfLength);
        
        [xWidth,yWidth,zWidth]=determineTrapBeamProfile(croppedData,'','','',magnification,zStepSize);

        disp(strcat('Widths of Beam (',num2str(trapIdx),'):'))
        disp(strcat('x width = ',num2str(xWidth)))
        disp(strcat('y width = ',num2str(yWidth)))
        disp(strcat('z width = ',num2str(zWidth)))




    end
    
end