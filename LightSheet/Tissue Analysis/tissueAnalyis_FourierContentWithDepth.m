%%% Function:           tissueAnalyis_FourierContentWithDepth
%%% Author:             Jonathan Nylk
%%% Created:            19/20/2016
%%% Description:        
%%%
%%% Inputs:             fileNames:  Cell containing all datasets of same
%%%                                 beam type and lateral tissue location
%%%                     
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function tissueAnalyis_FourierContentWithDepth(fileNames) 
    
    %variables:
    normalisation_FLAG=1; %normalise input data to [0 1] range if true
    lambda=532e-9;
    numericalAperture=0.42; %these can be pulled from the datafiles in future
    fftThreshold=0.05; % make function input

    if nargin<1
        % files order with surface first, then increasing depth
        fileNames={'./TestData/Airy_01.tif','./TestData/Airy_02.tif','./TestData/Airy_03.tif'};
    end
    
    % ensure input fileNames is in cell format
    if exist('fileNames','var')
        if ischar(fileNames)
            fileNames={fileNames};
        end
    end
    
    % allocate memory for results
    projections=zeros([2000,2048,length(fileNames)],'single');
    xStepSize=zeros([1,length(fileNames)],'single');
    zStepSize=zeros([1,length(fileNames)],'single');
    
    % load test data
    for n=1:length(fileNames)
        importedImage=single(imread(fileNames{n}));
        xRange(:,n)=[1:size(importedImage,2)]*0.08*1e-6;
        zRange(:,n)=[1:size(importedImage,1)]*0.08*1e-6;
        projections(:,:,n)=importedImage;
        clear importedImage;
        if normalisation_FLAG==1
            projections(:,:,n)=projections(:,:,n)/max(max(projections(:,:,n)));
        end
    end
    
    % allocate more memory for results
    rotProjections=zeros(size(projections));
    rotZRange=zeros(size(zRange));
    rotXRange=zeros(size(xRange));
    relativeCoords=zeros(2,length(fileNames));
    absRotZRange=zeros(size(rotZRange));
    k_rotXRange=zeros(size(rotXRange));
    fftProjections=zeros(size(rotProjections));
    thresholdedFftProjections=zeros(size(rotProjections));
    
    % create zero centred coordinate system for data regardless of actual
    % coordinate system
    for n=1:length(fileNames)
        xStepSize(n)=xRange(2,n)-xRange(1,n);
        xRange(:,n)=([1:size(projections(:,:,n),2)]-floor(size(projections(:,:,n),2)/2))*xStepSize(n);
        zStepSize(n)=zRange(2,n)-zRange(1,n);
        zRange(:,n)=([1:size(projections(:,:,n),1)]-floor(size(projections(:,:,n),1)/2))*zStepSize(n);
    end
    
    % rotate projections and coordinates by 45deg clockwise so tissue surface is
    % horizontal
    for n=1:length(fileNames)
        [rotProjections(:,:,n),rotZRange(:,n),rotXRange(:,n)]=rotate2DArray(projections(:,:,n),45/360*2*pi,zRange(:,n),xRange(:,n),0,0);
        rotProjections(:,:,n)=rotProjections(:,:,n)/max(max(rotProjections(:,:,n)));
    end
    
    % view top layer to select surface, then identify common features in
    % each consecutive pair of datasets and determine absolute tissue depth
    for n=1:length(fileNames)
        if n==1
            figure(1);imagesc(rotXRange(:,n)*1e6,rotZRange(:,n)*1e6,rotProjections(:,:,n));axis image;
            title('Surface layer, - Click point on surface');
            colormap hsv;
            drawnow;shg;
            [~,yClick]=ginput(1);
            relativeCoords(:,1)=[0 yClick*1e-6];
        else
            leftImage=subplot(1,2,1);imagesc(rotXRange(:,n-1)*1e6,rotZRange(:,n-1)*1e6,rotProjections(:,:,n-1));axis image;
            title(strcat('Layer (',num2str(n-1),') - click a point in this image'));
            colormap hsv;
            rightImage=subplot(1,2,2);imagesc(rotXRange(:,n)*1e6,rotZRange(:,n)*1e6,rotProjections(:,:,n));axis image;
            title(strcat('Layer (',num2str(n),')'));
            colormap hsv;
            drawnow;shg;
            axis(leftImage);
            [~,yClick]=ginput(1);
            relativeCoords(1,n)=yClick*1e-6;
            leftImage=subplot(1,2,1);imagesc(rotXRange(:,n-1)*1e6,rotZRange(:,n-1)*1e6,rotProjections(:,:,n-1));axis image;
            title(strcat('Layer (',num2str(n-1),')'));
            colormap hsv;
            rightImage=subplot(1,2,2);imagesc(rotXRange(:,n)*1e6,rotZRange(:,n)*1e6,rotProjections(:,:,n));axis image;
            title(strcat('Layer (',num2str(n),') - now click the same point in this image'));
            colormap hsv;
            drawnow;shg;
            axis(rightImage);
            [~,yClick]=ginput(1);
            relativeCoords(2,n)=yClick*1e-6;
        end
    end
    close(1);
    
    figure();
    for n=1:length(fileNames)
        % convert to absolute tissue depth coordinate
        absRotZRange(:,n)=rotZRange(:,n);
        for m=1:n
            absRotZRange(:,n)=absRotZRange(:,n)-relativeCoords(2,m)+relativeCoords(1,m);
        end
        subplot(1,3,n);imagesc(rotXRange(:,n)*1e6,absRotZRange(:,n)*1e6,rotProjections(:,:,n));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
    end
    
    figure();
    % Fourier transform at different tissue depths
    for n=1:length(fileNames)
        % Fourier coordinate system
%         k_rotXRange(:,n)=([1:length(rotXRange(:,n))]-floor(length(rotXRange(:,n))/2))/(length(rotXRange(:,n))*xStepSize(n)/2);
        k_rotXRange(:,n)=[1:length(rotXRange(:,n))]/length(rotXRange(:,n))/xStepSize(n)*2;
        k_rotXRange(:,n)=k_rotXRange(:,n)-k_rotXRange(1,n);
        k_rotXRange(:,n)=k_rotXRange(:,n)*lambda/2/numericalAperture;
        
        % FFT projections
        for m=1:size(rotProjections(:,1,n))
%             fftProjections(m,:,n)=fftshift(fft(squeeze(rotProjections(m,:,n))));
            fftProjections(m,:,n)=fft(squeeze(rotProjections(m,:,n)));
        end
        fftProjections=abs(fftProjections);
        fftProjections=fftProjections/max(fftProjections(:));
        thresholdedFftProjections=fftProjections>fftThreshold;
        subplot(2,3,n);imagesc(k_rotXRange(:,n),absRotZRange(:,n)*1e6,fftProjections(:,:,n));
        xlabel('k_{x+z}');ylabel('Tissue depth (x-z) [um]');
        xlim([0 1]);
        subplot(2,3,3+n);imagesc(k_rotXRange(:,n),absRotZRange(:,n)*1e6,thresholdedFftProjections(:,:,n));
        xlabel('k_{x+z}');ylabel('Tissue depth (x-z) [um]');
        xlim([0 1]);
    end
    
    thresholdedFftProjections=thresholdedFftProjections(:,floor(end/2)+1:end,:);
    [~,thresholdFftLinesIndex]=max(thresholdedFftProjections,[],2);
    thresholdFftLinesIndex=squeeze(thresholdFftLinesIndex);
    for n=1:length(fileNames)
        for m=1:size(thresholdedFftProjections,1)
            thresholdFftLines(m,n)=k_rotXRange(floor(end/2)+thresholdFftLinesIndex(m,n));
        end
    end
    figure(4);plot(absRotZRange(:,1)*1e6,thresholdFftLines(:,1)...
        ,absRotZRange(:,2)*1e6,thresholdFftLines(:,2)...
        ,absRotZRange(:,3)*1e6,thresholdFftLines(:,3));
end