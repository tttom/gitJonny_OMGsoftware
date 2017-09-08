%%% cross-correlation analysis of SIM nephrotic disease data (very much in
%%% development)

function [healthyData,mcdData] = nephroticCrossCorrelationAnalysis(healthyDir,mcdDir)

    if nargin < 1
        healthyDir = 'C:\Users\Jonathan Nylk\Desktop\Podocyte Images\Healthy';
    end

    if nargin < 2
        mcdDir = 'C:\Users\Jonathan Nylk\Desktop\Podocyte Images\MCD';
    end

    codeDir = pwd;
    
    
    %initialise output structures
    healthyData = [];
    healthyData.baseImage = [];
    healthyData.xCorrImage = [];
    healthyData.baseImageMeanIntensity = [];
    healthyData.xCorrBaseImageMeanIntensity = [];
    healthyData.xCorrMax = [];
    
    mcdData = [];
    mcdData.baseImage = [];
    mcdData.xCorrImage = [];
    mcdData.baseImageMeanIntensity = [];
    mcdData.xCorrBaseImageMeanIntensity = [];
    mcdData.xCorrMax = [];
    
    % load healthy data
    cd(healthyDir);
    fileList = dir('*.png');
    for fileIdx = 1:length(fileList)
        healthyData(fileIdx).baseImage = single(imread(fileList(fileIdx).name)) / ((2^16) - 1);
        healthyData(fileIdx).baseImageMeanIntensity = mean(healthyData(fileIdx).baseImage(:));
    end
    clear fileList;
    
    % load mcd data
    cd(mcdDir);
    fileList = dir('*.png');
    for fileIdx = 1:length(fileList)
        mcdData(fileIdx).baseImage = single(imread(fileList(fileIdx).name)) / ((2^16) - 1);
        mcdData(fileIdx).baseImageMeanIntensity = mean(mcdData(fileIdx).baseImage(:));
    end
    clear fileList
    
    cd(codeDir);
    
    %determine cross-correlation between datasets
    tic
    [outputDataStructure] = computeXCorrData(healthyData,healthyData);
    toc
    
    figure;
    for n = 1:length(outputDataStructure)
        scatter(outputDataStructure(n).baseImageMeanIntensity .* outputDataStructure(n).xCorrBaseImageMeanIntensity,outputDataStructure(n).xCorrMax);
        hold on
    end
    xlabel('product of mean image intensities');
    ylabel('max. xcorr');
    title('xcorr v. mean image brightness');
    
    figure;
    randBaseImageNo = randi(length(outputDataStructure),[1,1]);
    subplot(5,6,1);
    imagesc(outputDataStructure(randBaseImageNo).baseImage);axis image;
    title(strcat('base image:',num2str(randBaseImageNo)));
    for n = 2:30
        randProbeImageNo = randi(length(outputDataStructure(72).xCorrMax),[1,1]);
        subplot(5,6,n);
        imagesc(squeeze(outputDataStructure(randBaseImageNo).xCorrImage(:,:,randProbeImageNo)));axis image;
        title(num2str(randProbeImageNo));
    end
    clear outputDataStructure;
    
    
    tic
    [outputDataStructure] = computeXCorrData(mcdData,mcdData);
    toc
    
    figure;
    for n = 1:length(outputDataStructure)
        scatter(outputDataStructure(n).baseImageMeanIntensity .* outputDataStructure(n).xCorrBaseImageMeanIntensity,outputDataStructure(n).xCorrMax);
        hold on
    end
    xlabel('product of mean image intensities');
    ylabel('max. xcorr');
    title('xcorr v. mean image brightness');
    
    figure;
    randBaseImageNo = randi(length(outputDataStructure),[1,1]);
    subplot(5,6,1);
    imagesc(outputDataStructure(randBaseImageNo).baseImage);axis image;
    title(strcat('base image:',num2str(randBaseImageNo)));
    for n = 2:30
        randProbeImageNo = randi(length(outputDataStructure(72).xCorrMax),[1,1]);
        subplot(5,6,n);
        imagesc(squeeze(outputDataStructure(randBaseImageNo).xCorrImage(:,:,randProbeImageNo)));axis image;
        title(num2str(randProbeImageNo));
    end
    clear outputDataStructure;
    
    
    tic
    [outputDataStructure] = computeXCorrData(healthyData,mcdData);
    toc
    
    figure;
    for n = 1:length(outputDataStructure)
        scatter(outputDataStructure(n).baseImageMeanIntensity .* outputDataStructure(n).xCorrBaseImageMeanIntensity,outputDataStructure(n).xCorrMax);
        hold on
    end
    xlabel('product of mean image intensities');
    ylabel('max. xcorr');
    title('xcorr v. mean image brightness');
    
    figure;
    randBaseImageNo = randi(length(outputDataStructure),[1,1]);
    subplot(5,6,1);
    imagesc(outputDataStructure(randBaseImageNo).baseImage);axis image;
    title(strcat('base image:',num2str(randBaseImageNo)));
    for n = 2:30
        randProbeImageNo = randi(length(outputDataStructure(72).xCorrMax),[1,1]);
        subplot(5,6,n);
        imagesc(squeeze(outputDataStructure(randBaseImageNo).xCorrImage(:,:,randProbeImageNo)));axis image;
        title(num2str(randProbeImageNo));
    end
    clear outputDataStructure;
    
    
    tic
    [outputDataStructure] = computeXCorrData(mcdData,healthyData);
    toc
    
    figure;
    for n = 1:length(outputDataStructure)
        scatter(outputDataStructure(n).baseImageMeanIntensity .* outputDataStructure(n).xCorrBaseImageMeanIntensity,outputDataStructure(n).xCorrMax);
        hold on
    end
    xlabel('product of mean image intensities');
    ylabel('max. xcorr');
    title('xcorr v. mean image brightness');
    
    figure;
    randBaseImageNo = randi(length(outputDataStructure),[1,1]);
    subplot(5,6,1);
    imagesc(outputDataStructure(randBaseImageNo).baseImage);axis image;
    title(strcat('base image:',num2str(randBaseImageNo)));
    for n = 2:30
        randProbeImageNo = randi(length(outputDataStructure(72).xCorrMax),[1,1]);
        subplot(5,6,n);
        imagesc(squeeze(outputDataStructure(randBaseImageNo).xCorrImage(:,:,randProbeImageNo)));axis image;
        title(num2str(randProbeImageNo));
    end
    clear outputDataStructure;
    
    

end




function [outputDataStructure] = computeXCorrData(inputdataStructure1,inputdataStructure2)

    for fileIdx = 1:length(inputdataStructure1)
        for xIdx = 1:length(inputdataStructure2)
            [fileIdx,xIdx]
            inputdataStructure1(fileIdx).xCorrBaseImageMeanIntensity(xIdx) = inputdataStructure2(xIdx).baseImageMeanIntensity;
            
%             % gpu implementation (doesn't work)
%             GPU_input1 = gpuArray(inputdataStructure1(fileIdx).baseImage);
%             GPU_input2 = gpuArray(inputdataStructure2(xIdx).baseImage);
%             GPU_output = xcorr2(GPU_input1,GPU_input2);
%             GPU_output = gather(GPU_output);
%             inputdataStructure1(fileIdx).xCorrImage...
%                 = cat(3,inputdataStructure1(fileIdx).xCorrImage,GPU_output);
            
%             inputdataStructure1(fileIdx).xCorrImage...
%                 = cat(3,inputdataStructure1(fileIdx).xCorrImage,xcorr2(inputdataStructure1(fileIdx).baseImage,inputdataStructure2(xIdx).baseImage));

            inputdataStructure1(fileIdx).xCorrImage...
                = cat(3,inputdataStructure1(fileIdx).xCorrImage,normxcorr2(inputdataStructure1(fileIdx).baseImage,inputdataStructure2(xIdx).baseImage));
            inputdataStructure1(fileIdx).xCorrMax(xIdx) = max(max(inputdataStructure1(fileIdx).xCorrImage(:,:,xIdx)));
        end
    end
    
    outputDataStructure = inputdataStructure1;
    clear inputdataStructure1 inputdataStructure2;

end