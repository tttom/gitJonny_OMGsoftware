%%% texture analysis of SIM nephrotic disease data (very much in
%%% development)

function [healthyData,mcdData] = nephroticTextureAnalysis(healthyDir,mcdDir)

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
    healthyData.glcms = [];
    healthyData.SI = [];
    healthyData.stats = [];
    healthyData.stats.Contrast = [];
    healthyData.stats.Correlation = [];
    healthyData.stats.Energy = [];
    healthyData.stats.Homogeneity = [];
    
    mcdData = [];
    mcdData.baseImage = [];
    mcdData.glcms = [];
    mcdData.SI = [];
    mcdData.stats = [];
    mcdData.stats.Contrast = [];
    mcdData.stats.Correlation = [];
    mcdData.stats.Energy = [];
    mcdData.stats.Homogeneity = [];
    
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
    
    %texture analysis
    for n = 1:length(healthyData)
        [healthyData(n).glcms,healthyData(n).SI]...
            = graycomatrix(healthyData(n).baseImage,'NumLevels',8,'Offset',[0 1; -1 1;-1 0;-1 -1],'Symmetric',true);
%         [healthyData(n).glcms,healthyData(n).SI]...
%             = graycomatrix(healthyData(n).baseImage,'NumLevels',8,'GrayLimits',[0 1],'Offset',[0 1; -1 1;-1 0;-1 -1],'Symmetric',true);
        stats = graycoprops(healthyData(n).glcms);
        healthyData(n).stats.Contrast = stats.Contrast;
        healthyData(n).stats.Correlation = stats.Correlation;
        healthyData(n).stats.Energy = stats.Energy;
        healthyData(n).stats.Homogeneity = stats.Homogeneity;
    end
    
    for n = 1:length(mcdData)
        [mcdData(n).glcms,mcdData(n).SI]...
            = graycomatrix(mcdData(n).baseImage,'NumLevels',8,'Offset',[0 1; -1 1;-1 0;-1 -1],'Symmetric',true);
%         [mcdData(n).glcms,mcdData(n).SI]...
%             = graycomatrix(mcdData(n).baseImage,'NumLevels',8,'GrayLimits',[0 1],'Offset',[0 1; -1 1;-1 0;-1 -1],'Symmetric',true);
        stats = graycoprops(mcdData(n).glcms);
        mcdData(n).stats.Contrast = stats.Contrast;
        mcdData(n).stats.Correlation = stats.Correlation;
        mcdData(n).stats.Energy = stats.Energy;
        mcdData(n).stats.Homogeneity = stats.Homogeneity;
    end
    
%     figure;
%     for n = 1:length(healthyData)
%         subplot(3,1,1);
%         imagesc(healthyData(n).baseImage);axis image;
%         subplot(3,1,2);
%         imagesc(squeeze(healthyData(n).SI(:,:,1)));axis image;
%         subplot(3,1,3);
%         imagesc(squeeze(healthyData(n).glcms(:,:,1)));axis image;
%         colormap gray;
%         drawnow;shg;
%         pause(0.1);
%     end
%     
%     figure;
%     for n = 1:length(mcdData)
%         subplot(3,1,1);
%         imagesc(mcdData(n).baseImage);axis image;
%         subplot(3,1,2);
%         imagesc(squeeze(mcdData(n).SI(:,:,1)));axis image;
%         subplot(3,1,3);
%         imagesc(squeeze(mcdData(n).glcms(:,:,1)));axis image;
%         colormap gray;
%         drawnow;shg;
%         pause(0.1);
%     end
    
    for n = 1:length(healthyData)
        healthyContrast(n,:) = healthyData(n).stats.Contrast;
        healthyCorrelation(n,:) = healthyData(n).stats.Correlation;
        healthyEnergy(n,:) = healthyData(n).stats.Energy;
        healthyHomogeneity(n,:) = healthyData(n).stats.Homogeneity;
    end

    figure;
    subplot(4,1,1);
    plot(1:length(healthyData),healthyContrast(:,1),1:length(healthyData),healthyContrast(:,2)...
        ,1:length(healthyData),healthyContrast(:,3),1:length(healthyData),healthyContrast(:,4));
    title('Contrast');
    subplot(4,1,2);
    plot(1:length(healthyData),healthyCorrelation(:,1),1:length(healthyData),healthyCorrelation(:,2)...
        ,1:length(healthyData),healthyCorrelation(:,3),1:length(healthyData),healthyCorrelation(:,4));
    title('Correlation');
    subplot(4,1,3);
    plot(1:length(healthyData),healthyEnergy(:,1),1:length(healthyData),healthyEnergy(:,2)...
        ,1:length(healthyData),healthyEnergy(:,3),1:length(healthyData),healthyEnergy(:,4));
    title('Energy');
    subplot(4,1,4);
    plot(1:length(healthyData),healthyHomogeneity(:,1),1:length(healthyData),healthyHomogeneity(:,2)...
        ,1:length(healthyData),healthyHomogeneity(:,3),1:length(healthyData),healthyHomogeneity(:,4));
    title('Homogeneity');
    
    
    for n = 1:length(mcdData)
        mcdContrast(n,:) = mcdData(n).stats.Contrast;
        mcdCorrelation(n,:) = mcdData(n).stats.Correlation;
        mcdEnergy(n,:) = mcdData(n).stats.Energy;
        mcdHomogeneity(n,:) = mcdData(n).stats.Homogeneity;
    end

    figure;
    subplot(4,1,1);
    plot(1:length(mcdData),mcdContrast(:,1),1:length(mcdData),mcdContrast(:,2)...
        ,1:length(mcdData),mcdContrast(:,3),1:length(mcdData),mcdContrast(:,4));
    title('Contrast');
    subplot(4,1,2);
    plot(1:length(mcdData),mcdCorrelation(:,1),1:length(mcdData),mcdCorrelation(:,2)...
        ,1:length(mcdData),mcdCorrelation(:,3),1:length(mcdData),mcdCorrelation(:,4));
    title('Correlation');
    subplot(4,1,3);
    plot(1:length(mcdData),mcdEnergy(:,1),1:length(mcdData),mcdEnergy(:,2)...
        ,1:length(mcdData),mcdEnergy(:,3),1:length(mcdData),mcdEnergy(:,4));
    title('Energy');
    subplot(4,1,4);
    plot(1:length(mcdData),mcdHomogeneity(:,1),1:length(mcdData),mcdHomogeneity(:,2)...
        ,1:length(mcdData),mcdHomogeneity(:,3),1:length(mcdData),mcdHomogeneity(:,4));
    title('Homogeneity');
    
    %%% PCA
    
    PCAInput = vertcat(horzcat(healthyContrast,healthyCorrelation,healthyEnergy,healthyHomogeneity)...;
        ,horzcat(mcdContrast,mcdCorrelation,mcdEnergy,mcdHomogeneity));
    PCAInput(isnan(PCAInput)) = 1;
    
    diseaseState = horzcat(repmat('b',1,108),repmat('r',1,144));
    
    mdat = mean(PCAInput,2);
    vdat = PCAInput - repmat(mdat,1,size(PCAInput,2));
    mat = vdat' * vdat;
    [vec, ~] = eigs(mat,size(PCAInput,2));
    
    pc1 = vdat*vec(:,1);
    pc2 = vdat*vec(:,2);
    pc3 = vdat*vec(:,3);
    
    figure;
    for t = 1:length(diseaseState)
        subplot(1,3,1);
        plot(pc1(t,1),pc2(t,1),'*','color',diseaseState(t)); hold on;
        xlabel('PC1');ylabel('PC2');
        subplot(1,3,2);
        plot(pc1(t,1),pc3(t,1),'*','color',diseaseState(t)); hold on;
        xlabel('PC1');ylabel('PC3');
        subplot(1,3,3);
        plot(pc2(t,1),pc3(t,1),'*','color',diseaseState(t)); hold on;
        xlabel('PC2');ylabel('PC3');
    end
    
    figure;
    for t = 1:length(diseaseState)
        plot3(pc1(t,1),pc2(t,1),pc3(t,1),'*','color',diseaseState(t)); hold on;
        xlabel('PC2');ylabel('PC3');zlabel('PC3');
    end
    
end