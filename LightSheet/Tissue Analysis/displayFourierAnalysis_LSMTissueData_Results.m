%%% Function:           displayFourierAnalysis_LSMTissueData_Results.m
%%% Author:             Jonathan Nylk
%%% Created:            09/03/2016
%%% Description:        This function displays results from
%%%                     "FourierAnalysis_LSMTissueData.m" and determines an
%%%                     absolute depth axis for these reuslts.
%%%
%%%                     The datasets relating to a specific transverse
%%%                     tissue location (but different depth locations) are
%%%                     loaded and a GUI instructs the user how to specify
%%%                     the z-position of the tissue surface, and to slect
%%%                     common features in datasets from different depths
%%%                     to define absolute depth.
%%%
%%%                     The new z-axis is applied and the images and
%%%                     Fourier content plots displayed. Again the user is
%%%                     asked for inputs. The height and width of the
%%%                     apodized images where features are present in the
%%%                     images is required (the Fourier content graphs are
%%%                     meaningless otherwise). The graphs are then
%%%                     replotted over the user specified range.
%%%
%%% Inputs:             FolderPath: String containing folder path to
%%%                     results folder containing all results files for a
%%%                     given tissue region (all depths).
%%%                     
%%%
%%% Updates (latest first):
%%%         28/07/2016: Properly documented code.
%%%
%%%
%%% END

function displayFourierAnalysis_LSMTissueData_Results(FolderPath)

    % Default inputs
    if nargin<1
        FolderPath = 'Full file path to results folder';
    end

    fileList = dir(strcat(FolderPath,'\results_*.mat'));
    % Load results data
    for fileIdx = 1:length(fileList)
        fileName = strcat(FolderPath,'\',fileList(fileIdx).name);
        load(fileName);
        resultsAiryFile{fileIdx} = results.Airy.fileName;
        resultsGaussianFile{fileIdx} = results.Gaussian.fileName;
        fullAiryProjection(:,:,fileIdx) = results.Airy.fullProjection.';
        ROIAiryProjection(:,:,fileIdx) = results.Airy.ROIProjection.';
        fullGaussianProjection(:,:,fileIdx) = results.Gaussian.fullProjection.';
        ROIGaussianProjection(:,:,fileIdx) = results.Gaussian.ROIProjection.';
        rotYRange(:,fileIdx) = results.Airy.rotYRange;
        rotZRange(:,fileIdx) = results.Airy.rotZRange;
        AiryFourierContent(:,:,fileIdx) = results.Airy.fourierContent;
        GaussianFourierContent(:,:,fileIdx) = results.Gaussian.fourierContent;
        clear results;
    end
    
    % Display projections of data volumes and have user click to define top
    % surface and relate common features in each data volume for definition
    % of absolute tissue depth
    for fileIdx=1:length(fileList)
        % Display top surface and have user click on surface
        if fileIdx == 1
            figure(999);imagesc(rotYRange(:,fileIdx) * 1e6,rotZRange(:,fileIdx) * 1e6,fullAiryProjection(:,:,fileIdx));axis image;
            title('Surface layer, - Click point on surface');
            colormap hsv;   % increase contrast for ease of feature identification
            drawnow;shg;
            [~,yClick] = ginput(1);   % take y-coordinate (z-axis) position from user click
            relativeCoords(:,1) = [0 yClick * 1e-6];  % surface location
        else
            leftImage = subplot(1,2,1);imagesc(rotYRange(:,fileIdx - 1) * 1e6,rotZRange(:,fileIdx - 1) * 1e6...
                ,fullAiryProjection(:,:,fileIdx - 1));axis image;
            title(strcat('Layer (',num2str(fileIdx - 1),') - click a point in this image'));
            colormap hsv;
            rightImage = subplot(1,2,2);imagesc(rotYRange(:,fileIdx) * 1e6,rotZRange(:,fileIdx) * 1e6...
                ,fullAiryProjection(:,:,fileIdx));axis image;
            title(strcat('Layer (',num2str(fileIdx),')'));
            colormap hsv;
            drawnow;shg;
            axis(leftImage);
            [~,yClick] = ginput(1);
            relativeCoords(1,fileIdx) = yClick * 1e-6;  % fileIdx-th feature location
            leftImage = subplot(1,2,1);imagesc(rotYRange(:,fileIdx - 1) * 1e6,rotZRange(:,fileIdx - 1) * 1e6...
                ,fullAiryProjection(:,:,fileIdx - 1));axis image;
            title(strcat('Layer (',num2str(fileIdx - 1),')'));
            colormap hsv;
            rightImage = subplot(1,2,2);imagesc(rotYRange(:,fileIdx) * 1e6,rotZRange(:,fileIdx) * 1e6...
                ,fullAiryProjection(:,:,fileIdx));axis image;
            title(strcat('Layer (',num2str(fileIdx),') - now click the same point in this image'));
            colormap hsv;
            drawnow;shg;
            axis(rightImage);
            [~,yClick] = ginput(1);
            relativeCoords(2,fileIdx) = yClick * 1e-6;  % fileIdx-th feature location
        end
    end
    close(999);
    
    % Determine absolute z-axis range for each data volume
    for fileIdx=1:length(fileList)
        % Convert to absolute tissue depth coordinate
        absRotZRange(:,fileIdx) = rotZRange(:,fileIdx);
        for m = 1:fileIdx
            absRotZRange(:,fileIdx) = absRotZRange(:,fileIdx) - relativeCoords(2,m) + relativeCoords(1,m);
        end
    end
    
    % Display results
    for fileIdx=1:length(fileList)
        
        % Print file path for original datasets
        disp(strcat('Airy file:',resultsAiryFile{fileIdx}));
        disp(strcat('Gaussian file:',resultsGaussianFile{fileIdx}));
        
        % Display results over full data volume
        figure;
        % Full Gaussian projection
        subplot(5,4,[1 5]);imagesc(rotYRange(:,fileIdx) * 1e6...
            ,absRotZRange(:,fileIdx) * 1e6...
            ,fullGaussianProjection(:,:,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        title('Gaussian full image');
        % Full Airy projection
        subplot(5,4,[2 6]);imagesc(rotYRange(:,fileIdx) * 1e6...
            ,absRotZRange(:,fileIdx) * 1e6...
            ,fullAiryProjection(:,:,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        title('Airy full image');
        % Apodized Gaussian projection
        subplot(5,4,[3 7]);imagesc(rotYRange(:,fileIdx) * 1e6...
            ,absRotZRange(:,fileIdx) * 1e6...
            ,ROIGaussianProjection(:,:,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        title('Gaussian cropped image');
        % Apodized Airy projection
        subplot(5,4,[4 8]);imagesc(rotYRange(:,fileIdx) * 1e6...
            ,absRotZRange(:,fileIdx) * 1e6...
            ,ROIAiryProjection(:,:,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        title('Airy cropped image');
        % Plot Fourier content plots (left axis) and ratios (right axis)
        for fftIdx = 1:size(AiryFourierContent,1)
            subplot(5,4,8 + fftIdx);[dualAxes,plotLeft,plotright] = plotyy(absRotZRange(:,fileIdx) * 1e6...
                ,[squeeze(GaussianFourierContent(fftIdx,:,fileIdx));squeeze(AiryFourierContent(fftIdx,:,fileIdx))]...
                ,absRotZRange(:,fileIdx) * 1e6...
                ,squeeze(AiryFourierContent(fftIdx,:,fileIdx)) ./ squeeze(GaussianFourierContent(fftIdx,:,fileIdx)));
            axes(dualAxes(1));
            xlabel('Depth [um]');ylabel('Magnitude [a.u.]');
            title(strcat(num2str((fftIdx - 1) * 10),'-',num2str(fftIdx * 10),'% k_m_a_x'));
            axes(dualAxes(2));
            xlabel('Depth [um]');ylabel('Ratio [a.u.]');
            if fftIdx == 1
               legend([plotLeft;plotright],'Gaussian','Airy','Ratio') 
            end
        end
        % Plot high frequency average ratio
        highFreqAverage(:,fileIdx)...
            = mean(AiryFourierContent(6:10,:,fileIdx) ./ GaussianFourierContent(6:10,:,fileIdx),1);
        % Linear fit to high frequency average ratio
        linFitParams = polyfit(absRotZRange(:,fileIdx) * 1e6,highFreqAverage(:,fileIdx),1);
        linFit = polyval(linFitParams,absRotZRange(:,fileIdx) * 1e6);
        subplot(5,4,[19 20]);plot(absRotZRange(:,fileIdx) * 1e6...
            ,highFreqAverage(:,fileIdx),'r'...
            ,absRotZRange(:,fileIdx) * 1e6,linFit,'k');
        xlabel('Depth [um]');ylabel('Ratio [a.u.]');
        title('60-100% k_m_a_x average');
        drawnow;shg;
        clear highFreqAverage linFitParams linFit;
        
        
        
       
        % Stop code and request user inputs
        keyboard;
        %%% USER INPUTS REQUIRED
        %%% Uncomment the following lines and input the region of the
        %%% apodized volume over which features are present within the
        %%% image. The Fourier content graphs are only accurate if there
        %%% are features in the image volume.
        
%             % Plot limits in um
%             figDepthMin = 0;
%             figDepthMax = 35;
%             imageWidthMin = 0;
%             imageWidthMax = 47.8;

        %%% Once values have been entered, type "return" to continue code.
        
        
        
        
        % Convert plot limits to array indices
        [~,figDepthMin_pix] = max(absRotZRange(:,fileIdx) * 1e6 > figDepthMin);
        [~,figDepthMax_pix] = max(absRotZRange(:,fileIdx) * 1e6 > figDepthMax);
        figDepthMin_pix = figDepthMin_pix - 1;
        [~,imageWidthMin_pix] = max(rotYRange(:,fileIdx) * 1e6 > imageWidthMin);
        [~,imageWidthMax_pix] = max(rotYRange(:,fileIdx) * 1e6 > imageWidthMax);
        imageWidthMin_pix = imageWidthMin_pix - 1;
        
        % Display results over user specified subregion of data volume
        figure;
        % Full Gaussian projection
        subplot(5,4,[1 5]);imagesc(rotYRange(imageWidthMin_pix:imageWidthMax_pix,fileIdx) * 1e6...
            ,absRotZRange(1:figDepthMax_pix,fileIdx) * 1e6...
            ,fullGaussianProjection(1:figDepthMax_pix,imageWidthMin_pix:imageWidthMax_pix,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        xlim([imageWidthMin imageWidthMax]);ylim([figDepthMin - 10 figDepthMax]);
        title('Gaussian full image');
        % Full Airy projection
        subplot(5,4,[2 6]);imagesc(rotYRange(imageWidthMin_pix:imageWidthMax_pix,fileIdx) * 1e6...
            ,absRotZRange(1:figDepthMax_pix,fileIdx) * 1e6...
            ,fullAiryProjection(1:figDepthMax_pix,imageWidthMin_pix:imageWidthMax_pix,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        xlim([imageWidthMin imageWidthMax]);ylim([figDepthMin - 10 figDepthMax]);
        title('Airy full image');
        % Apodized Gaussian projection
        subplot(5,4,[3 7]);imagesc(rotYRange(imageWidthMin_pix:imageWidthMax_pix,fileIdx) * 1e6...
            ,absRotZRange(1:figDepthMax_pix,fileIdx) * 1e6...
            ,ROIGaussianProjection(1:figDepthMax_pix,imageWidthMin_pix:imageWidthMax_pix,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        xlim([imageWidthMin imageWidthMax]);ylim([figDepthMin - 10 figDepthMax]);
        title('Gaussian cropped image');
        % Apodized Airy projection
        subplot(5,4,[4 8]);imagesc(rotYRange(imageWidthMin_pix:imageWidthMax_pix,fileIdx) * 1e6...
            ,absRotZRange(1:figDepthMax_pix,fileIdx) * 1e6...
            ,ROIAiryProjection(1:figDepthMax_pix,imageWidthMin_pix:imageWidthMax_pix,fileIdx));axis image;
        xlabel('x+z [um]');ylabel('Tissue depth (x-z) [um]');
        xlim([imageWidthMin imageWidthMax]);ylim([figDepthMin - 10 figDepthMax]);
        title('Airy cropped image');
        colormap gray;
        % Plot Fourier content plots (left axis) and ratios (right axis)
        for fftIdx = 1:size(AiryFourierContent,1)
            subplot(5,4,8 + fftIdx);[dualAxes,plotLeft,plotright] = plotyy(absRotZRange(:,fileIdx) * 1e6...
                ,[squeeze(GaussianFourierContent(fftIdx,:,fileIdx));squeeze(AiryFourierContent(fftIdx,:,fileIdx))]...
                ,absRotZRange(:,fileIdx) * 1e6...
                ,squeeze(AiryFourierContent(fftIdx,:,fileIdx)) ./ squeeze(GaussianFourierContent(fftIdx,:,fileIdx)));
            % Increase linewidth of plots
            axes(dualAxes(1));
            Childs = get(gca,'Children');
            for n = 1:length(Childs)
                set(Childs(n),'LineWidth',3);
            end
            % Set regular plot limits and tickmarks
            ylim([min(min(squeeze(GaussianFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx)))...
                ,min(squeeze(AiryFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx)))) ...
                max(max(squeeze(GaussianFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx)))...
                ,max(squeeze(AiryFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx))))]);
            ylimits = get(gca,'YLim');
            ylimits(1) = 0;
            ylimits(2) = ceil(ylimits(2) * 10^(-1 * ceil(log10(ylimits(2))) + 1)) / 10^(-1 * ceil(log10(ylimits(2))) + 1);
            ylim(ylimits);
            yinc = (ylimits(2) - ylimits(1)) / 2;
            set(gca,'YTick',[ylimits(1):yinc:ylimits(2)]);
            xlim([figDepthMin figDepthMax]);
            xlabel('Depth [um]');ylabel('Magnitude [a.u.]');
            title(strcat(num2str((fftIdx - 1) * 10),'-',num2str(fftIdx * 10),'% k_m_a_x'));
            % Increase linewidth of plots
            axes(dualAxes(2));
            Childs = get(gca,'Children');
            for n = 1:length(Childs)
                set(Childs(n),'LineWidth',3);
            end
            % Set regular plot limits and tickmarks
            xlim([figDepthMin figDepthMax]);
            ylim([min(squeeze(AiryFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx))...
                ./ squeeze(GaussianFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx))) ...
                max(squeeze(AiryFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx))...
                ./ squeeze(GaussianFourierContent(fftIdx,figDepthMin_pix:figDepthMax_pix,fileIdx)))]);
            ylimits = get(gca,'YLim');
            ylimits(1) = floor(ylimits(1));
            ylimits(2) = ceil(ylimits(2));
            ylim(ylimits);
            yinc = (ylimits(2) - ylimits(1)) / 2;
            set(gca,'YTick',[ylimits(1):yinc:ylimits(2)]);
            xlabel('Depth [um]');ylabel('Ratio [a.u.]');
            if fftIdx == 1
               legend([plotLeft;plotright],'Gaussian','Airy','Ratio') 
            end
        end
        % Plot high frequency average ratio
        highFreqAverage...
            =mean(AiryFourierContent(6:10,figDepthMin_pix:figDepthMax_pix,fileIdx)...
            ./ GaussianFourierContent(6:10,figDepthMin_pix:figDepthMax_pix,fileIdx),1);
        highFreqAverage = highFreqAverage.';
        % Linear fit to high frequency average ratio
        [linFitParams,linFitError] = polyfit(absRotZRange(figDepthMin_pix:figDepthMax_pix,fileIdx) * 1e6,highFreqAverage,1);
        linFit = polyval(linFitParams,absRotZRange(figDepthMin_pix:figDepthMax_pix,fileIdx) * 1e6);
        subplot(5,4,[19 20]);plot(absRotZRange(figDepthMin_pix:figDepthMax_pix,fileIdx) * 1e6...
            ,highFreqAverage,'r' ...
            ,absRotZRange(figDepthMin_pix:figDepthMax_pix,fileIdx) * 1e6,linFit,'--k','LineWidth',3);
        xlim([figDepthMin figDepthMax]);
        ylim([min(highFreqAverage) max(highFreqAverage)]);
        ylimits = get(gca,'YLim');
        ylimits(1) = floor(ylimits(1));
        ylimits(2) = ceil(ylimits(2));
        ylim(ylimits);
        yinc = (ylimits(2) - ylimits(1)) / 2;
        set(gca,'YTick',[ylimits(1):yinc:ylimits(2)]);
        xlabel('Depth [um]');ylabel('Ratio [a.u.]');
        title('60-100% k_m_a_x average');
        text(figDepthMin + 10,ylimits(1) + 0.5,strcat('y=',num2str(linFitParams(1)),'x+',num2str(linFitParams(2))));
        drawnow;shg;
        %%%
        
        clear highFreqAverage;
    end   
end


