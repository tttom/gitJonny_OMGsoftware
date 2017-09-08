function makeAttnCompBioFigures_Opercula

    color_map = 'gray';
    
%         % scan 13
%         folderName = 'H:\Stored Files\NEW_LSM_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_13\2017-04-14 16_32_05.751';
% 
%         fig_handle = figure;
%         subplot_vert_total = 3;
%         subplot_horiz_total = 4;
%         
%         projDim = 3;
%         
%             %big plots
%             xRange_Start = -40;
%             xRange_End = 40;
%             yRange_Start = -100;
%             yRange_End = 150;
%             zRange_Start = -20;
%             zRange_End = 20;
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = [1 2];
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim ...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.71_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = [5 6];
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = [9 10];
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
%         
%             %small plots
%             xRange_Start = -5;
%             xRange_End = 20;
%             yRange_Start = 78;
%             yRange_End = 104;
%             zRange_Start = -20;
%             zRange_End = 20;
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = 3;
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.71_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = 7;
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%             fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');
% 
%             subplot_position = 11;
% 
%             makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
%                 ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

    % scan 13
        folderName = 'H:\Stored Files\NEW_LSM_SYSTEM_RESULTS\2017-04-14_Opercula\2Opercula\Scan_06\2017-04-14 15_51_53.940';

        fig_handle = figure;
        subplot_vert_total = 3;
        subplot_horiz_total = 4;
        
        projDim = 3;
        
            %big plots
            xRange_Start = -30;
            xRange_End = 30;
            yRange_Start = -125;
            yRange_End = 100;
            zRange_Start = -95;
            zRange_End = 20;

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');

            subplot_position = [1 2];

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim ...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat');

            subplot_position = [5 6];

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');

            subplot_position = [9 10];

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
        
            %small plots
            xRange_Start = -20;
            xRange_End = 10;
            yRange_Start = 65;
            yRange_End = 90;
            zRange_Start = -95;
            zRange_End = 20;

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');

            subplot_position = 3;

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat');

            subplot_position = 7;

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

            fileName = strcat(folderName,'\recording0_lambda532nm_alpha20_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');

            subplot_position = 11;

            makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim...
                ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

end

% actual function doing all the work
function makeAttnCompBioImage(fileName,xRange_Start,xRange_End,yRange_Start,yRange_End,zRange_Start,zRange_End,projDim,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map)

    % load data from file
    load(fileName,'restoredDataCube','xRange','yRange','zRange');
    
    % scale from [m] to [um]
    xRange = xRange * 1e6;
    yRange = yRange * 1e6;
    zRange = zRange * 1e6;
    
    xRange_start = find(xRange>=xRange_Start,1);
    xRange_end = find(xRange>=xRange_End,1);
    yRange_start = find(yRange>=yRange_Start,1);
    yRange_end = find(yRange>=yRange_End,1);
    zRange_start = find(zRange>=zRange_Start,1);
    zRange_end = find(zRange>=zRange_End,1);
    
    % make max projection
    if projDim == 3
        proj = squeeze(max(...
            restoredDataCube(xRange_start:xRange_end,yRange_start:yRange_end,zRange_start:zRange_end),[],projDim));
    else
        proj = squeeze(max(...
            restoredDataCube(xRange_start:xRange_end,yRange_start:yRange_end,zRange_start:zRange_end),[],projDim)).';
    end
    proj = proj / max(proj(:));
    
    clear restoredDataCube
    
    % plot
    figure(fig_handle);
    subplot(subplot_vert_total,subplot_horiz_total,subplot_position);
    if projDim == 1
        imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),proj);
        xlabel('x-axis [um]');
        ylabel('z-axis [um]');
    elseif projDim == 2
        imagesc(xRange(xRange_start:xRange_end),zRange(zRange_start:zRange_end),proj);
        xlabel('y-axis [um]');
        ylabel('z-axis [um]');
    else
        imagesc(yRange(yRange_start:yRange_end),xRange(xRange_start:xRange_end),proj);
        xlabel('x-axis [um]');
        ylabel('y-axis [um]');
    end
    axis image;
    colormap(color_map);
    drawnow;shg;
    
    clear xRange yRange zRange
end