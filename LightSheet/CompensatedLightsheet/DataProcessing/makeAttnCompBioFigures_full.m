function makeAttnCompBioFigures_full

    %just a script to call the subfunction below lots 
    
    xRange_crop = 10;
    color_map = 'hot';
    
    % sample 2 scan 5 (77cm^-1)
        folderName = 'F:\Stored Files\DataForJonny_2017-03-30_Ferrier_Opercula\20170327_DFopercula_attnComp\sample2\Cabs_77cm-1\scan5\2017-03-27 18_15_27.067';

        fig_handle = figure;
        subplot_vert_total = 3;
        subplot_horiz_total = 1;

        fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');

        subplot_position = 1;

        yRange_start = 742;
        yRange_end = 1760;
        zRange_start = 76;
        zRange_end = 500;

        makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
            ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

        fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_ASG0.71_kSG1.41_sigmaSG8.00.mat');

        subplot_position = 2;

        makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
            ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
        
        fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');

        subplot_position = 3;

        makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
            ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
        
%         
%         folderName = 'F:\Stored Files\DataForJonny_2017-03-30_Ferrier_Opercula\20170327_DFopercula_attnComp\sample2\Cabs_50cm-1\scan3\2017-03-27 17_53_13.354';
% 
%         fig_handle = figure;
%         subplot_vert_total = 3;
%         subplot_horiz_total = 4;
% 
%         fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_ASG0.48_kSG1.41_sigmaSG8.00.mat');
% 
%         subplot_position = [1 2];
% 
%         yRange_start = 799;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 351;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 3;
% 
%         yRange_start = 884;
%         yRange_end = 1025;
%         zRange_start = 76;
%         zRange_end = 201;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 4;
% 
%         yRange_start = 1392;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 126;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.33_ASG0.77_kSG1.41_sigmaSG8.00.mat');
% 
%         subplot_position = [5 6];
% 
%         yRange_start = 799;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 351;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 7;
% 
%         yRange_start = 884;
%         yRange_end = 1025;
%         zRange_start = 76;
%         zRange_end = 201;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 8;
% 
%         yRange_start = 1392;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 126;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         fileName = strcat(folderName,'\recording0_lambda532nm_alpha7_beta100_sigmaU0.00_sigmaV0.66_ASG1.31_kSG1.41_sigmaSG8.00.mat');
% 
%         subplot_position = [9 10];
% 
%         yRange_start = 799;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 351;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 11;
% 
%         yRange_start = 884;
%         yRange_end = 1025;
%         zRange_start = 76;
%         zRange_end = 201;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);
% 
%         subplot_position = 12;
% 
%         yRange_start = 1392;
%         yRange_end = 1534;
%         zRange_start = 1;
%         zRange_end = 126;
% 
%         makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
%             ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map);

end




% actual function doing all the work
function makeAttnCompBioImage(fileName,xRange_crop,yRange_start,yRange_end,zRange_start,zRange_end...
    ,fig_handle,subplot_vert_total,subplot_horiz_total,subplot_position,color_map)

    % load data from file
    load(fileName,'restoredDataCube','yRange','zRange');
    
    % scale from [m] to [um]
    yRange = yRange * 1e6;
    zRange = zRange * 1e6;
    
    % make xz (yz in matlab coords) max projection
    xz_proj = squeeze(max(...
        restoredDataCube(xRange_crop + 1:size(restoredDataCube,1) - xRange_crop,yRange_start:yRange_end,zRange_start:zRange_end),[],1)).';
    % ensure positivity
    xz_proj = xz_proj / max(xz_proj(:));
    
    clear restoredDataCube
    
    % plot
    figure(fig_handle);
    subplot(subplot_vert_total,subplot_horiz_total,subplot_position);
    imagesc(yRange(yRange_start:yRange_end),zRange(zRange_start:zRange_end),xz_proj);
    axis image;
    xlabel('x-axis [um]');
    ylabel('z-axis [um]');
    colormap(color_map);
    drawnow;shg;
    
    clear yRange zRange
end