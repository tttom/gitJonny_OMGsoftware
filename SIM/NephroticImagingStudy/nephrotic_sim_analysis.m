function nephrotic_sim_analysis(tissue_type)

    if nargin < 1
        tissue_type = 'normal';
    end
    
    
    switch tissue_type
        case 'normal'
            data_folder = 'G:\Stored Files\SIM_PodocytePaper\DataUsedInPaper\ImageAnalysisDatasets\NormalTissue';
        case 'mcd'
            data_folder = 'G:\Stored Files\SIM_PodocytePaper\DataUsedInPaper\ImageAnalysisDatasets\MCDTissue';
        otherwise
            data_folder = 'G:\Stored Files\SIM_PodocytePaper\DataUsedInPaper\ImageAnalysisDatasets\NormalTissue';
    end
    
    % load data stack from folder
    current_directory = pwd;
    cd(data_folder);
    file_list = dir(strcat(data_folder,'\*.tif'));
    file_info = imfinfo(file_list(1).name);
    yRange = ([1:file_info.Height] - floor(file_info.Height / 2)) / file_info.YResolution * 1e-6;
    xRange = ([1:file_info.Width] - floor(file_info.Width / 2)) / file_info.XResolution * 1e-6;
    zRange = ([1:length(file_list)] - floor(length(file_list) / 2)) * 0.12 * 1e-6;  % z-plane spacing can't be determined from stack fo tiffs
    data_stack = zeros([length(yRange),length(xRange),length(zRange)],'single');
    for fileIdx = 1:length(file_list)
        data_stack(:,:,fileIdx) = single(imread(file_list(fileIdx).name));
    end
    data_stack = data_stack / (2^16 - 1);
    cd(current_directory);
    
    % define Fourier coords and FFT data
    kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) * 1 / 2 / xRange(end);
    kyRange = ([1:length(yRange)] - floor(length(yRange) / 2) - 1) * 1 / 2 / yRange(end);
    kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end);
    [kX,kY,kZ] = meshgrid(kxRange,kyRange,kzRange);
    data_FFT = fftshift(fftn(data_stack));
    
    % filter edge
    kx0 = max(kxRange) / 100;
    ky0 = max(kyRange) / 100;
    kz0 = max(kzRange) / 25;
    
    % low pass filter
    low_pass_FFT_filter = exp(-1 .* kX.^2 / 2 / kx0^2)...
        .* exp(-1 .* kY.^2 / 2 / ky0^2)...
        .* exp(-1 .* kZ.^2 / 2 / kz0^2);
    data_FFT = data_FFT .* low_pass_FFT_filter;
    
    data_stack_low_pass = real(ifftn(ifftshift(data_FFT)));
    data_stack_low_pass = data_stack_low_pass .* (data_stack_low_pass > 0); % all positive
    data_stack_low_pass = data_stack_low_pass - min(data_stack_low_pass(:)); % minimum = 0
    data_stack_low_pass = data_stack_low_pass / max(data_stack_low_pass(:)); % maximum = 1
    
    figure;
    
    subplot(4,2,[1 3]);
    imagesc(xRange * 1e6,yRange * 1e6,squeeze(max(data_stack,[],3)));axis image;
    xlabel('x [um]');ylabel('y [um]');
    title('High-res SIM projections')
    subplot(4,2,5);
    imagesc(xRange * 1e6,zRange * 1e6,squeeze(max(data_stack,[],1)).');axis image;
    xlabel('x [um]');ylabel('z [um]');
    subplot(4,2,7);
    imagesc(yRange * 1e6,zRange * 1e6,squeeze(max(data_stack,[],2)).');axis image;
    xlabel('y [um]');ylabel('z [um]');
    
    subplot(4,2,[2 4]);
    imagesc(xRange * 1e6,yRange * 1e6,squeeze(max(data_stack_low_pass,[],3)));axis image;
    xlabel('x [um]');ylabel('y [um]');
    ('Low-pass images');
    subplot(4,2,6);
    imagesc(xRange * 1e6,zRange * 1e6,squeeze(max(data_stack_low_pass,[],1)).');axis image;
    xlabel('x [um]');ylabel('z [um]');
    subplot(4,2,8);
    imagesc(yRange * 1e6,zRange * 1e6,squeeze(max(data_stack_low_pass,[],2)).');axis image;
    xlabel('y [um]');ylabel('z [um]');
    
    
    % analysis
    % 1. Take only the super-res image within the region defined by the
    % thresholded low-pass image.
    % 2. Normalise total intensity in this image to 1 and segment into 1um
    % x 1um x 1um cubes.
    % Within each cube, sum the intensity from the super-res thresholded
    % image, and multiply by the value of the low-pass image within the
    % cube.
    
    data_stack = data_stack / sum(data_stack(:));   % total intensity = 1
    data_stack_low_pass = data_stack_low_pass / sum(data_stack_low_pass(:));    % total intesntiy = 1
    
    low_pass_threshold = max(data_stack_low_pass(:)) *0.25;
    data_stack_threshold = data_stack_low_pass > low_pass_threshold;
    data_stack_thresholded = data_stack .* data_stack_threshold;
    
    figure;
    
    subplot(4,2,[1 3]);
    imagesc(xRange * 1e6,yRange * 1e6,squeeze(max(data_stack,[],3)));axis image;
    xlabel('x [um]');ylabel('y [um]');
    title('High-res SIM projections')
    subplot(4,2,5);
    imagesc(xRange * 1e6,zRange * 1e6,squeeze(max(data_stack,[],1)).');axis image;
    xlabel('x [um]');ylabel('z [um]');
    subplot(4,2,7);
    imagesc(yRange * 1e6,zRange * 1e6,squeeze(max(data_stack,[],2)).');axis image;
    xlabel('y [um]');ylabel('z [um]');
    
    subplot(4,2,[2 4]);
    imagesc(xRange * 1e6,yRange * 1e6,squeeze(max(data_stack_thresholded,[],3)));axis image;
    xlabel('x [um]');ylabel('y [um]');
    title('Thresholded images images');
    subplot(4,2,6);
    imagesc(xRange * 1e6,zRange * 1e6,squeeze(max(data_stack_thresholded,[],1)).');axis image;
    xlabel('x [um]');ylabel('z [um]');
    subplot(4,2,8);
    imagesc(yRange * 1e6,zRange * 1e6,squeeze(max(data_stack_thresholded,[],2)).');axis image;
    xlabel('y [um]');ylabel('z [um]');
    
    fraction_data_stack_within_threshold = sum(sum(sum(data_stack_threshold == 1)))...
        / size(data_stack,1) / size(data_stack,2) / size(data_stack,3);
    volume_data_stack_within_threshold = fraction_data_stack_within_threshold...
        * (xRange(end) - xRange(1)) * 1e6...
        * (yRange(end) - yRange(1)) * 1e6...
        * (zRange(end) - zRange(1)) * 1e6; % [um^3]
    
    % section into 1um^3 cubes
    x_start =  ceil(xRange(1) * 1e6);
    x_end =  floor(xRange(end) * 1e6) - 1;
    x_guage = [x_start:1:x_end] * 1e-6;
    y_start =  ceil(yRange(1) * 1e6);
    y_end =  floor(yRange(end) * 1e6) - 1;
    y_guage = [y_start:1:y_end] * 1e-6;
    z_start =  ceil(zRange(1) * 1e6);
    z_end =  floor(zRange(end) * 1e6) - 1;
    z_guage = [z_start:1:z_end] * 1e-6;
    
    results_storage = [];
    results_storage.high_freq_subvol_content = zeros([length(y_guage),length(x_guage),length(z_guage)],'single');
    results_storage.low_freq_subvol_content = zeros([length(y_guage),length(x_guage),length(z_guage)],'single');
    results_storage.threshold_subvol_content = zeros([length(y_guage),length(x_guage),length(z_guage)],'single');
    
    for x_idx = 1:length(x_guage)
        x_subvol_start = find(xRange > x_guage(x_idx),1,'first');
        x_subvol_end = find(xRange < x_guage(x_idx) + 1e-6,1,'last');
        for y_idx = 1:length(y_guage)
            y_subvol_start = find(yRange > y_guage(y_idx),1,'first');
            y_subvol_end = find(yRange < y_guage(y_idx) + 1e-6,1,'last');
            for z_idx = 1:length(z_guage)
                z_subvol_start = find(zRange > z_guage(z_idx),1,'first');
                z_subvol_end = find(zRange < z_guage(z_idx) + 1e-6,1,'last');
                results_storage.high_freq_subvol_content(y_idx,x_idx,z_idx)...
                    = sum(sum(sum(data_stack_thresholded(y_subvol_start:y_subvol_end...
                    ,x_subvol_start:x_subvol_end,z_subvol_start:z_subvol_end))));
                results_storage.low_freq_subvol_content(y_idx,x_idx,z_idx)...
                    = sum(sum(sum(data_stack_low_pass(y_subvol_start:y_subvol_end...
                    ,x_subvol_start:x_subvol_end,z_subvol_start:z_subvol_end))));
                results_storage.threshold_subvol_content(y_idx,x_idx,z_idx)...
                    = sum(sum(sum(data_stack_threshold(y_subvol_start:y_subvol_end...
                    ,x_subvol_start:x_subvol_end,z_subvol_start:z_subvol_end))));
            end
        end
    end
    
    low_pass_weighted = mean(mean(mean(results_storage.high_freq_subvol_content...
        .* results_storage.low_freq_subvol_content)));
    threshold_weighted = mean(mean(mean(results_storage.high_freq_subvol_content...
        .* results_storage.threshold_subvol_content)));
    
    [low_pass_weighted]
    [threshold_weighted]
    [volume_data_stack_within_threshold]
    
        
        
        
        
end