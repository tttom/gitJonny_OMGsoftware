function [data_structure,outmat] = spotfinder3D(inmatrix,threshold,yRange,xRange,zRange,ySearchLimits)
%inmatrix is a recorded light sheet volume

%threshold is the minimum pixel intensity the algorithm should look for.
%threshold should be high enough that it is only found in bead PSFs, not
%the background.

if nargin < 2
    threshold = 0.5;
end

if nargin < 3
    mag = 40/36.74;
    xscale = 0.1709*mag; %microns/pixel
    xRange = ([1:size(inmatrix,1)] - floor(size(inmatrix) / 1)) * xscale; %microns (centred at zero)
end

if nargin < 4
    mag = 40/36.74;
    xscale = 0.1709*mag; %microns/pixel
    yRange = ([1:size(inmatrix,2)] - floor(size(inmatrix) / 2)) * xscale; %microns (centred at zero)
end

if nargin < 5
    zscale = 0.1996; %microns/pixel [magnification doesn't change scan depth]
    zRange = [1:size(inmatrix,3)] * zscale; % microns
end

if nargin < 6
    ySearchLimits = [0 0];
    ySearchLimits(1) = yRange(floor(length(yRange) / 4)); % microns
    ySearchLimits(2) = yRange(floor(3 * length(yRange) / 4)); % microns
end


%outmat is an output array containing the (r,c,z) indices of every spot
%center found, the z depth in microns of that spot, its FWHMs in the
%r, c, and z directions, and the RMSE values of the Gaussian fits.

xscale = yRange(2) - yRange(1);
zscale = zRange(2) - zRange(1);

yStartPixel = find(yRange > ySearchLimits(1),1);
yEndPixel = find(yRange > ySearchLimits(2),1);

outmat = {};
coords = [];
radii = [];

%Make a yz projection
% inmatrix = inmatrix(21:492,641:1408,186:315);
inmatrix = inmatrix(1+20:size(inmatrix,1)-20,yStartPixel:yEndPixel,1+20:size(inmatrix,3)-20);
inmatrix = inmatrix .* (inmatrix > 0); % ensure positive values only
maxIntensity = max(inmatrix(:));
inmatrix = inmatrix / maxIntensity; % scale between [0:1];
h = size(inmatrix,3);
xyproj = squeeze(max(inmatrix,[],3)); %y is in laboratory coordinates
xzproj = squeeze(max(inmatrix,[],1)).'; %y is in laboratory coordinates
[r,c] = size(xyproj); %rows (r) are y, columns (c) are x.

projectionFig = figure();
xyproj_Fig = subplot(2,1,1);imshow(xyproj);
caxis([min(min(xyproj)),max(max(xyproj))]);
xlabel('x /pixels');
ylabel('y /pixels');
xzproj_Fig = subplot(2,1,2);imshow(xzproj);
caxis([min(min(xzproj)),max(max(xzproj))]);
xlabel('x /pixels');
ylabel('z /pixels');
drawnow;shg;

projectionFig = figure();
xyproj_Fig = subplot(2,1,1);imshow(xyproj > threshold);
caxis([0,1]);
xlabel('x /pixels');
ylabel('y /pixels');
xzproj_Fig = subplot(2,1,2);imshow(xzproj > threshold);
caxis([0,1]);
xlabel('x /pixels');
ylabel('z /pixels');
drawnow;shg;

no_elements = numel(inmatrix);
no_elements_above_threshold = sum(sum(sum(inmatrix >= threshold)));
fprintf('Total pixels in volume: %d.\n',no_elements);
fprintf('Number of pixels above threshold: %d.\n',no_elements_above_threshold);
fprintf('Number of pixels above threshold: %2.2f %%.\n',no_elements_above_threshold / no_elements * 100);

projectionFig_good_bad = figure();
xyproj_Fig_good_bad = subplot(2,1,1);imshow(xyproj);
caxis([min(min(xyproj)),max(max(xyproj))]);
xlabel('x /pixels');
ylabel('y /pixels');
xzproj_Fig_good_bad = subplot(2,1,2);imshow(xzproj);
caxis([min(min(xzproj)),max(max(xzproj))]);
xlabel('x /pixels');
ylabel('z /pixels');
drawnow;shg;

projectionFig_good_only = figure();
xyproj_Fig_good_only = subplot(2,1,1);imshow(xyproj);
caxis([min(min(xyproj)),max(max(xyproj))]);
xlabel('x /pixels');
ylabel('y /pixels');
xzproj_Fig_good_only = subplot(2,1,2);imshow(xzproj);
caxis([min(min(xzproj)),max(max(xzproj))]);
xlabel('x /pixels');
ylabel('z /pixels');
drawnow;shg;

% threshold = input('Enter numerical threshold value: ');
% threshold = 0.00025;
% threshold = 0.0001;

%Raster scan along every nth row. Ignore the first and last 10 rows and columns.
timeForSearch = tic;
for i = 11:3:r-10 %Every 8th row means the (old) 9x9x9 boxes just overlap.
	for j = 11:c-10
		if xyproj(i,j) >= threshold
			%Go to the same (r,z) coordinates in inmatrix and scan along the third coordinate, c.
			for m = 11:h-10
				if inmatrix(i,j,m) > threshold
                    PSFspot = inmatrix(i-10:i+10,j-10:j+10,m-10:m+10); %a 9x9x9 box is ~1.4 microns/side, just larger than diffraction-limited
                    %Search for pixels of higher value. If they are present, re-center on that pixel, then draw 						
                    %a 21x21x21 box around the central pixel.
                    [~,I] = max(PSFspot(:));
                    [d1, d2, d3] = ind2sub([21,21,21],I);
                    newi = i-(11-d1);
                    %logi = inrange(newi,[outmat{end,1}(1,1)-1,outmat{end,1}(1,1)+1]);
                    newj = j-(11-d2);
                    %logj = inrange(newj,[outmat{end,1}(1,2)-1,outmat{end,1}(1,2)+1]);
                    newm = m-(11-d3);
                    %logm = inrange(newm,[outmat{end,1}(1,3)-1,outmat{end,1}(1,3)+1]);
                    if isempty(outmat) || (inrange(newi,[outmat{end,1}(1,1)-1,outmat{end,1}(1,1)+1]) && inrange(newj,[outmat{end,1}(1,2)-1,outmat{end,1}(1,2)+1]) && inrange(newm,[outmat{end,1}(1,3)-1,outmat{end,1}(1,3)+1]))
                        dist = [newi,newj,newm,r-newi,c-newj,h-newm];
                        %Check whether ALL values in dist are >= 11
                        log = dist >= 11;
                        rng = all(log);
                        %Skip spots closer to any edge than 10 pixels.
                        if rng == 1
                            sep = 10;
                            PSFspot = inmatrix(newi-10:newi+10,newj-10:newj+10,newm-10:newm+10);
                        else
                            sep = 0;
                        end
                        %Fit with Gaussian distributions along r, c, and z.
                        if sep ~= 0
            %                             axis_ = -sep:sep;
            %                             axis_ = axis_';
            %                             [FWHMr,RMSEr,yOffsetr,Ampr] = FWHMcalc(axis_,PSFspot(sep+1,:,sep+1)');
            %                             [FWHMc,RMSEc,yOffsetc,Ampc] = FWHMcalc(axis_,PSFspot(:,sep+1,sep+1));
            %                             [FWHMh,RMSEh,yOffseth,Amph] = FWHMcalc(axis_,squeeze(PSFspot(sep+1,sep+1,:)));
            %                             yOffset = mean([yOffsetr yOffsetc yOffseth]);
            %                             yOffset_variance = var([yOffsetr yOffsetc yOffseth]);
            %                             peakI = max(PSFspot(:)) - yOffset; %Height of peak in arbitrary units;
            %                             Amp = mean([Ampr,Ampc,Amph]); %Should match peakI;
            %                             Amp_variance = var([Ampr,Ampc,Amph]);
            %                             localSNR = peakI/yOffset;
            %                     		%Calculate FWHMs and convert to microns.
            %                         	FWHMr = FWHMr*xscale;
            %                             FWHMc = FWHMc*xscale;
            %         					FWHMh = FWHMh*zscale;
            %             				%Convert the z coordinate to microns.
            %                 			depth = m*zscale;
            %                             %For Airy pixel reassignment data, convert the
            %                             %x coordinate to microns, too.
            %                             disp = j*xscale;
            %                     		%Store (r,c,z), depth, and FWHMs to outmat.
            %                             data{1,1} = [newi newj newm];
            %                             data{1,2} = depth;
            %                             data{1,3} = FWHMr;
            %                             data{1,4} = FWHMc;
            %                             data{1,5} = FWHMh;
            %                             data{1,6} = [RMSEr RMSEc RMSEh];
            % %                             data{1,7} = yOffset;
            % %                             data{1,8} = peakI;
            %                             data{1,7} = localSNR;
            %                             data{1,8} = disp; % relative x-axis position
            %                             data{1,9} = yOffset;
            %                             data{1,10} = yOffset_variance;
            %                             data{1,11} = Amp;
            %                             data{1,12} = Amp_variance;
            %                             data{1,13} = disp + yRange(yStartPixel - 1); % absolute x-axis position in full datacube
            %                             data{1,14} = peakI;
            %                             outmat = cat(1,outmat,data);
            %                             %Store x,y coordinates and FWHMx in vectors to plot.
            % %                             coords = [coords; newj newi];
            % %                             radii = [radii; FWHMc/(2*xscale)];
                            [FWHMr,FWHMc,FWHMh,yOffset,Amp] = FWHMcalc3D(PSFspot);
                            % Calculate FWHMs and convert to microns.
                            FWHMr = FWHMr*xscale;
                            FWHMc = FWHMc*xscale;
                            FWHMh = FWHMh*zscale;
                            %Convert the z coordinate to microns.
                            depth = m*zscale;
                            % Calculate SNR
                            localSNR = Amp / yOffset;
                            disp = j*xscale;
                            peakI = max(PSFspot(:));
                            %Store (r,c,z), depth, and FWHMs to outmat.
                            data{1,1} = [newi newj newm];
                            data{1,2} = depth;
                            data{1,3} = FWHMr;
                            data{1,4} = FWHMc;
                            data{1,5} = FWHMh;
%                             data{1,6} = [RMSEr RMSEc RMSEh];
                            data{1,6} = [];
%                             data{1,7} = yOffset;
%                             data{1,8} = peakI;
                            data{1,7} = localSNR;
                            data{1,8} = disp; % relative x-axis position
                            data{1,9} = yOffset;
%                             data{1,10} = yOffset_variance;
                            data{1,10} = [];
                            data{1,11} = Amp;
%                             data{1,12} = Amp_variance;
                            data{1,12} = [];
                            data{1,13} = disp + yRange(yStartPixel - 1); % absolute x-axis position in full datacube
                            data{1,14} = peakI;
                            outmat = cat(1,outmat,data);
                        end
                    end
                end
            end
        end
    end
end

elapsedTime = toc(timeForSearch);
fprintf('Elapsed time is %d seconds.\n',round(elapsedTime));

figure(projectionFig_good_bad);
coords = cell2mat(outmat(:,1));
xycoords = circshift(coords(:,1:2),[0 1]);
xzcoords = coords(:,2:3);
radii = cell2mat(outmat(:,4)) / 2 / xscale;
axes(xyproj_Fig_good_bad);
viscircles(xycoords,radii,'EdgeColor','r');
axes(xzproj_Fig_good_bad);
viscircles(xzcoords,radii,'EdgeColor','r');
drawnow;shg;

%remove duplicate entries from outmat
coordinates = outmat(:,1);
for i=1:length(coordinates)
    coordinates{i}=num2str(coordinates{i});
end
[~,IA,~]=unique(coordinates,'stable'); %Outputs are the matrix of unique values (suppressed), the indices in the original matrix, and the indices in the new matrix (suppressed).
outmat=outmat(IA,:);

figure;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,3)),'b');hold on;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,4)),'g');
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,5)),'r');
xlabel('x [um]');
ylabel('FWHM [um]');
title('FWHM in row, column, and depth (before filtering)');

figure;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,7)),'b');
xlabel('x [um]');
ylabel('SNR/SBR [a.u.]');
title('SNR/SBR (before filtering)');

% filter spurious results
    fprintf('Initial number of spots found is %d.\n',length(outmat));
    % lateral resolution outwith [1 1.5]um
    keep_indices = find(cell2mat(outmat(:,3)) >= 1 & cell2mat(outmat(:,3)) <= 1.5); % FWHMr
    outmat = outmat(keep_indices,:);
    fprintf('Filtering spurious FWHM in "r" direction.\nNumber of spots remaining is %d.\n',length(outmat));
    keep_indices = find(cell2mat(outmat(:,4)) >= 1 & cell2mat(outmat(:,4)) <= 1.5); % FWHMc
    outmat = outmat(keep_indices,:);
    fprintf('Filtering spurious FWHM in "c" direction.\nNumber of spots remaining is %d.\n',length(outmat));
%     % lateral resolution outwith [1 2]um
%     keep_indices = find(cell2mat(outmat(:,3)) >= 1 & cell2mat(outmat(:,3)) <= 2); % FWHMr
%     outmat = outmat(keep_indices,:);
%     fprintf('Filtering spurious FWHM in "r" direction.\nNumber of spots remaining is %d.\n',length(outmat));
%     keep_indices = find(cell2mat(outmat(:,4)) >= 1 & cell2mat(outmat(:,4)) <= 2); % FWHMc
%     outmat = outmat(keep_indices,:);
%     fprintf('Filtering spurious FWHM in "c" direction.\nNumber of spots remaining is %d.\n',length(outmat));
    
    %axial resolution outwith [0.5 3]um
    keep_indices = find(cell2mat(outmat(:,5)) >= 0.5 & cell2mat(outmat(:,5)) <= 2.5); % FWHMh
    outmat = outmat(keep_indices,:);
    fprintf('Filtering spurious FWHM in "h" direction.\nNumber of spots remaining is %d.\n',length(outmat));
    
figure;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,3)),'b');hold on;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,4)),'g');
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,5)),'r');
xlabel('x [um]');
ylabel('FWHM [um]');
title('FWHM in row, column, and depth (after filtering)');

figure;
scatter(cell2mat(outmat(:,13)),cell2mat(outmat(:,7)),'b');
xlabel('x [um]');
ylabel('SNR/SBR [a.u.]');
title('SNR/SBR (after filtering)');

    % bin FWHM in 10um x-axis intervals
    num_bins = floor((yRange(yEndPixel) - yRange(yStartPixel)) / 10);
    bin_step_size = floor(length(yRange(yStartPixel:yEndPixel)) / num_bins) - 1;
    y_mean = zeros([1,num_bins]);
    FWHMr_mean = zeros([1,num_bins]);
    FWHMc_mean = zeros([1,num_bins]);
    FWHMh_mean = zeros([1,num_bins]);
    FWHMr_std = zeros([1,num_bins]);
    FWHMc_std = zeros([1,num_bins]);
    FWHMh_std = zeros([1,num_bins]);
    SNR1_mean = zeros([1,num_bins]);
    SNR1_std = zeros([1,num_bins]);
    SNR2_mean = zeros([1,num_bins]);
    SNR2_std = zeros([1,num_bins]);
    peakI_mean = zeros([1,num_bins]);
    peakI_std = zeros([1,num_bins]);
    for bin_idx = 1:num_bins
        yRange_start = max(yRange(yStartPixel + bin_step_size * (bin_idx - 1) + 1),yRange(yStartPixel));
        yRange_end = min(yRange(yStartPixel + bin_step_size * (bin_idx)),yRange(yStartPixel) + length(yRange(yStartPixel:yEndPixel)));
        y_mean(bin_idx) = mean([yRange_start yRange_end]);
        bin_spot_locations = find(cell2mat(outmat(:,13)) >= yRange_start & cell2mat(outmat(:,13)) <= yRange_end);
        FWHMr_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,3)));
        FWHMc_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,4)));
        FWHMh_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,5)));
        FWHMr_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,3)));
        FWHMc_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,4)));
        FWHMh_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,5)));
        SNR1_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,7)));
        SNR1_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,7)));
        SNR2_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,11)) ./ cell2mat(outmat(bin_spot_locations,9)));
        SNR2_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,11)) ./ cell2mat(outmat(bin_spot_locations,9)));
        peakI_mean(bin_idx) = mean(cell2mat(outmat(bin_spot_locations,14)));
        peakI_std(bin_idx) = std(cell2mat(outmat(bin_spot_locations,14)));
    end

figure;
errorbar(y_mean,FWHMr_mean,FWHMr_std,'o');
xlabel('x [um]');
ylabel('FWHM [um]');
title('FWHM - row (x?)');
figure;
errorbar(y_mean,FWHMc_mean,FWHMc_std,'o');
xlabel('x [um]');
ylabel('FWHM [um]');
title('FWHM - row (y?)');
figure;
errorbar(y_mean,FWHMh_mean,FWHMh_std,'o');
xlabel('x [um]');
ylabel('FWHM [um]');
title('FWHM - row (z)');
figure;
errorbar(y_mean,SNR1_mean,SNR1_std,'o');
xlabel('x [um]');
ylabel('SNR/SBR 1 [a.u.]');
title('SNR/SBR 1');
errorbar(y_mean,SNR2_mean,SNR2_std,'o');
xlabel('x [um]');
ylabel('SNR/SBR 2 [a.u.]');
title('SNR/SBR 2');

figure;
errorbar(y_mean,peakI_mean,peakI_std,'o');
xlabel('x [um]');
ylabel('peak intensity [um]');
title('peak intensity (x?)');
    
figure(projectionFig_good_bad);
coords = cell2mat(outmat(:,1));
xycoords = circshift(coords(:,1:2),[0 1]);
xzcoords = coords(:,2:3);
radii = cell2mat(outmat(:,4)) / 2 / xscale;
axes(xyproj_Fig_good_bad);
viscircles(xycoords,radii,'EdgeColor','g');
axes(xzproj_Fig_good_bad);
viscircles(xzcoords,radii,'EdgeColor','g');
drawnow;shg;

figure(projectionFig_good_only);
coords = cell2mat(outmat(:,1));
xycoords = circshift(coords(:,1:2),[0 1]);
xzcoords = coords(:,2:3);
radii = cell2mat(outmat(:,4)) / 2 / xscale;
axes(xyproj_Fig_good_only);
viscircles(xycoords,radii,'EdgeColor','g');
axes(xzproj_Fig_good_only);
viscircles(xzcoords,radii,'EdgeColor','g');
drawnow;shg;

data_structure = [];
data_structure.y_mean = y_mean;
data_structure.FWHMr_mean = FWHMr_mean;
data_structure.FWHMr_std = FWHMr_std;
data_structure.FWHMc_mean = FWHMc_mean;
data_structure.FWHMc_std = FWHMc_std;
data_structure.FWHMh_mean = FWHMh_mean;
data_structure.FWHMh_std = FWHMh_std;
data_structure.SNR1_mean = SNR1_mean;
data_structure.SNR1_std = SNR1_std;
data_structure.SNR2_mean = SNR2_mean;
data_structure.SNR2_std = SNR2_std;
data_structure.peakI_mean = peakI_mean;
data_structure.peakI_std = peakI_std;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logout = inrange(number,range)
    %Checks whether number is within the range defined by range = [u,v],
    %inclusive, and outputs a logical 1 or 0. NOTE: currently it returns 1
    %if number is OUTSIDE of range.
    if number < range(1) || number > range(2)
        logout = 1;
    else
        logout = 0;
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FWHM,RMSE,yOffset,Amp] = FWHMcalc(ax,vals)
    %First, subtract that background value from vals to try to eliminate
    %the huge errors I've been getting.
    vals = double(vals-min(vals));
    %Now try scaling to the range 0-1 since the values are very small.
    %It'll be a constant multiplier, so the FWHM shouldn't change--the
    %original half max will still be half the (new) max.
    mv = max(vals);
    vals = vals/mv;
    [parameters,~,RMSE] = GsnFit(ax,vals); %parameters = [amplitude, sigma (not ^2!), center, yOffset]
    FWHM = 2*sqrt(2*log(2))*parameters(2);
    yOffset = parameters(4)*mv; %Scaled back to the original intensity range.
    Amp = parameters(1)*mv;
    %The ~ outputs are the adjusted R^2 value of the fit and its RMSE.
    %Write them to variables if you need more fit quality information.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FWHMr,FWHMc,FWHMh,yOffset,Amp] = FWHMcalc3D(input_matrix)
    % Relies on getting a 3d array as input matrix
    
    % Define zero-centred coordiante system (in pixels)
    dim1Range = [1:size(input_matrix,1)] - floor(size(input_matrix,1) / 2) - 1;
    dim2Range = [1:size(input_matrix,2)] - floor(size(input_matrix,2) / 2) - 1;
    dim3Range = [1:size(input_matrix,3)] - floor(size(input_matrix,3) / 2) - 1;
    [dim1Coords,dim2Coords,dim3Coords] = ndgrid(dim1Range,dim2Range,dim3Range);
    
    % Define 3D Gaussian function
    % params = [(1) amplitude, (2) dim1_centre, (3) dim2_centre, (4) dim3_centre, (5) dim1_width, (6) dim2_width, (7) dim3_width, (8) y_offset];
    functor_3DGaussian = @(dim1Coords,dim2Coords,dim3Coords,params) params(1) ...
        .* exp(-1 .* ((dim1Coords - params(2)).^2) ./ 2 ./ (params(5)).^2) ...
        .* exp(-1 .* ((dim2Coords - params(3)).^2) ./ 2 ./ (params(6)).^2) ...
        .* exp(-1 .* ((dim3Coords - params(4)).^2) ./ 2 ./ (params(7)).^2) ...
        + params(8);
    % Function to minimise to fit 3D Gaussian to data
    functor_3DGaussian_minimisation = @(dim1Coords,dim2Coords,dim3Coords,params,data)...
        sum(sum(sum((functor_3DGaussian(dim1Coords,dim2Coords,dim3Coords,params) - data).^2)));
    % input guess params
    params0 = zeros([1,8]);
    params0(1) = max(input_matrix(:)); % guess Gaussian amplitude
    params0(2:4) = 0; % guess at dim1/2/3 centre offsets (assumed to be centred)
    params0(5:7) = 3; % guess of dim1/2/3 widths (based on average FWHM values of 1.25um converted back into widths in pixels
    params0(8) = min(input_matrix(:)); % guess y_offset
    % fitting
    [params,~,exitflag] = fminsearch(@(params) functor_3DGaussian_minimisation(dim1Coords,dim2Coords,dim3Coords,params,input_matrix),params0);
    if exitflag ==1
        % if converged on a solution, use these values or closest
        % p
        FWHMr = 2 * sqrt(2 * log(2)) * max(params(5),-1 * params(5));
        FWHMc = 2 * sqrt(2 * log(2)) * max(params(6),-1 * params(6));
        FWHMh = 2 * sqrt(2 * log(2)) * max(params(7),-1 * params(7));
        yOffset = params(8);
%         yOffset = max(0,params(8));
        Amp = params(1);
%         Amp = min(params(1),max(input_matrix(:)));
    else
        % if can't fit, make sure results are rejected by setting to zero
        % (must be non-zero for circle plotting - make really small)
        FWHMr = 0.001;
        FWHMc = 0.001;
        FWHMh = 0.001;
        yOffset = 0;
        Amp = 0;
    end
end


