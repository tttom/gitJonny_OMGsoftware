function outmat = spotfinder3D(inmatrix,threshold)
%inmatrix is a recorded light sheet volume

%threshold (user input in line 26) is the minimum pixel intensity the algorithm should look for.
%threshold should be high enough that it is only found in bead PSFs, not
%the background.

%outmat is an output array containing the (r,c,z) indices of every spot
%center found, the z depth in microns of that spot, its FWHMs in the
%r, c, and z directions, and the RMSE values of the Gaussian fits.

mag = 40/36.74;
xscale = 0.1709*mag; %microns/pixel
zscale = 0.1996*mag; %microns/pixel

outmat = {};
coords = [];
radii = [];

%Make a yz projection
inmatrix = inmatrix(21:492,641:1408,186:315);
h = size(inmatrix,3);
xyproj = squeeze(max(inmatrix,[],3)); %y is in laboratory coordinates
[r,c] = size(xyproj); %rows (r) are y, columns (c) are x.
imshow(xyproj);
caxis([min(min(xyproj)),max(max(xyproj))]);
xlabel('Distance /microns');
ylabel('Depth /microns');
% threshold = input('Enter numerical threshold value: ');
% threshold = 0.00025;
% threshold = 0.0001;

%Raster scan along every nth row. Ignore the first and last 10 rows and columns.
for i = 11:2:r-10 %Every 8th row means the (old) 9x9x9 boxes just overlap.
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
                            axis = -sep:sep;
                            axis = axis';
                            [FWHMr,RMSEr,yOffsetr] = FWHMcalc(axis,PSFspot(sep+1,:,sep+1)');
                            [FWHMc,RMSEc,yOffsetc] = FWHMcalc(axis,PSFspot(:,sep+1,sep+1));
                            [FWHMh,RMSEh,yOffseth] = FWHMcalc(axis,squeeze(PSFspot(sep+1,sep+1,:)));
                            yOffset = mean([yOffsetr yOffsetc yOffseth]); %%%add var to check
                            peakI = max(PSFspot) - yOffset; %Height of peak in arbitrary units; comp to A
                            localSNR = peakI/yOffset;
                    		%Calculate FWHMs and convert to microns.
                        	FWHMr = FWHMr*xscale;
                            FWHMc = FWHMc*xscale;
        					FWHMh = FWHMh*zscale;
            				%Convert the z coordinate to microns.
                			depth = m*zscale;
                            %For Airy pixel reassignment data, convert the
                            %x coordinate to microns, too.
                            disp = j*xscale;
                    		%Store (r,c,z), depth, and FWHMs to outmat.
                            data{1,1} = [newi newj newm];
                            data{1,2} = depth;
                            data{1,3} = FWHMr;
                            data{1,4} = FWHMc;
                            data{1,5} = FWHMh;
                            data{1,6} = [RMSEr RMSEc RMSEh];
%                             data{1,7} = yOffset;
%                             data{1,8} = peakI;
                            data{1,7} = localSNR;
                            data{1,8} = disp;
                            outmat = cat(1,outmat,data);
                            %Store x,y coordinates and FWHMx in vectors to plot.
                            coords = [coords; newj newi];
                            radii = [radii; FWHMc/(2*xscale)];
                        end
                    end
                end
            end
        end
    end
end
viscircles(coords,radii,'Color','r');

%remove duplicate entries from outmat
coordinates = outmat(:,1);
for i=1:length(coordinates)
    coordinates{i}=num2str(coordinates{i});
end
[~,IA,~]=unique(coordinates,'stable'); %Outputs are the matrix of unique values (suppressed), the indices in the original matrix, and the indices in the new matrix (suppressed).
outmat=outmat(IA,:);

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
function [FWHM,RMSE,yOffset] = FWHMcalc(ax,vals)
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
    %The ~ outputs are the adjusted R^2 value of the fit and its RMSE.
    %Write them to variables if you need more fit quality information.
end