% Generate PCA Test Data

imHeight = 60;
imWidth = 60;
nFrames = 1000;

[xCoords,yCoords] = meshgrid(1:imWidth,1:imHeight);
noiseFreeImageArray = zeros([imHeight,imWidth,nFrames],'double');

wobbleStrength = 3;


crossArray = zeros([size(noiseFreeImageArray,1),size(noiseFreeImageArray,2)],'double');
crossArray(floor(imHeight / 4):ceil(3 * imHeight / 4),round(imWidth / 2) - 3:round(imWidth / 2) + 3) = 1;
crossArray(round(imHeight / 2) - 3:round(imHeight / 2) + 3,floor(imWidth / 4):ceil(3 * imWidth / 4)) = 1;

simpleObjects = 1;


% rng(random_seed);


if simpleObjects
    % "Simple" objects
    for idx = 1:nFrames
        % Single ring
        noiseFreeImageArray(:,:,idx) = exp(-1 .* ((sqrt((xCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2 ...
            + (yCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2) - 15) / sqrt(2) / 3).^2);
%         % Single Gaussian
%         noiseFreeImageArray(:,:,idx) = exp(-1 .* (((xCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2 ...
%             + (yCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2) / sqrt(2) / 20).^2);
%         % Two Gaussians
%         noiseFreeImageArray(:,:,idx) = exp(-1 .* (((xCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2 ...
%             + (yCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2) / sqrt(2) / 20).^2)...
%             + exp(-1 .* (((xCoords - (imWidth / 4) - wobbleStrength * (rand(1) - 0.5)).^2 ...
%             + (yCoords - (imWidth / 4) - 2 * wobbleStrength * (rand(1) - 0.5)).^2) / sqrt(2) / 20).^4);
%         % Ring with small dot (independent)
%         noiseFreeImageArray(:,:,idx) = exp(-1 .* ((sqrt((xCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2 ...
%             + (yCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2) - 15) / sqrt(2) / 3).^2)...
%             + exp(-1 .* (sqrt((xCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2 ...
%             + (yCoords - (imWidth / 2) - wobbleStrength * (rand(1) - 0.5)).^2) / sqrt(2) / 1).^20);
    end
else
    % "Complex" objects
    for idx = 1:nFrames
        % Cross Array
        noiseFreeImageArray(:,:,idx) = interp2(xCoords,yCoords...
            ,crossArray,xCoords - wobbleStrength * (rand(1) - 0.5)...
            ,yCoords - wobbleStrength * (rand(1) - 0.5),'linear',0);
    end
end    
noiseFreeImageArray = noiseFreeImageArray / max(noiseFreeImageArray(:));
imageArray = noiseFreeImageArray + 1 .* rand(size(noiseFreeImageArray));
imageArray = imageArray / max(imageArray(:));
% for idx = 1:nFrames
%     imageArray(:,:,idx) = imageArray(:,:,idx) .* (1 + 0.1 .* randn(1)); % global fluctuations (noise and signal)
% end

[filteredArray,E,PP,val0] = function_noiseRemoval(imageArray, 20, 1:1);

val0 = val0 / sum(val0(:));
PCA_Strength = max(val0,[],1);
figure(1);plot(PCA_Strength);

figure(2);
for k = 1:20
    subplot(5,4,k);imagesc(E(:,:,k));
    xlabel('x [pixels]');
    ylabel('y [pixels]');
    axis image;
    colormap jet;
end

% for display
noiseFreeImageArray = noiseFreeImageArray - min(noiseFreeImageArray(:));
noiseFreeImageArray = noiseFreeImageArray / max(noiseFreeImageArray(:));
imageArray = imageArray - min(imageArray(:));
imageArray = imageArray / max(imageArray(:));
filteredArray = filteredArray - min(filteredArray(:));
filteredArray = filteredArray / max(filteredArray(:));

if 1
    figure(3);
    for idx = 1:nFrames
        subplot(1,3,1);
        imagesc(noiseFreeImageArray(:,:,idx));
        axis image;
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        colormap gray;
        title(strcat(num2str(idx),': Original'));
        subplot(1,3,2);
        imagesc(imageArray(:,:,idx));
        axis image;
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        colormap gray;
        title('+ noise');
        subplot(1,3,3);
        imagesc(filteredArray(:,:,idx));
        axis image;
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        colormap gray;
        title('Filtered');
        drawnow;shg;
    end
end