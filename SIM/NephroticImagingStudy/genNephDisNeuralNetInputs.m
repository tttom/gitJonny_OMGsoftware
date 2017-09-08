%%% Takes large input tif files and segments these into smaller datasets
%%% for neural network analysis.

function genNephDisNeuralNetInputs(parentFolder)

    if nargin < 1
%         parentFolder = 'C:\Users\Jonathan Nylk\Desktop\Podocyte Images\Healthy';
        parentFolder = 'C:\Users\Jonathan Nylk\Desktop\Podocyte Images\MCD';
    end
    
    outputFileCounter = 1;
    
    currentFolder = pwd;
    
    cd(parentFolder)
    fileList = dir('.\MasterImages\*.tif');

    for fileIdx = 1:length(fileList)
        inputImage = double(imread(strcat('.\MasterImages\',fileList(fileIdx).name)));
        inputImage = inputImage / max(inputImage(:));
        % 1px = 0.0318um
        
        %cut into 5um squares
        for rIdx = 1:6
            for cIdx = 1:6
                subImage = inputImage(41 + (cIdx-1) * 157: 40 + cIdx * 157, 41 + (rIdx-1) * 157: 40 + rIdx * 157);
%                 totalImageIntensity(outputFileCounter) = sum(subImage(:));
%                 maxImageIntensity(outputFileCounter) = max(subImage(:));
%                 noPixelsAboveThreshold(outputFileCounter) = sum(sum(subImage >= 0.1));
                imageNumber = num2str(1000 + outputFileCounter);
                outputFileCounter = outputFileCounter + 1;
                imwrite(subImage,strcat('image_',imageNumber(2:end),'.png'),'bitdepth',16);
                
            end
        end
    end
    
%     figure;
%     subplot(2,1,1);plot(totalImageIntensity);
%     subplot(2,1,2);plot(totalImageIntensity / max(totalImageIntensity));
%     
%     figure;
%     plot(maxImageIntensity);
%     
%     figure;
%     plot(noPixelsAboveThreshold);
    
    cd(currentFolder);
end