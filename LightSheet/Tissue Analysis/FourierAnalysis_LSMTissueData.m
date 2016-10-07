%%% Function:           FourierAnalysis_LSMTissueData.m
%%% Author:             Jonathan Nylk (University of St Andrews)
%%% Created:            07/03/2016
%%% Description:        This function performs analysis of Fourier content
%%%                     versus depth into a 3D image acquired using
%%%                     light-sheet microscopy as described in "Enhancement
%%%                     of image quality and imaging depth with Airy
%%%                     light-sheet microscopy in cleared and non-cleared
%%%                     neural tissue" by J. Nylk et al (2016) and is used
%%%                     to perform the analysis in Figures 3 and 5 from
%%%                     this publication.
%%%                     
%%%                     The data volumes corresponding to a single lateral
%%%                     tissue location, but increasing depth within the
%%%                     tissue are sequentially loaded, and rotated to
%%%                     align the volume with the tissue coordinates.
%%%                     Each x-y plane of the volume is then Fourier
%%%                     transformed and the signal within specified
%%%                     spectral widnows is summed.

%%%                     The resulting data is stored in file and can be
%%%                     dispalyed using the function
%%%                     "displayFourierAnalysis_LSMTissueData_Results.m".
%%%
%%%
%%% Inputs:             FilePathCell: Cell containing filepaths to
%%%                     dataset folders, surface dataset
%%%                     first, then datasets going deeper into tissue.
%%%                     OutputFolder: String containing folder path for
%%%                     saving data to.
%%%                     
%%%
%%% Updates (latest first):
%%%         28/07/2016: Properly documented code.
%%%
%%% END

function FourierAnalysis_LSMTissueData(FilePathCell,OutputFolder,Flag_flippedImage)


    % Default inputs
    if nargin < 3
        Flag_flippedImage = 0;  % Flag for datasets acquired with adapted imaging pathway
                                % (Image is flipped left-to-right)
                                % (Not required for data in Figures 2-5).
    end
    
    if nargin < 2
        % Generate unqiue output folder name based on current time
        exactTime = clock;
        OutputFolder = strcat(pwd,'/',num2str(exactTime(1)),'_'...
            ,num2str(exactTime(2)),'_',num2str(exactTime(3)),'_'...
            ,num2str(exactTime(4)),'_',num2str(exactTime(5)),'_'...
            ,num2str(exactTime(6)));
    end
    
    if nargin < 1
        % Default input folders, surface dataset first
        FilePathCell = {'First full file path'...
            ,'Second full file path'...
            ,'Third full filepath'};
    end



    for fileIdx = 1:length(FilePathCell)
        results = struct();
        % Load data
        for typeIdx = 1:2 %(Gaussian/Airy)
            folderName = FilePathCell{fileIdx};
            if typeIdx == 1
                % Load Gaussian datasets sequentially
                fileName = strcat(folderName,'/recording0_lambda532nm_alpha0_beta100.mat');
                disp(strcat('Loading file:',fileName));
                load(fileName,'recordedImageStack','xRange','yRange','zRange');
                disp('File loaded');
                restoredDataCube = recordedImageStack;
                clear recordedImageStack;
            elseif typeIdx == 2
                 % Load Airy datasets sequentially
                 if Flag_flippedImage == 1
                    fileName = strcat(folderName,'/recording0_lambda532nm_alpha-7_beta100.mat');
                 else
                    fileName = strcat(folderName,'/recording0_lambda532nm_alpha7_beta100.mat');
                 end
                    disp(strcat('Loading file:',fileName));
                    load(fileName,'restoredDataCube','xRange','yRange','zRange');
                    disp('File loaded');
            end
            
            % Notes on MATLAB coordinates (xRange,yRange,zRange) and
            % real-space coordinates.
                % xRange is y-dimension [metres]
                % yRange is x-dimension [metres]
                % zRange is actually z-dimension [metres]
            
            % Flip camera images left-to-right (along x-axis (yRange))
            if Flag_flippedImage == 1
               restoredDataCube = flipdim(restoredDataCube,2);
            end
            
            % Set edge values to zero (remove edge artefacts)
                restoredDataCube(1:10,:,:) = 0;
                restoredDataCube(end-9:end,:,:) = 0;
                restoredDataCube(:,1:10,:) = 0;
                restoredDataCube(:,end-9:end,:) = 0;
                restoredDataCube(:,:,1:10) = 0;
                restoredDataCube(:,:,end-9:end) = 0;
                
            % x-z projecton of data for feature matching and depth
            % determination
                fullProjection = squeeze(max(restoredDataCube,[],1));
                [fullProjection,rotYRange,rotZRange] = rotate2DArray(fullProjection,-45/360*2*pi,yRange,zRange,0,0);
                
            % Apply Gaussian filter to (smoothly) crop data around Gaussian focus
                [Y,~] = meshgrid(yRange,xRange);
                dataWidth = 8e-6;   % Gaussian half-Field-of-View [metres]
                dataFilter = exp(-(Y / 2 / dataWidth).^2);
                for n=1:size(restoredDataCube,3)
                    restoredDataCube(:,:,n) = restoredDataCube(:,:,n) .* dataFilter;
                end
                
            % Rotating dataset by 45 deg
                disp('Rotating data')
                for xIdx = 1:length(xRange)
                    [restoredDataCube(xIdx,:,:),~,~] = rotate2DArray(squeeze(restoredDataCube(xIdx,:,:)),-45 / 360 * 2 * pi,yRange,zRange,0,0);
                    % Notify every 100 iterations
                        if mod(xIdx,100) == 0
                            disp(strcat('Rotating plane (',num2str(xIdx),') of(',num2str(length(xRange)),')'));
                        end
                end
                disp('Data rotated');
                
            % x-z projection of ROI'ed (apodized) data
                ROIProjection = squeeze(max(restoredDataCube,[],1));
                
            % Define normalised Fourier space coordinates (normalised to lambda/2NA)
                kxRange = ([1:length(xRange)] - floor(length(xRange) / 2) - 1) / length(xRange) / (xRange(2) - xRange(1))...
                    * 532e-9 / 2 / 0.42 * 2;
                kyRange = ([1:length(rotYRange)] - floor(length(rotYRange) / 2) - 1) / length(rotYRange) / (rotYRange(2) - rotYRange(1))...
                    * 532e-9 / 2 / 0.42 * 2;
                kzRange = ([1:length(rotZRange)] - floor(length(rotZRange) / 2) - 1) / length(rotZRange) / (rotZRange(2) - rotZRange(1))...
                    * 532e-9 / 2 / 0.42 * 2;
                [kY,kX] = meshgrid(kyRange,kxRange);
                
            % Generate 10% band filters from 0% to 100% of diffraction limit
                fourierFilter = zeros([length(kxRange),length(kyRange),10],'single');
                for nIdx = 1:10
                    fourierFilter(:,:,nIdx) = (sqrt(kY.^2 + kX.^2) > (nIdx - 1) / 10) & (sqrt(kY.^2 + kX.^2) <= (nIdx) / 10);
                    fourierFilter(:,:,nIdx) = fourierFilter(:,:,nIdx) / sum(sum(fourierFilter(:,:,nIdx)));
                end
                
            % Apply filter and sum Fourier content within
            % (Repeat for each z-plane)
                disp('Calculating Fourier content');
                fourierContent = zeros([10,length(kzRange)],'single');
                for zIdx = 1:length(kzRange)
                    fourierPlane = abs(fftshift(fft2(squeeze(restoredDataCube(:,:,zIdx)))));
                    for nIdx = 1:size(fourierContent,1)
                        fourierContent(nIdx,zIdx) = sum(sum(fourierPlane.*squeeze(fourierFilter(:,:,nIdx))));
                    end
                    % Notify every 100 iterations
                        if mod(zIdx,100) == 0
                            disp(strcat('Calculating Fourier content for plane (',num2str(zIdx),') of(',num2str(length(kzRange)),')'));
                        end
                end
                disp('Fourier content calculated');
                
            % Save results
                if typeIdx == 1
                    results.Gaussian.fileName = fileName;
                    results.Gaussian.fullProjection = fullProjection;
                    results.Gaussian.ROIProjection = ROIProjection;
                    results.Gaussian.fourierContent = fourierContent;
                    results.Gaussian.xRange = xRange;
                    results.Gaussian.yRange = yRange;
                    results.Gaussian.zRange = zRange;
                    results.Gaussian.rotYRange = rotYRange;
                    results.Gaussian.rotZRange = rotZRange;
                    results.Gaussian.kxRange = kxRange;
                    results.Gaussian.kyRange = kyRange;
                    results.Gaussian.kzRange = kzRange;
                elseif typeIdx == 2
                    results.Airy.fileName = fileName;
                    results.Airy.fullProjection = fullProjection;
                    results.Airy.ROIProjection = ROIProjection;
                    results.Airy.fourierContent = fourierContent;
                    results.Airy.xRange = xRange;
                    results.Airy.yRange = yRange;
                    results.Airy.zRange = zRange;
                    results.Airy.rotYRange = rotYRange;
                    results.Airy.rotZRange = rotZRange;
                    results.Airy.kxRange = kxRange;
                    results.Airy.kyRange = kyRange;
                    results.Airy.kzRange = kzRange;
                end
            % Clear "restoredDataCube" from memory
                clear restoredDataCube;
        end
        mkdir(OutputFolder)
        save(strcat(OutputFolder,'/results_',num2str(fileIdx),'.mat'),'results');
    end
end