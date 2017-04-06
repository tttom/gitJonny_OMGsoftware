% deconvolvedImg = deconvolveCubic(img,pixelPitchInSample,wavelength,objectiveNumericalAperture,refractiveIndexOfSample,alpha)
%
% Example cubic phase mask deconvolution in 2D.
% 
% img: 2D matrix of monochrome intensity values
% pixelPitch: the (vertical, horizontal) distance between camera pixels,
%      back projected into the sample.
% objectiveNumericalAperture: the numerical aperture of the detection objective
% refractiveIndexOfSample: default water: 1.33
% alpha: the number of wavelengths of cubic modulation from tbe center of the aperture to
%      the edge of the aperture in one dimension.
%
function deconvolvedImg = deconvolveCubic(img,pixelPitchInSample,wavelength,objectiveNumericalAperture,refractiveIndexOfSample,alpha)
    if nargin<1 || isempty(img)
%         img=readDataCubeFromPngFolder('G:\Data_28-07-15\Beads 500ms 1um\2015-07-28_16-59-06.png');
        img=imread('G:\Data_28-07-15\Beads 500ms 1um\2015-07-28_16-59-06.png');
        img=single(img([1:512]+256,[1:512]+200))./256;
    end
    if nargin<2 || isempty(pixelPitchInSample),
        pixelPitchInSample = [1 1]*6.5e-6./40;
    end
    if numel(pixelPitchInSample)<2,
        pixelPitchInSample(2)=pixelPitchInSample(1);
    end
    if nargin<3 || isempty(wavelength),
        wavelength=532e-9;
    end
    if nargin<4 || isempty(objectiveNumericalAperture),
        objectiveNumericalAperture=0.80;
    end
    if nargin<5 || isempty(refractiveIndexOfSample),
        refractiveIndexOfSample=1.33;
    end
    if nargin<6 || isempty(alpha),
        alpha=-20;
    end
    
    kSNR=5.0;
    
    pupilFunctor=@(U,V) exp(alpha*2i*pi*(U.^3+V.^3)); % cubic
    %pupilFunctor=@(U,V) exp(alpha*2i*pi*(U.^3+V.^3-3*(U.^2.*V+U.*V.^2))); % gen. cubic
    
    % define the sample grid
    imgSize=size(img); imgSize(end+1:3)=1;
    xRange=pixelPitchInSample(1)*([1:imgSize(1)]-floor(imgSize(1)./2)-1);
    yRange=pixelPitchInSample(2)*([1:imgSize(2)]-floor(imgSize(2)./2)-1);
    
    % determine the PSF and the corresponding Wiener filter
    psf=calcVectorialPsf(xRange,yRange,0,wavelength,...
        pupilFunctor,@(U,V) 1i*pupilFunctor(U,V),... % circular polarization
        objectiveNumericalAperture,refractiveIndexOfSample);
    otf=fft2(ifftshift(psf));
    cutOffSpatialFrequency=2*objectiveNumericalAperture./wavelength;
    [~,~,fRel]=calcOtfGridFromSampleFrequencies(1./pixelPitchInSample,imgSize,cutOffSpatialFrequency);
    NSR=ifftshift(fRel)./kSNR;
    clear fRel;
    filter=conj(otf)./(abs(otf).^2+NSR.^2); % Wiener filter
    
    % Deconvolve 2D
    deconvolvedImg=fft2(img); % in place calculation to save memory
    for zIdx=1:imgSize(3), % loop to save a bit more memory
        deconvolvedImg(:,:,zIdx) = ifft2(deconvolvedImg(:,:,zIdx).*filter,'symmetric');
    end
    
    if nargout<1,
        close all;
        
        figure('Position',[100 100 800 600]);
        for zIdx=1:imgSize(3),
            axs(1)=subplot(1,2,1);
            imagesc(yRange*1e6,xRange*1e6,max(0,img(:,:,zIdx))); axis equal tight;
            xlabel('x [\mum]'); ylabel('y [\mum]');
            colorbar();
            axs(2)=subplot(1,2,2);
            imagesc(yRange*1e6,xRange*1e6,max(0,deconvolvedImg(:,:,zIdx)));
            xlabel('x [\mum]'); ylabel('y [\mum]'); axis equal tight;
            colorbar();

            linkaxes(axs);
            
            drawnow();
            
            pause(.1);
        end
        
        clear deconvolvedImg;
    end
end