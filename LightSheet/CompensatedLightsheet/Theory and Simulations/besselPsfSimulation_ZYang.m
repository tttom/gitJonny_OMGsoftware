
function [beam_CrossSection,lightSheet_CrossSection,transverse_Range,longitudinal_Range,fullPsf] = besselPsfSimulation_ZYang(cartesianPupilFunction,z)

% define the optical system
lambda = 0.532e-6;
% FL1 = 0.45;             % focal length of the first lens, after the axicon
% FL2 = 0.018;            % focal length of the objective lens
FL2 = 0.0084;            % focal length of the objective lens (gives NA = 0.4203)
NA = 0.4203;
% alpha = 2*pi/360;       % angle of the axicon is 1 degree
% n = 1.4533;             % refractive index of the axicon material
% NA = 0.5;
% w0=2.5e-3;                % radius of the incoming Bessel beam
n_sample = 1.33;         % refractive index of sample medium

% z = [-100:5:100] * 1e-6;	%propagation

if nargin < 1
    
%     S = load('NormalBessel_501x501px.mat','cartesianPupilFunction');
%     S = load('Compensated_Bessel(0.25)_501x501px.mat','cartesianPupilFunction');
    S = load('FlatTopBessel_501x501px.mat','cartesianPupilFunction');


    cartesianPupilFunction = S.cartesianPupilFunction;
    clear S;
    
%     [X,Y] = meshgrid(-250:1:250,-250:1:250);
%     spiral = atan2(Y,X);
%     delta = pi/126;
%     pieSlice = ((abs(spiral) >= (pi/6 - delta)) & (abs(spiral) <= (pi/6 + delta)))...;
% + ((abs(spiral) >= (pi/2 - delta)) & (abs(spiral) <= (pi/2 + delta)))...
% +((abs(spiral) >= (5*pi/6 - delta)) & (abs(spiral) <= (5*pi/6 + delta)));
% 
%     cartesianPupilFunction = cartesianPupilFunction .* pieSlice;
    
end

if nargin < 2
   z = [-100:1:100] * 1e-6; % propagation [metres] 
%    z = [-5:5:5] * 1e-6; % propagation [metres] 
%     z = [-40:20:40] * 1e-6; % propagation [metres] 
end
    

Flag_2D_Simulation = 1;

if ~Flag_2D_Simulation
    fullPsf = 0;
end

%     w0 = nn/2*1e-3;
% %     
% %     % ============= calculate the parameter of ring at back apeture [NOT REQURIED] ===========
%     % based on paper Pierre-Andre Belanger 1978.
%     theoryThickness = 3.3*lambda*FL1/pi/w0;     % Thickness of the ring
%     theoryRadius = (n-1)*alpha*FL1;             % diameter of the ring
%     n1=1;
%     trueNA = sin(atan(theoryRadius/FL2));       % Resulting NA with the diametre of ring
%     beta = theoryThickness/theoryRadius;        % Size of the ring
%     theoryFOV = lambda/n1 / (2*(1-sqrt(1-(trueNA/n1)^2))*beta); % FOV based on Airy paper
%     
    % =========================== define grid  ========================
    W1 = 0.1;               % sample region (m) - Make large to get good resolution in Fourier plane
    M = 5000;               % Half-sampling range - make this an integer times larger than the pupil function exported from Bessel pupil simulation.
    dx1 = W1/M;             % pixel size
    x1 = -W1/2:dx1:W1/2;    % define coordinate
    [X,Y] = meshgrid(x1,x1);
    [~, r1] = cart2pol(X,Y);% radial coords not needed
    
    % ======================== beam at back aperture ====================
%     g = exp(-(r1-theoryRadius).^2/(theoryThickness/2)^2);
%     Ig = abs(g).^2;
    
%     load('NormalBessel_501x501px.mat','cartesianPupilFunction');
%     load('FlatTopBessel_501x501px.mat','cartesianPupilFunction');
%     load('Compensated_Bessel(0.25)_501x501px.mat','cartesianPupilFunction');
%     load('Compensated_Bessel(0.5)_501x501px.mat','cartesianPupilFunction');
%     load('Compensated_Bessel(0.75)_501x501px.mat','cartesianPupilFunction');
%     load('Compensated_Bessel(1.0)_501x501px.mat','cartesianPupilFunction');
    
    cartesianPupilFunction(isnan(cartesianPupilFunction)) = 0;
%     cartesianComplexPupilFunction = cartesianPupilFunction;
%     clear cartesianPupilFunction;
%     cartesianPupilFunction = (real(cartesianComplexPupilFunction) >= 0) .* cartesianComplexPupilFunction...
%         + (real(cartesianComplexPupilFunction) < 0) * -1 .* cartesianComplexPupilFunction * exp(1i * pi);
%     figure;
%     subplot(2,2,1);imagesc(real(cartesianComplexPupilFunction));axis image;
%     subplot(2,2,2);imagesc(imag(cartesianComplexPupilFunction));axis image;
%     subplot(2,2,3);imagesc(real(cartesianPupilFunction));axis image;
%     subplot(2,2,4);imagesc(imag(cartesianPupilFunction));axis image;
    
    
    g = zeros(size(r1));
    g([1:size(cartesianPupilFunction,1)] + (size(g,1) - 1) / 2 - (size(cartesianPupilFunction,1) - 1) / 2 ...
    ,[1:size(cartesianPupilFunction,2)] + (size(g,2) - 1) / 2 - (size(cartesianPupilFunction,2) - 1) / 2)...
        = cartesianPupilFunction / max(cartesianPupilFunction(:));
    
    
    Ig = abs(g).^2;
    
    % ======================= display the original plane ======================
        displayLength = 10e-3;
        pixelNo = displayLength/dx1;
%         figure
%         subplot(2,2,1);
%         displayCoord = M / 2 + 1 - pixelNo / 2:M / 2 + 1 + pixelNo / 2;
%         imagesc(x1(displayCoord) * 1e3, x1(displayCoord) * 1e3,abs(g(displayCoord,displayCoord)));
%         axis image;
%         xlabel('mm');
%         ylabel('mm');
%         title('Pupil function amplitude.');
%         subplot(2,2,3);
%         imagesc(x1(displayCoord) * 1e3, x1(displayCoord) * 1e3,angle(g(displayCoord,displayCoord)));
%         axis image;
%         xlabel('mm');
%         ylabel('mm');
%         title('Pupil function phase.');
    
    
    % ================ FFT to focal plane    =============
    g = fftshift(g);
    g2 = fft2(g);
    g2  = fftshift(g2);
    Ig2 = abs(g2).^2;
    
    dx2 = 1/W1*FL2*lambda;
    x2 = -1/(2*dx1):1/W1:1/(2*dx1);
    x2 = x2 *FL2*lambda;
    y2 = x2;
    W2 = dx2 *M;
    
    % ============== display the focal plane Bessel beam ==================
        [~,maxIdx] = max(Ig2(:));
        [maxIdx,~] = ind2sub(size(Ig2),maxIdx);
    
        displayLength2 = 20e-6;
        pixelNo = floor(displayLength2/dx2/2);
        displayCoord = maxIdx-pixelNo:maxIdx+pixelNo;
%         subplot(2,2,2);
%         imagesc(x2(displayCoord)*1e6, x2(displayCoord)*1e6,...
%             Ig2(displayCoord,displayCoord));
% %         imagesc(Ig2(displayCoord,displayCoord));
%         axis image;
%         xlabel('um');
%         ylabel('um');
%         title('Intensity of the Bessel beam on the focal plane');
%         subplot(2,2,4);
%         plot(x2(displayCoord)*1e6,Ig2(maxIdx,displayCoord));
% %         plot(Ig2(M/2+1,displayCoord));
%         xlabel('um');
%         title('Profile of the  Bessel beam on focal plane');
    
%     %================ propagation simulation ==========================
%     
%     fx2 = -1/(2*dx2):1/W2:1/(2*dx2);    % frequence coords
%     [FX2,FY2] = meshgrid(fx2,fx2);
%     [~,r2] = cart2pol(FX2,FY2);
%     
%     z = 10e-6;                           % propogation step
%     propagationDistance = 500e-6;       % simulation distance
%     
%     H = exp(-1i*pi*lambda*z*r2.^2);     % transfer function
%     H = fftshift(H);
%     G2 = fft2(fftshift(g2));
%     centreLine = g2(M/2+1,:);
%     Psf(1,:) = abs(centreLine).^2;
%     
%     integratedPsf(1,:) = mean(abs(g2).^2);
%     
%     
%     % ==================== propagation ================================
%     for n=1:propagationDistance/z
%         n
%         U2 = H.*G2;
%         u2 = ifftshift(ifft2(U2));
%         centreLine = u2(M/2+1,:);               % Centre line of the field
%         Psf(n+1,:) = abs(centreLine).^2;        % Intensity of the centre Line to PSF
%         integratedPsf(n+1,:) = mean(abs(u2).^2);
%         G2 = U2;
%     end

    %================ propagation simulation ==========================
    
    fx2 = -1/(2*dx2):1/W2:1/(2*dx2);    % frequence coords
    [FX2,FY2] = meshgrid(fx2,fx2);
    [~,r2] = cart2pol(FX2,FY2);
    
%     z = [-100:5:100] * 1e-6;	%propagation
    
    Psf = zeros(size(z,2),size(r2,2));
    integratedPsf = zeros(size(Psf));
    
    if Flag_2D_Simulation
        fullPsf = zeros(size(r2,1),size(r2,1),size(z,2));
    end
    
    for n=1:length(z)
%         n
        H = exp(-1i * pi * lambda / n_sample * z(n) * r2.^2);     % transfer function
        H = fftshift(H);
        U2 = H.*g;
        u2 = ifftshift(ifft2(U2));
        centreLine = u2(maxIdx,:);               % Centre line of the field
        Psf(n,:) = abs(centreLine).^2;        % Intensity of the centre Line to PSF
        integratedPsf(n,:) = sum(abs(u2).^2,1);
        if Flag_2D_Simulation
            fullPsf(:,:,n) = abs(u2).^2;
        end
        G2 = U2;
    end

    % normalisation
    Psf = Psf / max(Psf(:));
    integratedPsf = integratedPsf / max(integratedPsf(:));
    
    beam_CrossSection = Psf;
    lightSheet_CrossSection = integratedPsf;
    
    transverse_Range = x2;
    longitudinal_Range = z;
    
%     C_abs = 64.95 * 100;   % [metres^-1]
%     absorption_decay = repmat(exp(-C_abs * z).',[1 length(x2)]);
%     
%     psf = Psf;
%     % normalise
%     psf = psf ./ max(psf(:));
%     % take a copy of the PSF before applying absorption
%     psf_before_absorption = psf;
%     % simulate "perfect" absorption
%     psf = psf .* absorption_decay;
%     % re-normalise
%     psf = psf ./ max(psf(:));
%     Psf = psf;
%     
%     integratedpsf = integratedPsf;
%     % normalise
%     integratedpsf = integratedpsf ./ max(integratedpsf(:));
%     % take a copy of the PSF before applying absorption
%     integratedpsf_before_absorption = integratedpsf;
%     % simulate "perfect" absorption
%     integratedpsf = integratedpsf .* absorption_decay;
%     % re-normalise
%     integratedpsf = integratedpsf ./ max(integratedpsf(:));
%     integratedPsf = integratedpsf;

% 
%     
%     %     normPsf = normr(Psf);
% %     Psf2F = Psf.^2;
% %     Psf3F = Psf.^3;
%     displayLength2 = 100e-6;
%     pixelNo = floor(displayLength2/dx2/2);
%     displayCoord = M/2+1-pixelNo:M/2+1+pixelNo;
%     
%     figure;
% %     imagesc([0:z:displayX]*1e6,x2(displayCoord)*1e6,integratedPsf(1:(displayX/z),displayCoord)')
%     imagesc(z*1e6,x2(displayCoord)*1e6,psf_before_absorption(:,displayCoord)')
%     axis image;
%     xlabel('\mum');
%     ylabel('\mum');
%     title('Single photon Intensity of the PSF - centre line');
%     
%     figure;
% %     imagesc([0:z:displayX]*1e6,x2(displayCoord)*1e6,integratedPsf(1:(displayX/z),displayCoord)')
%     imagesc(z*1e6,x2(displayCoord)*1e6,Psf(:,displayCoord)')
%     axis image;
%     xlabel('\mum');
%     ylabel('\mum');
%     title('Single photon Intensity of the PSF - centre line');
% %     subplot(3,3,3)
% %     plot(x2(displayCoord)*1e6,integratedPsf(1,displayCoord)/max(integratedPsf(:)));
% %     view(90,90) % swap the x and y axis
% %     ylabel('Normalized intensity');
% %     xlabel('\mum');
% %     title('Intensity profile at 0 (focus)');
%     
%     figure;
% %     imagesc([0:z:displayX]*1e6,x2(displayCoord)*1e6,integratedPsf(1:(displayX/z),displayCoord)')
%     imagesc(z*1e6,x2(displayCoord)*1e6,integratedpsf_before_absorption(:,displayCoord)')
%     axis image;
%     xlabel('\mum');
%     ylabel('\mum');
%     title('Single photon Intensity of the PSF - integrated');
% 
%     figure;
% %     imagesc([0:z:displayX]*1e6,x2(displayCoord)*1e6,integratedPsf(1:(displayX/z),displayCoord)')
%     imagesc(z*1e6,x2(displayCoord)*1e6,integratedPsf(:,displayCoord)')
%     axis image;
%     xlabel('\mum');
%     ylabel('\mum');
%     title('Single photon Intensity of the PSF - integrated');
%     
%     [~,z_zero] = max(z == 0);
%     [~,maxIdx_I] = max(Psf(z_zero,:));
%     figure;
%     plot(z*1e6,Psf(:,maxIdx_I) / max(Psf(:,maxIdx_I)),z*1e6,integratedPsf(:,maxIdx_I) / max(integratedPsf(:,maxIdx_I)));
%     xlabel('\mum');
%     ylabel('Intensity [a.u.]');
%     title('On axis intensity (beam/sheet)');
%     
%     
%     %%% MTF plots (partly for calibration)
%     % code switching to microscope coordinates: [x: propagation, y, z: axial]
%     zRange = x2;
%     xRange = z;
% %     kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end); % not normalised coords
%     kzRange = ([1:length(zRange)] - floor(length(zRange) / 2) - 1) * 1 / 2 / zRange(end) / 2 / NA * lambda;
%     kzRange = kzRange(ceil(length(kzRange) / 2):end);
%     
%     MTF_z = abs(fftshift(fft(Psf,[],2),2));
%     % normalisation
%     MTF_z = MTF_z ./ max(MTF_z(:));
%     MTF_z = MTF_z(:,ceil(size(MTF_z,2) / 2):end); % range, kz = [0:1]
%     
%     figure;
%     plot(kzRange,MTF_z(floor(length(xRange) / 2) + 1,:));
%     xlabel('k_z');
%     ylabel('Magnitude [a.u.]');
%     title('MTF (x = 0)');
%     xlim([0 1]);
% %     
% %     figure;
% %     imagesc(xRange * 1e6,kzRange,MTF_z.');
% %     ylim([0 1]);
% %     xlabel('x [um]');
% %     ylabel('k_z');
% %     
% %     figure;
% %     imagesc(xRange * 1e6,kzRange,flipud(MTF_z >=0.05).');
% %     ylim([0 1]);
% %     colormap gray;
% %     xlabel('x [um]');
% %     ylabel('k_z');
%     
%     
%     
% %     figure;
% %     imagesc(z*1e6,x2(displayCoord)*1e6,Psf2F(:,displayCoord)')
% % %     imagesc(Psf2F(-1:(displayX/z),displayCoord)')
% %     axis image;
% %     xlabel('\mum');
% %     ylabel('\mum');
% %     title('2 photon Intensity of the PSF - integrated');
% % %     subplot(3,3,6)
% % %     plot(x2(displayCoord)*1e6,Psf2F(1,displayCoord)/max(Psf2F(:)));
% % %     view(90,90) % swap the x and y axis
% % %     ylabel('Normalized intensity');
% % %     xlabel('\mum');
% % %     title('Intensity profile at 0 (focus)');
% %     figure;
% %     imagesc(z*1e6,x2(displayCoord)*1e6,Psf3F(:,displayCoord)')
% % %     imagesc(Psf3F(1:(displayX/z),displayCoord)')
% %     axis image;
% %     xlabel('\mum');
% %     ylabel('\mum');
% %     title('3 photon Intensity of the PSF - integrated');
% % %     subplot(3,3,9)
% %     plot(x2(displayCoord)*1e6,Psf3F(1,displayCoord)/max(Psf3F(:)));
% % %     view(90,90) % swap the x and y axis
% %     ylabel('Normalized intensity');
% %     xlabel('\mum');
% %     title('Intensity profile at 0 (focus)');
%     
%     
%     
%     %     fileName = ['psf',num2str(w0*1e3),'mm'];
%     %     save(strcat(fileName,'.mat'),'Psf','dx2','M','z','propagationDistance','x2');
%     %     saveas(h,strcat(fileName,'.svg'),'svg');
%     %     saveas(h,strcat(fileName,'.png'),'png')
end

