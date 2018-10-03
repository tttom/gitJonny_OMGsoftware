function generateBeamCrossSections

beamTypes = {'Gaussian','Bessel','Airy'};
% beamTypes = {'Bessel'};

lambda = 488e-9; % [metres]
NA = 0.4;
ref_index = 1;
illumination_mag = 40;
illumination_tube_length = 0.2; % [metres]

xR = 2.4e-6; % rayleigh range

bitDepth = 8;

for bIdx = 1:length(beamTypes)
    beamType = beamTypes{bIdx};

    switch beamType
        case 'Gaussian'
            pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
        case 'Bessel'
            beta = 0.1;
            pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* (sqrt(U.^2 + V.^2) >= (1 - beta));
        case 'Airy'
            alpha = 2;
            pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) .* exp(2* pi * 1i * alpha .* (U.^3 + V.^3));
        otherwise
            pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1);
    end




    yRange = [-7.5:0.01:12.5] * 1e-6; % [metres]
    zRange = [-7.5:0.01:12.5] * 1e-6; % [metres]
%     xRange = [-2 * xR, -1 * xR, -0.5 * xR, 0, 0.5 * xR, 1 * xR, 2 * xR, 10 * xR]; % propagation axis [metres]
    xRange = [0, 0.5 * xR, 1 * xR, 2 * xR, 5 * xR]; % propagation axis [metres]

    psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
        ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
        ,NA,ref_index,illumination_mag,illumination_tube_length,[]);

    psf = psf / max(psf(:)) .* (2^bitDepth - 1);
    
    figure;
    for xIdx = 1:length(xRange)
        subplot(2,length(xRange),xIdx);
        image(yRange * 1e6,zRange * 1e6,squeeze(psf(:,:,xIdx)));axis image;
        xlabel('z [um]');
        ylabel('y [um]');
        switch beamType
            case 'Airy'
                xlim([-7.5 7.5] + 5);
                ylim([-7.5 7.5] + 5);
            otherwise
                xlim([-7.5 7.5]);
                ylim([-7.5 7.5]);
        end
        title(strcat('x = ',num2str(xRange(xIdx)/xR),'xR'));
        colormap(jet(2^bitDepth - 1));
        drawnow;shg;
    end
    
    
    yRange = [-7.5:0.01:12.5] * 1e-6; % [metres]
    zRange = [0] * 1e-6; % [metres]
%     xRange = [-1.1 * xR:0.05 * 1e-6:7.6 * xR]; % propagation axis [metres]
    xRange = [-7.6 * xR:0.05 * 1e-6:7.6 * xR]; % propagation axis [metres]
    
    clear psf;
    
    psf = calcVectorialPsf(yRange,zRange,xRange,lambda...
        ,@(U,V) pupilFunctor(U,V) / sqrt(2),@(U,V) 1i * pupilFunctor(U,V) / sqrt(2)...
        ,NA,ref_index,illumination_mag,illumination_tube_length,[]);

    psf = psf / max(psf(:)) .* (2^bitDepth - 1);
    
    
    subplot(2,length(xRange),[(length(xRange) + 1) (length(xRange) * 2)]);
        image(xRange * 1e6, yRange * 1e6,squeeze(psf));axis image;
        xlabel('x [um]');
        ylabel('y [um]');
%         xlim([(-1 * xR) (7.5 * xR)] * 1e6);
        xlim([(-7.5 * xR) (7.5 * xR)] * 1e6);
        switch beamType
            case 'Airy'
                ylim([-7.5 7.5] + 5);
            otherwise
                ylim([-7.5 7.5]);
        end
        colormap(jet(2^bitDepth - 1));
        drawnow;shg;    
end


end