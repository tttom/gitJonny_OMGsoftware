function [psfField,xRange,yRange,zRange] = polarisationEFieldVectorialCalculations(polarisation_state,wavelength)
% psfSliceIntensity=sum(abs(psfSlice).^2,4); % Intensity from field

    if nargin < 1
        polarisation_state = 'C_LH';
    end

    if nargin < 2
        wavelength = 532e-9;
    end
    
%     xRange = [-1:0.0025:1]*1e-6;
%     yRange = [-1:0.0025:1]*1e-6;
%     zRange = [-1:0.005:1]*1e-6;
    xRange = [-1:0.01:1]*1e-6;
    yRange = [-1:0.01:1]*1e-6;
    zRange = [-2:0.01:2]*1e-6;
    
    objectiveNumericalAperture = 0.9;
    refractiveIndexOfSample = 1;
    objectiveMagnification = 100;
    objectiveTubeLength = 0.2;
    projectionDimensions = [];
    
    pupilFunctor = @(U,V) (sqrt(U.^2 + V.^2) <=1) * 1;
    
    switch polarisation_state
        case 'C_LH'
            % Circular (Left Handed)
            pupilFunctorH = @(U,V) pupilFunctor(U,V) / sqrt(2);
            pupilFunctorV = @(U,V) pupilFunctor(U,V) / sqrt(2) * 1i;
        case 'C_RH'
            % Circular (Reft Handed)
            pupilFunctorH = @(U,V) pupilFunctor(U,V) / sqrt(2);
            pupilFunctorV = @(U,V) pupilFunctor(U,V) / sqrt(2) * -1i;
        case 'L_H'
            % Linear (Horizontal)
            pupilFunctorH = @(U,V) pupilFunctor(U,V);
            pupilFunctorV = @(U,V) 0;
        case 'L_V'
            % Linear (Vertical)
            pupilFunctorH = @(U,V) 0;
            pupilFunctorV = @(U,V) pupilFunctor(U,V);
        case 'L_H+V'
            % Linear (+45deg)
            pupilFunctorH = @(U,V) pupilFunctor(U,V) / sqrt(2);
            pupilFunctorV = @(U,V) pupilFunctor(U,V) / sqrt(2);
        case 'L_H-V'
            % Linear (-45deg)
            pupilFunctorH = @(U,V) pupilFunctor(U,V) / sqrt(2);
            pupilFunctorV = @(U,V) pupilFunctor(U,V) / sqrt(2) * -1;
        otherwise
            disp('Polarisation state not recognised. Resorting to "Circular (Left Handed)" (C_LH).')
            pupilFunctorH = @(U,V) pupilFunctor(U,V) / sqrt(2);
            pupilFunctorV = @(U,V) pupilFunctor(U,V) / sqrt(2) * 1i;
    end

    % Compute 3D proifle then take 2D cross-sections
    tic
    [psf,psfField] = calcVectorialPsf(xRange,yRange,zRange,wavelength...
        ,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture...
        ,refractiveIndexOfSample,objectiveMagnification,objectiveTubeLength...
        ,projectionDimensions);
    toc
    
    plotEFieldCartesianCrossSections(xRange,yRange,zRange,psf,psfField)
    
end
    
