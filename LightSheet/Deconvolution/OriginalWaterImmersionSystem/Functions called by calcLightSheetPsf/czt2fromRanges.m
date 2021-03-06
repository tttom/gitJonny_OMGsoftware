function f=czt2fromRanges(x,xRange,yRange)
% Prepares data from calcVectorialPsf for calculation of its partial
% spectrum, ultimately using the built-in function czt (2D chirp z transform).
% The c.z.t. is the same as f=fftshift(fft2(ifftshift(x)))*sqrt(prod(size(x))) for:
% xRange=[0:size(x,1)-1]-floor(size(x,1)/2);yRange=[0:size(x,2)-1]-floor(size(x,2)/2);
% The L2-norm is scaled so that it would equal the L2-norm of the input for large output sizes.
%
% The units of x/yRange are thus twice the highest possible spatial frequency
% of the field PSF when the pupil disc fits the sample area of input x.
% Nyquist sampling is done for x/yRange=[from:0.5:to];
%
% x should not be ifftshifted, implicit zero padding to the right!
    
    inputSize=size(x);
    outputRanges={xRange,yRange};
    nbRanges=length(outputRanges);
    
    deltaRng=0.5*ones(1,nbRanges);
    M=zeros(1,nbRanges); A=M;
    for (dimIdx=1:nbRanges)
        rng=outputRanges{dimIdx};
        M(dimIdx)=length(rng); % The output length of the transform
        A(dimIdx)=exp(2i*pi*rng(1+floor(end/2))/inputSize(dimIdx)); % Offset on contour (fftshifted, hence~/2)
        if (M(dimIdx)>1)
            deltaRng(dimIdx)=diff(rng(1:2));
        end
    end
    deltaOmegaInCycles=deltaRng./inputSize(1:2);
    W=exp(-2i*pi*deltaOmegaInCycles); % The ratio between the points on the spiral contour
    
    if (any(deltaOmegaInCycles.*M>1))
        logMessage('Warning: undersampling pupil by a factor of (%0.3f, %0.3f). Image replication will occur! Reduce the number of pixels in the lateral dimension.',deltaOmegaInCycles.*M);
    end
    
    %Call function cztn to decompose x into vectorial components for
    %separate chirped z transforms.
    f=cztn(x,M,W,A,'centered'); 
    
    % The energy is spread over more sample points when sub-sampling. The next
    % line makes sure that the input and output norm are the same when the output space size M goes to infinity.
    f=f*prod(sqrt(deltaOmegaInCycles(deltaOmegaInCycles~=0)));
end
