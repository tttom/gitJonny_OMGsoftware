function [A,fval,residual] = BesselCrossSectionDataOptimisation(referenceCrossSection,inputCrossSection)

    %sum of squares of difference between scaled inputs
    optimisationFunctor = @(referenceCrossSection,inputCrossSection,A) sum((referenceCrossSection - A .* inputCrossSection).^2);

    %initial guess for A
    A0 = max(referenceCrossSection) / max(inputCrossSection);
    
    %optimisation (solve for A)
    [A,fval] = fminsearch(@(A) optimisationFunctor(referenceCrossSection,inputCrossSection,A),A0);

    % determine residual
    residual = referenceCrossSection - A .* inputCrossSection;    
    
end