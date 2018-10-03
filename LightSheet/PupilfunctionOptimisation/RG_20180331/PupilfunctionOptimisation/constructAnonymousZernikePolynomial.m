function [fullZernikePolynomial] = constructAnonymousZernikePolynomial(anonZernikeTerms,coeffs)

    if nargin < 1
        anonZernikeTerms = generateAnonymousZernikeTerms([0,1,1],[0,-1,1],'U','V');
%         anonZernikeTerms = generateAnonymousZernikeTerms([1],[-1],'U','V');
    end

    if nargin < 2
        coeffRangeMin = -10;
        coeffRangeMax = 10;
        coeffs = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* rand(size(anonZernikeTerms));
    end
    
    % check dimensions match
    if size(anonZernikeTerms) ~= size(coeffs)
        disp('Size of "anonZernikeTerms" and "coeffs" do not match. Changing "coeffs" to random numbers in range: -10:10.');
        coeffRangeMin = -10;
        coeffRangeMax = 10;
        coeffs = coeffRangeMin + (coeffRangeMax - coeffRangeMin) .* rand(size(anonZernikeTerms));
    end
    
    %initialise polynomial
    fullZernikePolynomial = @(U,V) 0;
    
    % construct anonymous Zernike polynomial
    for p_idx = 1:length(anonZernikeTerms)
        term_p = anonZernikeTerms{p_idx};
        term_p_coeff = coeffs(p_idx);
        fullZernikePolynomial = @(U,V) fullZernikePolynomial(U,V) + term_p_coeff .* term_p(U,V);
    end


end