function [anonZernikeTerms] = generateAnonymousZernikeTerms(ns,ms,horizCoord_string,vertCoord_string)

    if nargin < 1
        ns = [0,1,1];
    end

    if nargin < 2
        ms = [0,-1,1];
    end
    
    if nargin < 3
        horizCoord_string = 'U';
    end
    
    if nargin < 4
        vertCoord_string = 'V';
    end

    %initialise anonymous function cell
    anonZernikeTerms = cell([1,length(ns)]);
    
    % generate anonymous Zernike polynomial terms
    for p_idx = 1:length(ns)
        anonZernikeTerms{p_idx} = anonymousZernike(horizCoord_string,vertCoord_string,ns(p_idx),ms(p_idx));
    end

end
    