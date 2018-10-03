%%% Function:           anonymousPolynomial
%%% Author:             Jonathan Nylk
%%% Created:            02/05/2018
%%% Description:        Generates an anonymous function form of a
%%%                     zernike polynomial function with given coefficients.
%%%
%%% Inputs:             
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function anonZernike = anonymousZernike(horizCoord_string,vertCoord_string,n,m)

    if nargin < 4
        m = 0;
    end
    if nargin < 3
        n = 2;
    end
    if nargin < 2
        vertCoord_string = 'V';
    end
    if nargin < 1
        horizCoord_string = 'U';
    end
    
    % Check n is non-negative
    if ~(n >= 0)
        disp('n must be non-negative. Setting n = 0.');
        n = 0;
    end
    
    % Check n,m are integers
    if mod(n,1) ~= 0
        disp('n must be an integer. Setting n = 0.');
        n = 0;
    end
    if mod(m,1) ~= 0
        disp('n must be an integer. Setting m = 0.');
        m = 0;
    end
    
    % Check n,m are both odd or both even
    if mod(n,2) == 0
        if mod(m,2) ~= 0
            disp('n is even. m must be even. Setting m = 0.');
            m = 0;
        end
    else
        if mod(m,2) == 0
            disp('n is odd. m must be odd. Setting m = 1.');
            m = 1;
        end
    end
    
    % Check n > m
    if abs(m) > n
        if mod(n,2) == 0
            disp('Cannot have m > n. n is even. Setting m = 0.');
            m = 0;
        else
            disp('Cannot have m > n. n is odd. Setting m = 1.');
            m = 1;
        end
    end
    
    % generate base string (1)
    anonZernike_string = strcat('@(',horizCoord_string,',',vertCoord_string,') 1');
    
    % include azimuthal dependence
%     if mod(n,2) == 0
    if m >= 0
        anonZernike_string = strcat(anonZernike_string,' .* cos('...
            ,num2str(m),'.* atan2(',vertCoord_string,',',horizCoord_string,')) .* (0');
    else
        anonZernike_string = strcat(anonZernike_string,' .* sin('...
            ,num2str(m),'.* atan2(',vertCoord_string,',',horizCoord_string,')) .* (0');
    end
    % generate all radial terms in the n,m-th Zernike polynomial
    for k = 0: (n - abs(m)) / 2
        anonZernike_string = strcat(anonZernike_string,' + '...
            ,num2str((-1).^(k) .* factorial(n - k) ./ factorial(k) ./ factorial(((n + m) ./ 2) - k) ./ factorial(((n - m) ./ 2) - k))...
            ,' .* sqrt(',horizCoord_string,'.^2 + ',vertCoord_string,'.^2).^',num2str(n - (2 .* k)));
    end
    
    % close string and apodize with unit circle
    anonZernike_string = strcat(anonZernike_string,') .* (sqrt(',horizCoord_string,'.^2 + ',vertCoord_string,'.^2) <= 1)');
    
    anonZernike = str2func(anonZernike_string);

    

    
    
end