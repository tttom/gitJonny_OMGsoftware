%%% Function:           anonymousPolynomial
%%% Author:             Jonathan Nylk
%%% Created:            14/03/2018
%%% Description:        Generates an anonymous function form of a
%%%                     polynomial function with given coefficient.
%%%
%%% Inputs:             'var_string':
%%%                         A string containing the variable name for the
%%%                         anonymous function (default 'x').
%%%                     'polyCoeffs':
%%%                         For a n-th order polynomial, this is a size n+1
%%%                         vector containin the polynomial coefficients in
%%%                         descending order (x^n,x^n-1,...n^0).
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function anonPoly = anonymousPolynomialAsymmetric(var_string,polyCoeffs)

    if nargin < 2
        polyCoeffs = [7,0,0,0];
    end
    if nargin < 1
        var_string = 'x';
    end
    
    % generate base function (zero)
    anonPoly_string = strcat('@(',var_string,') 0');
    
    %add (n-1)th term to function
    for polyIdx = 1:length(polyCoeffs)
        anonPoly_string = strcat(anonPoly_string,'+'...
            ,num2str(polyCoeffs(polyIdx)),'.*'...
            ,'sign(',var_string,') .* abs('...
            ,var_string,'.^',num2str(length(polyCoeffs) - polyIdx),')');
    end
    
    anonPoly = str2func(anonPoly_string);
    
end