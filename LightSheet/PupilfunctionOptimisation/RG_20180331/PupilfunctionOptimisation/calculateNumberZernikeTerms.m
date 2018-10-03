function [ns,ms,number_zernike_terms] = calculateNumberZernikeTerms(n_max)

    if nargin < 1
        n_max = 10;
    end

    % initialise term counter
    number_zernike_terms = 0;

    % permute over all allowed n,m combinations
    for n = 0:n_max
        m = -1 * n;
        while abs(m) <= n
            number_zernike_terms = number_zernike_terms + 1;
            m = m + 2;
        end
    end

    % initialise ns and ms arrays
    ns = zeros([1,number_zernike_terms]);
    ms = zeros([1,number_zernike_terms]);

    % re-initialise term counter
    number_zernike_terms = 0;

    for n = 0:n_max
        m = -1 * n;
        while abs(m) <= n
            number_zernike_terms = number_zernike_terms + 1;
            ns(number_zernike_terms) = n;
            ms(number_zernike_terms) = m;
            m = m + 2;
        end
    end

end
