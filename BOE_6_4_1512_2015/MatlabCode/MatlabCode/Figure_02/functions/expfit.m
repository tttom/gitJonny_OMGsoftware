function [beta, resid, COV] = expfit(model, x, y, fPar)
% CREATED: Søren Gammelmark, Aarhus University, July 2011
% MODIFIED: Thue Bjerring Lindballe, Aarhus University, 2011
% MODIFIED: Martin Verner Gammelgaard Kristensen, University of St Andrews, October 2013

    %% Exponential fitting
    beta    = fminsearch(@(alpha) Likelihood(x, y, alpha), fPar); 
        % minimizing the output of Likelihood, i.e. testing/finding values
        % for the fitting parameters, starting out with the values in beta0
    
    N       = length(fPar);         % Number of parameters to fit
    I       = zeros([N, N]);        % Predefining for faster computation
    h       = 1e-5;                 % Small step for gradient calculation
    basis   = eye(N);               % Used in the computation of the covariance matrix

    grad    = zeros([N, length(x)]);
    for i = 1:N
        grad(i, :)      = (model(beta + h * basis(i, :), x) - model(beta, x)) / h;
    end     % Calculating the gradient
    
    lambda  = model(beta, x)';
    for i = 1:N
        for j=1:N
            I(i, j)     = sum(grad(i, :) .* grad(j, :) ./ (lambda .^ 2) );
        end % Calculating the inverse of the covariance matrix
    end
    
    COV     = inv(I);
    resid   = model(beta, x) - y;  	% Difference between the minimization result and data
    
    function L = Likelihood(x, y, alpha)
        lambda  = model(alpha, x);	% Model/theory
        L       = sum(log(lambda) + y ./ lambda );
            % Sum over the logarithm to the difference between model and
            % data weighted with the known (exponential distribution*) 
            % uncertainty = the value of the model
            % See "Least Squares as a Maximum Likelihood Estimator"
            % in Numerical Recipes by William Press et. al
            % * See Kirstine Berg-Sorensen et. al, Rev. Sci. Instrum. 2004
            % Probability of one measurement (value of "spectral density", 
            % P at a given frequency) in power spectrum: 
            %        p_i(P,lam_i)=1/lam_i * exp(-P/lam_i)
            % Probability of dataset: p(P,lam_i)=prod_i(p_i) gives L above 
            % when taking the logarithm (to make values larger and give a more
            % computationally friendly expresion)... Note that mean(p_i)=1/lambda_i
    end

end