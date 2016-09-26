function results = pat_regress(y,X)
    
    % Fit regression model to calcium signal.
    %
    % USAGE: results = pat_regress(y,X)
    %
    % INPUTS:
    %   y - [T x 1] time series of neural data
    %   X - [T x D] design matrix, where each column corresponds to a regressor (see pat_design.m)
    %
    % OUTPUTS:
    %   results - structure containing the following fields:
    %               .b - [D x 1] vector of maximum likelihood coefficient estimates
    %               .bint - [D x 2] 95% confidence intervals
    %               .stats - array of regression statistics (see regress.m)
    %               .s2 - maximum likelihood estimate of noise variance
    %               .loglik - log likelihood
    %               .bic - Bayesian information criterion
    %
    % Sam Gershman, Sep 2016
    
    [results.b,results.bint,res,~,results.stats] = regress(y,X);
    results.yhat = X*results.b;
    results.s2 = mean(res.^2);
    results.loglik = sum(log(normpdf(y,results.yhat,sqrt(results.s2))));
    results.bic = -2*results.loglik + length(results.b)*log(length(y));