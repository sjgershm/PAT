function p = pat_fit_crf(y,opts)
    
    % Fit calcium response function.
    %
    % USAGE: p = pat_fit_crf(y,opts)
    %
    % INPUTS:
    %   y - [T x 1] neural response to an event
    %   opts - options structure
    %
    % OUTPUTS:
    %   p - [1 x 3] parameters of the generalized extreme value distribution
    %
    % Sam Gershman, Sep 2016
    
    if nargin < 2; opts = []; end
    opts = pat_opts(opts);
    y = y';
    fun = @(p) sum((y-pat_crf(p,opts.samprate,opts.maxt)).^2);
    p = fminsearch(fun,opts.p);