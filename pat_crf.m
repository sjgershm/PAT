function crf = pat_crf(p,samprate,maxt)
    
    % Compute calcium response function, modeled using the generalized extreme value PDF.
    %
    % USAGE: crf = pat_crf(p,samprate,maxt)
    %
    % INPUTS:
    %   p - [1 x 3] parameters
    %   samprate - sampling rate (hz)
    %   maxt - maximum time (seconds)
    %
    % OUTPUTS:
    %   crf - calcium response function
    %
    % Sam Gershman, Sep 2016
    
    u = linspace(0,10,maxt*samprate);
    crf = gevpdf(u,p(1),p(2),p(3));
    crf = crf/sum(crf);