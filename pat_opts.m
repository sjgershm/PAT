function opts = pat_opts(opts)
    
    % Set default options.
    %
    % USAGE: opts = pat_opts([opts])
    %
    % INPUTS:
    %   opts - structure with any of the following fields (missing or empty fields will
    %   be set to defaults):
    %           .samprate - sampling rate (hz; default: 1000)
    %           .maxt - max time (seconds; default: 10)
    %           .p - parameters of the generalized extreme value distribution (default: [0.82 1.2 1.45])
    %
    % OUTPUTS:
    %   opts - completed options structure
    %
    % Sam Gershman, Sep 2016
    
    def_opts.samprate = 1000;
    def_opts.maxt = 7;
    def_opts.p = [0.82 1.2 1.45];
    
    if nargin < 1 || isempty(opts)
        opts = def_opts;
    else
        F = fieldnames(def_opts);
        for j = 1:length(F)
            if ~isfield(opts,F{j})
                opts.(F{j}) = def_opts.(F{j});
            end
        end
    end