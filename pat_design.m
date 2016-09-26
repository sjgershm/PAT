function [X, name] = pat_design(data,opts)
    
    % Construct design matrix.
    %
    % USAGE: [X, name] = pat_design(data,[opts])
    %
    % INPUTS:
    %   data - [1 x S] structure, where each element corresponds to a session. Each element has the following fields:
    %           .y - [T x 1] neural data
    %           .events - [1 x nEvents] structure where each element corresponds to an event type, with the following fields:
    %                     .ons - vector of onsets (seconds relative to session start time)
    %                     .name - string specifying event name
    %                     .pmod - [1 x nPmod] structure where each element corresponds to a parametric modulator for that event, with the following fields:
    %                             .param - vector of parameter values (same length as ons)
    %                             .name - name of parametric modulator
    %   opts (optional) - options structure (see pat_opts.m)
    %
    % OUTPUTS:
    %   X - [T x D] design matrix
    %   name - [1 x D] cell array of regressor names
    %
    % Sam Gershman, Sep 2016
    
    % make sure all options are defined
    if nargin < 2; opts = []; end
    opts = pat_opts(opts);
    
    S = length(data);
    X = [];
    for s = 1:S
        events = data(s).events;
        
        % if no data field, make dummy data
        if ~isfield(data,'y')
            maxons = 0;
            for i = 1:length(data(s).events)
                maxons = max(maxons,max(data(s).events(i).ons));
            end
            T = opts.samprate*(maxons + opts.maxt);
        else
            T = length(data(s).y);
        end
        
        % construct timepoints
        L = T/opts.samprate;
        t = linspace(0,L,T);
        
        % get calcium response function
        crf = pat_crf(opts.p,opts.samprate,opts.maxt);
        
        % build design matrix
        d = 0; Xs = [];
        for i = 1:length(events)
            d = d + 1;
            x = zeros(T,1);
            k = zeros(length(events(i).ons),1);
            for j = 1:length(events(i).ons)
                [~,k(j)] = min(abs(events(i).ons(j)-t));
            end
            x(k) = 1;
            name{d} = events(i).name;
            x = conv(x,crf); % convolve regressor with calcium response function
            Xe = x(1:T);
            
            % add parametric modulators
            if isfield(events(i),'pmod') && ~isempty(events(i).pmod)
                for j = 1:length(events(i).pmod)
                    d = d + 1;
                    x(k) = events(i).pmod(j).param;
                    name{d} = [events(i).name,' * ',events(i).pmod(j).name];
                    x = conv(x,crf); % convolve regressor with calcium response function
                    Xe(:,j+1) = x(1:T);
                end
            end
            
            % orthogonalize parametric modulators with respect to the event
            Xe = orthogonalize(Xe);
            
            % add to design matrix
            Xs = [Xs Xe];
        end
        
        % add intercept
        Xi = zeros(T,S);
        Xi(:,s) = 1;
        
        X = [X; Xi Xs];
    end
    
end

function X = orthogonalize(X)
    % Recursive Gram-Schmidt orthogonalisation 
    
    sw    = warning('off','all');
    [n,m] = size(X);
    X     = X(:, any(X));
    rankX = rank(full(X));
    try
        x     = X(:,1);
        j     = 1;
        for i = 2:size(X, 2)
            D = X(:,i);
            D = D - x*(pinv(x)*D);
            if norm(D,1) > exp(-32)
                x          = [x D];
                j(end + 1) = i;
            end
            if numel(j) == rankX, break, end
        end
    catch
        x     = zeros(n,0);
        j     = [];
    end
    warning(sw);
    X      = zeros(n,m);
    X(:,j) = x;
end