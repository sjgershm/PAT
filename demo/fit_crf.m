function p = fit_crf
    
    % Fit calcium response function GCaMP data from free reward trials (data courtesy of HyungGoo Kim).
    %
    % USAGE: p = fit_crf
    %
    % p - estimated parameters of the CRF.
    %
    % Sam Gershman, Sep 2016
    
    load gcamp_data1
    ix = (find(x==0)+1):size(g,2);
    y = mean(g(:,ix)); y = y-y(1); y = y./sum(y);
    opts = pat_opts;
    p = pat_fit_crf(y',opts);
    crf = pat_crf(p,opts.samprate,opts.maxt);
    
    % plot results
    plot(x(ix),y,'-k','LineWidth',4); hold on;
    plot(x(ix),crf,'-b','LineWidth',4)
    set(gca,'FontSize',25,'XLim',[min(x(ix)) max(x(ix))]);
    xlabel('Time (msec)','FontSize',25);
    ylabel('Normalized response','FontSize',25);
    legend({'Data','Model'},'FontSize',25);