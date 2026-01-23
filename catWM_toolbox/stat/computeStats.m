function [b, stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh)

% select time window
plottimeind  = timeseries >= -Inf;
plottime     = timeseries(plottimeind);

%set up stats struct
it      = plottime>=corrwin(1) & plottime<= corrwin(2);
sigtime = plottime(it);
stats   = struct;
stats.pthresh = pthresh;
stats.corrwin = corrwin;
stats.nperm   = nperm;
stats.plottesttype = 'Ttest';
stats.sigtime      = sigtime;
stats.pvals = [];
% data for t-test
np = size(plotdata, 3);
if np == 1
    testdata(:, :, 1) = plotdata;
elseif np > 1 && np<3
    testdata(:, :, :) = plotdata;
    testdata(:,:,end+1) = diff(testdata,[],3); %only use for more than one cond
elseif np >=3
    testdata(:, :, :) = plotdata;
end
nsig = size(testdata,3);
for isig = 1:nsig
    [~,praw] = ttest(testdata(:,it,isig));
    if nperm
        pcorr   = ClusterCorrection2(testdata(:,it,isig),nperm,0.05);
    else
        pcorr = praw;
    end
    b{isig} = bwconncomp(pcorr(:,:)<pthresh); % finds the connected components that are above the threshold set by pthresh
    stats.pvals(isig).praw    = praw;
    stats.pvals(isig).pcorr   = pcorr;
end
clear testdata

end