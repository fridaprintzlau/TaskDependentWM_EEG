% Supplementary Figure S8 in 'Task goals dynamically reconfigure working memory representations'
% by Frida Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2026

%Dependencies:
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- 'fukuda_layout.mat' file contains the EEG layout for creating neighbor
%    structure
%-- Matlab 'Statistics and Machine Learning toolbox'
%-- Cbrewer toolbox

% Decoding script
clear all
close all
clc

%Datapath
fs = filesep;
outputpath = '';  % output data path
inputpath = '';  % input data path
catWM_toolbox_path = '';
addpath(genpath(catWM_toolbox_path)) %decoding functions
figpath = ''; %path for saving figures

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);

%Load trialwise decoding
loadfile = ''; %path to trialwise decoding file
load(loadfile)

%%
for isub = 1:nsubs
    %Load data
    substrg  = sprintf('S%02d',sublist(isub));
    subpath  = sprintf('%s%s%s_WM',inputpath,fs,substrg);
    filename = sprintf('%s%s%s_ftEEGfinal.mat',subpath,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    time_orig = clean_data.time{1};

    sub_decoding = trialwise_dist{isub};
    n_blocks = max(new_behav.blockCount);

    first_half(isub,:) = squeeze(mean(sub_decoding(new_behav.blockCount<=n_blocks/2 & new_behav.acc==1,:),1));
    second_half(isub,:) = squeeze(mean(sub_decoding(new_behav.blockCount>n_blocks/2 & new_behav.acc==1,:),1));
end

%% Prepare for plotting
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('div', 'PiYG', 12);

% Smooth data
ftype   = 'gaussian'; % smooth the data with a gaussian filter
fsample = 250;
fsize   = (0.04*fsample); % 0.04*500 = 20 ms
smooth_early = filtfast(first_half,2,[],ftype,fsize);
smooth_late = filtfast(second_half,2,[],ftype,fsize);

%% Plot decoding over time
%Run cluster-corrected stats over time-resolved decoding
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2]; %cued and uncued
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = cat(3,smooth_early, smooth_late);
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = [cmap1([1,end],:);[0,0,0]];
ylims = [-0.002 0.010];
hold on

for i = 1:size(decod_b, 2)
    yheight = min(ylims) + range(ylims)*0.05*i;
    pheight = yheight + range(ylims)*0.02;
    stepsize = median(diff(decod_stats.sigtime))/4;
    for ib = 1:decod_b{i}.NumObjects
        currid = decod_b{i}.PixelIdxList{ib};
        if nnz(currid)>0
            xvec = [decod_stats.sigtime(currid(1))-stepsize decod_stats.sigtime(currid(end))+stepsize];
            plot(xvec,ones(length(xvec),1)*yheight,'k-','linewidth',4,'color',statcols(i,:));
            plotpval(gca,xvec,[1 1]*pheight,nanmean(decod_stats.pvals(i).pcorr(1,currid)),statcols(i,:),true);
            hold on
            decod_stats.cluster_pval{i}(ib) = nanmean(decod_stats.pvals(i).pcorr(1,currid));
            decod_stats.cluster_time{i}(ib,:) = xvec;
        end
    end
end

cfg = [];
hplot = plotpatch(plotdata,time_orig,cmap1([1,end],:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
xlim([-0.2 2.2])
ylim(ylims);
legend(hplot,{'first half', 'second half'},'FontSize',14)
ylabel('Decoding quality','FontSize',14)
xlabel('time from stimulus onset (s)','FontSize',14)

earlyTOI = squeeze(mean(plotdata(:,time_orig>0.2 & time_orig<1, :),2));
lateTOI = squeeze(mean(plotdata(:,time_orig>1 & time_orig<1.7, :),2));

N = nsubs;
% yMean_diff = squeeze(mean(earlyTOI,1));
% ySEM_diff = squeeze(std(earlyTOI)./sqrt(N));
% CI95 = tinv([0.025 0.975], N-1);
% yCI95_diff = bsxfun(@times, ySEM_diff, CI95(:));
% 
% figure('Color', 'white', 'Position', [100 100 300 300])
% ax1=scatter(1.1*ones(1,length(earlyTOI(:,1))), earlyTOI(:,1), 30, cmap1(1,:), 'filled'); hold on
% ax2=scatter(1.9*ones(1,length(earlyTOI(:,2))), earlyTOI(:,2), 30, cmap1(end,:), 'filled');
% ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
% errorbar(1, yMean_diff(1,1), yCI95_diff(2,1), '.', 'color', cmap1(1,:), 'linewidth', 2, 'MarkerSize', 30)
% errorbar(2, yMean_diff(1,2), yCI95_diff(2,2), '.', 'color', cmap1(end,:), 'linewidth', 2, 'MarkerSize', 30)
% plot([1.2,1.8], [earlyTOI(:,1), earlyTOI(:,2)], 'color', [0.7,0.7,0.7])
% plot(0.5:2.5, zeros(size(0.5:2.5)), 'k--')
% xlim([0.5, 2.5])
% xticks([1,2])
% xticklabels({'first half', 'second half'})
% title('decoding quality (early TOI)')

[h p stats ci] = ttest(earlyTOI(:,1),0,'tail','right')
[h p stats ci] = ttest(earlyTOI(:,2),0,'tail','right')
[h p stats ci] = ttest(earlyTOI(:,1),earlyTOI(:,2),'tail','left')

yMean_diff = squeeze(mean(lateTOI,1));
ySEM_diff = squeeze(std(lateTOI)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_diff = bsxfun(@times, ySEM_diff, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(lateTOI(:,1))), lateTOI(:,1), 30, cmap1(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(lateTOI(:,2))), lateTOI(:,2), 30, cmap1(end,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_diff(1,1), yCI95_diff(2,1), '.', 'color', cmap1(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_diff(1,2), yCI95_diff(2,2), '.', 'color', cmap1(end,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [lateTOI(:,1), lateTOI(:,2)], 'color', [0.7,0.7,0.7])
plot(0.5:2.5, zeros(size(0.5:2.5)), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'first half', 'second half'})
title('decoding quality (late TOI)')

[h p stats ci] = ttest(lateTOI(:,1),0,'tail','right')
[h p stats ci] = ttest(lateTOI(:,2),0,'tail','right')
[h p stats ci] = ttest(lateTOI(:,1),lateTOI(:,2))
