%% Decoding quality for correct and incorrect trials
% Supplementary Figure S1 in 'Categorisation dynamically
% reconfigures working memory in a task-dependent manner' by Frida
% Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2024

%Dependencies:
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- Matlab 'Signal Processing Toolbox' & 'Statistics and Machine Learning toolbox'

clear all
close all
clc

%Datapath
fs = filesep;
studypath = ''; %define path where the data folder sits
datapath = [studypath strcat(fs,'Data/EEG/Preprocessed')];  % input data
addpath(genpath('catWM_toolbox')) %required functions

%Load trialwise data
loadfile = sprintf('%s%sAnalysis_output/Shape_decoding_target_trlws_allchans_alltrials.mat', datapath,fs);
load(loadfile)

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);

for isub = 1:nsubs
    %Select trialwise decoding for each subject
    trialwise_quality = trialwise_dist{isub};

    %Load behavioural data
    substrg  = sprintf('S%02d',sublist(isub));
    subpath  = sprintf('%s%s%s',datapath,fs,substrg);
    filename = sprintf('%s%s%s_ftEEGfinal.mat',subpath,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename, 'new_behav');

    %Correct trials
    decoding_qual_m2c(isub,:,:) = squeeze(mean(trialwise_quality(new_behav.acc==1 & new_behav.task_type==2,:,:),1));
    decoding_qual_m2s(isub,:,:) = squeeze(mean(trialwise_quality(new_behav.acc==1 & new_behav.task_type==1,:,:),1));

    %Incorrect trials
    decoding_qual_m2c_incorr(isub,:,:) = squeeze(mean(trialwise_quality(new_behav.acc==0 & new_behav.task_type==2,:,:),1));
    decoding_qual_m2s_incorr(isub,:,:) = squeeze(mean(trialwise_quality(new_behav.acc==0 & new_behav.task_type==1,:,:),1));

end


%%
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('qual', 'Set1', 12);
cmap2 = cbrewer('qual', 'Accent', 12);
cmap_paired = cbrewer('qual', 'Paired', 12);

time_axis=-0.499:.004:3.2;

% Smooth data for each subject
ftype   = 'gaussian'; % smooth the data with a gaussian filter (for regression analysis, not ERP/TF)
fsample = 250;
fsize   = (0.04*fsample); % 0.02*500 = 10 ms
smooth_m2c = filtfast(squeeze(decoding_qual_m2c),2,[],ftype,fsize);
smooth_m2s = filtfast(squeeze(decoding_qual_m2s),2,[],ftype,fsize);
smooth_m2c_incorr = filtfast(squeeze(decoding_qual_m2c_incorr),2,[],ftype,fsize);
smooth_m2s_incorr = filtfast(squeeze(decoding_qual_m2s_incorr),2,[],ftype,fsize);

%Match-to-sample
figure()
cfg = [];
hplot = plotpatch(squeeze(smooth_m2s), time_axis, cmap_paired(4,:), cfg); hold on;
gplot = plotpatch(squeeze(smooth_m2s_incorr), time_axis, cmap_paired(3,:), cfg);
plot(time_axis, zeros(1,length(time_axis)), 'k--')
legend([hplot,gplot],{'correct', 'incorrect'})
xlim([-0.2, 2.2])
title('Decoding quality')

m2s_stat_dat = cat(3, squeeze(decoding_qual_m2s), squeeze(decoding_qual_m2s_incorr));
toi = time_axis>=0 & time_axis <=0.7;
m2s_early = squeeze(mean(m2s_stat_dat(:,toi,:),2));

N = nsubs;
yMean_early = squeeze(mean(m2s_early,1));
ySEM_early = squeeze(std(m2s_early)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_early = bsxfun(@times, ySEM_early, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(m2s_early(:,1))), m2s_early(:,1), 30, cmap_paired(4,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(m2s_early(:,2))), m2s_early(:,2), 30, cmap_paired(3,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_early(1,1), yCI95_early(2,1), '.', 'color', cmap_paired(4,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_early(1,2), yCI95_early(2,2), '.', 'color', cmap_paired(3,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [m2s_early(:,1), m2s_early(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'Correct', 'Incorrect'})
title('M2S decoding quality (early 0-0.7 s)')

[h p ci stat] = ttest(m2s_early(:,1),m2s_early(:,2))

toi = time_axis>=0.7 & time_axis <=1.7;
m2s_late = squeeze(mean(m2s_stat_dat(:,toi,:),2));

N = nsubs;
yMean_late = squeeze(mean(m2s_late,1));
ySEM_late = squeeze(std(m2s_late)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(m2s_late(:,1))), m2s_late(:,1), 30, cmap_paired(4,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(m2s_late(:,2))), m2s_late(:,2), 30, cmap_paired(3,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap_paired(4,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap_paired(3,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [m2s_late(:,1), m2s_late(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'Correct', 'Incorrect'})
title('M2S decoding quality (late 0.7-1.7 s)')

[h p ci stat] = ttest(m2s_late(:,1),m2s_late(:,2))

%Match-to-category
figure()
cfg = [];
hplot = plotpatch(squeeze(smooth_m2c), time_axis, cmap_paired(2,:), cfg); hold on;
gplot = plotpatch(squeeze(smooth_m2c_incorr), time_axis, cmap_paired(1,:), cfg);
plot(time_axis, zeros(1,length(time_axis)), 'k--')
legend([hplot,gplot],{'correct', 'incorrect'})
xlim([-0.2, 2.2])
title('Decoding quality')

m2c_stat_dat = cat(3, squeeze(decoding_qual_m2c), squeeze(decoding_qual_m2c_incorr));
toi = time_axis>=0 & time_axis <=0.7;
m2c_early = squeeze(mean(m2c_stat_dat(:,toi,:),2));

N = nsubs;
yMean_early = squeeze(mean(m2c_early,1));
ySEM_early = squeeze(std(m2c_early)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_early = bsxfun(@times, ySEM_early, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(m2c_early(:,1))), m2c_early(:,1), 30, cmap_paired(2,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(m2c_early(:,2))), m2c_early(:,2), 30, cmap_paired(1,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_early(1,1), yCI95_early(2,1), '.', 'color', cmap_paired(2,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_early(1,2), yCI95_early(2,2), '.', 'color', cmap_paired(1,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [m2c_early(:,1), m2c_early(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'Correct', 'Incorrect'})
title('M2C decoding quality (early 0-0.7 s)')

[h p ci stat] = ttest(m2c_early(:,1),m2c_early(:,2))

toi = time_axis>=0.7 & time_axis <=1.7;
m2c_late = squeeze(mean(m2c_stat_dat(:,toi,:),2));

N = nsubs;
yMean_late = squeeze(mean(m2c_late,1));
ySEM_late = squeeze(std(m2c_late)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(m2c_late(:,1))), m2c_late(:,1), 30, cmap_paired(2,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(m2c_late(:,2))), m2c_late(:,2), 30, cmap_paired(1,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap_paired(2,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap_paired(1,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [m2c_late(:,1), m2c_late(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'Correct', 'Incorrect'})
title('M2C decoding quality (late 0.7-1.7 s)')

[h p ci stat] = ttest(m2c_late(:,1),m2c_late(:,2))
