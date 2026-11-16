% Supplementary Figure S5 in 'Task goals dynamically reconfigure working memory representations'
% by Frida Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2026

%Dependencies:
%-- Fieldtrip toolbox available for free download here: https://www.fieldtriptoolbox.org/
%   Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs
%   Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, 
%   EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, 2011; 2011:156869.
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- 'fukuda_layout.mat' file contains the EEG layout for creating neighbor
%    structure
%-- Matlab 'Statistics and Machine Learning toolbox'

clear all
close all
clc

%Datapath
fs = filesep;
studypath = ''; %define path where the data folder sits
datapath = [studypath strcat(fs,'Data/EEG/Preprocessed')];  % input data
addpath(genpath('catWM_toolbox')) %required functions
figpath = ''; %define where to save the figures 

fieldtrippath = ''; %add fieldtrip path here
addpath(fieldtrippath);
ft_defaults;

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);

%Settings
minf = 1; 
maxf = 30;
tons = -1;
toff = +3.000;

%%
for isub = 1:nsubs %subject loop
    %% load data
    substrg  = sprintf('S%02d',sublist(isub));
    subpath  = sprintf('%s%s%s_WM',datapath,fs,substrg);
    filename = sprintf('%s%s_ftEEGfinal_TF%s.mat',subpath,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);
    data = clean_data;
    behav = new_behav;
    ntrials = size(data.trial,2); % number of trials AFTER artefact rejection
    time    = data.time{1};

    usable = ntrials;
    fprintf('\n\nS%02d: %03d usable trials (%1.2f percent).\n',sublist(isub),usable,100*usable/864);
    
    overwrite_TF = 0;
    fname = sprintf('%s%sHanningTF_reref.mat',subpath,fs);
    if ~exist(fname,'file') || overwrite_TF

        chanlist = find(ismember(clean_data.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
            'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'})); %skip EOG channels
        nchans   = length(chanlist);

        %select eeg channels
        cfg.channel = chanlist;
        data_eeg = ft_selectdata(cfg, data);

        %re-reference to average of all channels
        cfg = [];
        cfg.reref      = 'yes';
        cfg.refmethod  = 'avg';
        cfg.refchannel = 'all';
        data_ref = ft_preprocessing(cfg, data_eeg);

        chanoi = chanlist;
        timeoi = data.time{1}(data.time{1} >= tons & data.time{1} <= toff); % time samples of interest (s)
        dstep  = 10; %used to be 5!
        timeoi = timeoi(1:dstep:end);

        freqoi = minf:maxf;

        cfg            = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.keeptrials = 'yes';
        %cfg.feedback   = 'none';
        cfg.channel    = chanoi;
        cfg.foi        = freqoi;
        cfg.toi        = timeoi;
        cfg.taper      = 'hanning';
        cfg.t_ftimwin  = 4./cfg.foi; % time window (s) - 5 cycles

        % datatfr = ft_freqanalysis(cfg,data);
        datatfr = ft_freqanalysis(cfg,data_ref);

        powspctrm_m2s = datatfr.powspctrm(behav.task_type==1,:,:,:);
        powspctrm_m2c = datatfr.powspctrm(behav.task_type==2,:,:,:);
        save(fname,'datatfr','cfg','powspctrm_m2s','powspctrm_m2c', '-v7.3')

    else
        fprintf('\nTF analysis already done! Loading existing file');
        load(fname);
    end

    cfg = [];
    cfg.baseline     = [-0.2 -0.1];      % pre-stimulus window
    cfg.baselinetype = 'absolute';             % log ratio (decibel)
    freq{isub} = ft_freqbaseline(cfg, datatfr);

end

%% Time-frequency plots
%Channel selection
posterior_chans = {'O1', 'O2', 'POz', 'PO3', 'PO4', 'PO7', 'PO8', 'Pz', 'P3', 'P4'};
frontal_chans = {'AF3', 'AFz', 'AF4', 'F3', 'Fz', 'F4', 'FC1', 'FC2'};

%% Average over frequencies

%Compute frontal theta
cfg = [];
cfg.channel     = frontal_chans;
cfg.frequency   = [4 7];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
cfg.avgoverrpt = 'yes';

theta_m2s_pow = [];
theta_m2c_pow = [];
for isub = 1:nsubs
    cfg.trials = freq{isub}.trialinfo == 12; %M2S
    theta_M2S{isub}  = ft_selectdata(cfg, freq{isub});

    cfg.trials = freq{isub}.trialinfo == 22; %M2C
    theta_M2C{isub} = ft_selectdata(cfg, freq{isub});

    theta_m2s_pow(isub,:) = squeeze(theta_M2S{isub}.powspctrm);
    theta_m2c_pow(isub,:) = squeeze(theta_M2C{isub}.powspctrm);
end

%Compute posterior alpha
cfg = [];
cfg.channel     = posterior_chans;
cfg.frequency   = [8 14];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
cfg.avgoverrpt = 'yes';

alpha_m2s_pow = [];
alpha_m2c_pow = [];
for isub = 1:nsubs
    cfg.trials = freq{isub}.trialinfo == 12; %M2S
    alpha_M2S{isub}  = ft_selectdata(cfg, freq{isub});

    cfg.trials = freq{isub}.trialinfo == 22; %M2C
    alpha_M2C{isub} = ft_selectdata(cfg, freq{isub});

    alpha_m2s_pow(isub,:) = squeeze(alpha_M2S{isub}.powspctrm);
    alpha_m2c_pow(isub,:) = squeeze(alpha_M2C{isub}.powspctrm);
end

%% Stats

design = zeros(2,2*nsubs);
design(1,:) = [1:nsubs 1:nsubs];
design(2,:) = [ones(1,nsubs) 2*ones(1,nsubs)];

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = 5000;
cfg.design           = design;
cfg.uvar             = 1;
cfg.ivar             = 2;
cfg.latency          = [0.2 1.7];

stat_theta = ft_freqstatistics(cfg, theta_M2S{:}, theta_M2C{:});
stat_alpha = ft_freqstatistics(cfg, alpha_M2S{:}, alpha_M2C{:});

%% Plot
colorpath = '/Volumes/FridaPrintz/IncidentalWM_EEG_Exp1/scripts/toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('qual', 'Set1', 12);
cmap2 = cbrewer('qual', 'Accent', 12);
cmap3 = cbrewer('div', 'PRGn', 12);
cmap4 = cbrewer('qual', 'Set2', 8);

%Plot frontal theta
time_ind = time_axis>=-0.2 & time_axis<=2.5;
figure()
cfg = [];
hplot = plotpatch(theta_m2s_pow(:,time_ind),time_axis(time_ind),cmap2(1,:),cfg);
jplot = plotpatch(theta_m2c_pow(:,time_ind),time_axis(time_ind),cmap2(2,:),cfg);
plot(time_axis(time_ind), zeros(1,length(time_axis(time_ind))), 'k--')
ylims = ylim;
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
ylim(ylims);
xlim([-0.2,2.2])
ylabel('frontal theta','FontSize',12)
xlabel('time from stimulus onset (s)','FontSize',12)

theta_diff = theta_m2c_pow - theta_m2s_pow;
theta_clusters = stat_theta.time(squeeze(stat_theta.negclusterslabelmat)>0 | squeeze(stat_theta.posclusterslabelmat)>0);

figure()
cfg = [];
hplot = plotpatch(theta_diff(:,time_ind),time_axis(time_ind),cmap4(3,:),cfg);
plot(time_axis(time_ind), zeros(1,length(time_axis(time_ind))), 'k--'); hold on
xlim([-0.2 2.2])
ylims = ylim;
yheight = min(ylims) + range(ylims)*0.1;
plot(theta_clusters,ones(length(theta_clusters),1)*yheight,'k.','linewidth',4,'color',cmap4(3,:)); %plot significance bar
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
ylim(ylims);
xlim([-0.2,2.2])
ylabel('Frontal theta','FontSize',12)
xlabel('time from stimulus onset (s)','FontSize',12)


%Theta late TOI
toi = time_axis>0.7 & time_axis<1.7;
theta_toi = squeeze(mean(theta_diff(:,toi),2));

N = nsubs;
yMean_theta = squeeze(mean(theta_toi,1));
ySEM_theta = squeeze(std(theta_toi)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_theta = bsxfun(@times, ySEM_theta, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(theta_toi)), theta_toi, 30, cmap4(3,:), 'filled'); hold on
ax1.MarkerFaceAlpha=0.5;
errorbar(1, yMean_theta(1,1), yCI95_theta(2,1), '.', 'color', cmap4(3,:), 'linewidth', 2, 'MarkerSize', 30)
plot([0.5:1.5], zeros(size([0.5:1.5])), 'k--')
xlim([0.5, 1.5])
title('theta difference (late 0.7-1.7 s)')

[h,p,ci,stats] = ttest(theta_toi,0)

%Plot posterior alpha
figure()
cfg = [];
hplot = plotpatch(alpha_m2s_pow(:,time_ind),time_axis(time_ind),cmap2(1,:),cfg);
jplot = plotpatch(alpha_m2c_pow(:,time_ind),time_axis(time_ind),cmap2(2,:),cfg);
plot(time_axis(time_ind), zeros(1,length(time_axis(time_ind))), 'k--')
xlim([-0.2 2.5])
ylims = ylim;
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
ylim(ylims);
xlim([-0.2,2.2])
ylabel('Posterior alpha','FontSize',12)
xlabel('time from stimulus onset (s)','FontSize',12)

alpha_diff = alpha_m2c_pow - alpha_m2s_pow;
% alpha_clusters = stat_alpha.time(squeeze(stat_alpha.negclusterslabelmat)>0 | squeeze(stat_alpha.posclusterslabelmat)>0);

figure()
cfg = [];
hplot = plotpatch(alpha_diff(:,time_ind),time_axis(time_ind),cmap4(3,:),cfg);
plot(time_axis(time_ind), zeros(1,length(time_axis(time_ind))), 'k--')
xlim([-0.2 2.5])
ylims = ylim;
yheight = min(ylims) + range(ylims)*0.1;
% plot(alpha_clusters,ones(length(alpha_clusters),1)*yheight,'k.','linewidth',4,'color',cmap4(3,:)); %plot significance bar
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
ylim(ylims);
xlim([-0.2,2.2])
ylabel('Posterior alpha','FontSize',12)
xlabel('time from stimulus onset (s)','FontSize',12)

%Alpha early TOI
toi = time_axis>0 & time_axis<0.7;
alpha_toi = squeeze(mean(alpha_diff(:,toi),2));

N = nsubs;
yMean_theta = squeeze(mean(alpha_toi,1));
ySEM_theta = squeeze(std(alpha_toi)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_theta = bsxfun(@times, ySEM_theta, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(alpha_toi)), alpha_toi, 30, cmap4(3,:), 'filled'); hold on
ax1.MarkerFaceAlpha=0.5;
errorbar(1, yMean_theta(1,1), yCI95_theta(2,1), '.', 'color', cmap4(3,:), 'linewidth', 2, 'MarkerSize', 30)
plot([0.5:1.5], zeros(size([0.5:1.5])), 'k--')
xlim([0.5, 1.5])
title('alpha difference (early 0-0.7 s)')

[h,p,ci,stats] = ttest(alpha_toi,0)
