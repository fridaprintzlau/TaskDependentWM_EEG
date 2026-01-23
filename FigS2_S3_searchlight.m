%% Decoding searchlight analysis
% Supplementary Figure S2 and S3 in 'Categorisation dynamically
% reconfigures working memory in a task-dependent manner' by Frida
% Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2024

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

run_analysis = false; %run analysis or load saved results

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);

%% Define channel neighbors

load('fukuda_layout.mat')
select_chans = find(ismember(lay.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
    'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'}));
label = lay.label(select_chans);
chanpos = lay.pos(select_chans,:);
nsensors = length(label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
    dist(i,:) = sqrt(sum((chanpos(1:nsensors,:) - repmat(chanpos(i,:), nsensors, 1)).^2,2))';
end

neighbours = struct;
nneighbs   = 6;
for i = 1:nsensors
    [~,isort] = sort(dist(i,:));
    neighbours(i).label       = label{i};
    neighbours(i).neighblabel = label(isort(2:(nneighbs+1)));
end

%%
if run_analysis
    for isub = 1:nsubs
        %Load data
        substrg  = sprintf('S%02d',sublist(isub));
        subpath  = sprintf('%s%s%s',datapath,fs,substrg);
        filename = sprintf('%s%s%s_ftEEGfinal.mat',subpath,fs,substrg);
        if ~exist(filename,'file')
            error(sprintf('%s does not exist!',filename));
        end
        load(filename);

        %Combine channels into single trialxchannelxtime matrix
        ntrials = new_behav.ntrial;
        dat = [];
        for itrl = 1:ntrials
            dat(itrl,:,:) = clean_data.trial{itrl};
        end
        time_orig = clean_data.time{1};

        %Define shape bins
        bin_size = 30; % 30 deg for 12 bins between 1-360 deg.
        bin_angles_deg = 0:bin_size:360;
        to_add = 0;
        nsets = 2; %2 different bin sets
        nangs = length(bin_angles_deg) - 1;
        bin_angles = NaN(nsets,nangs+1);
        for iset = 1:nsets
            bin_angles(iset,:) = bin_angles_deg + to_add;
            to_add = to_add + (bin_size/nsets);
        end
        bin_angles = [zeros(nsets,1), bin_angles];
        bin_angles_rad = bin_angles.*pi/180;

        %Get target labels & assign to bins
        target_ang = new_behav.target.*pi/180; %convert to radians
        bins = NaN(nsets, ntrials);
        for iset = 1:nsets
            [N, edges, bins(iset,:)] = histcounts(target_ang, bin_angles_rad(iset,:));
        end
        bins(bins == nangs+1) = 1;

        labels = bin_angles_deg.*pi/180;
        binned_shapes = NaN(nsets, ntrials);
        for i_bin = 1:length(labels)
            binned_shapes(bins == i_bin) = labels(i_bin);
        end

        %folding variable (blocks)
        foldvar = new_behav.blockCount;

        %Baseline correction
        onset = 0;
        baseline_window = time_orig>=-0.2 & time_orig<=-0.1;
        dat = bsxfun(@minus,dat,mean(dat(:,:,baseline_window),3));

        %Which channels to use?
        eeg_chans = find(ismember(clean_data.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
            'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'}));

        %% Decoder
        params.nreps = 1; %100; use 1 to debug; how many times to run the decoder to average over results? Use 100
        params.nfold = 8;
        params.chans = eeg_chans; %posterior chans

        fprintf('\nDecoding...')
        for irep = 1:params.nreps
            fprintf('%02d,',irep);
            for iset = 1:nsets
                angs = binned_shapes(iset,:); %angles you are trying to decode
                itrl = new_behav.acc==1;

                for currchan = 1:length(params.chans) %for each channel
                    n       = [{neighbours(currchan).label}; neighbours(currchan).neighblabel];
                    ichan   = find(ismember(clean_data.label,n));
                    m = sprintf('\nRep %d/%d, Channel %d/%d ',irep,params.nreps,currchan,27);
                    fprintf(m)

                    for itask=1:length(unique(new_behav.task_type))
                        %within task
                        itrl_train = new_behav.task_type==itask & new_behav.acc==1;
                        [dist_cos_within, dist_within] = mahal_theta_block_b3(dat(itrl_train,ichan,:), angs(itrl_train), foldvar(itrl_train));

                        %across task
                        itrl_test = new_behav.task_type==abs(itask-3) & new_behav.acc==1;
                        [dist_cos_across, dist_across] = mahal_theta_b_ind(dat(itrl_train,ichan,:), angs(itrl_train), dat(itrl_test,ichan,:), angs(itrl_test));

                        Cproj_temp_within(:,currchan,itask,iset,irep) = squeeze(mean(dist_cos_within, 2)); %average over trials
                        Cproj_temp_across(:,currchan,itask,iset,irep) = squeeze(mean(dist_cos_across, 2)); %average over trials
                    end

                end

            end %for iset = 1:nsets
        end %for irep = 1:params.nreps

        %Average decoding accuracy over repetitions of the decoder
        Cproj_temp_within = mean(Cproj_temp_within, 5);
        Cproj_temp_across = mean(Cproj_temp_across, 5);

        %Average decoding accuracy over bin-sets and save for each subject
        Cproj_within(isub,:,:,:)  = mean(Cproj_temp_within, 4);
        Cproj_across(isub,:,:,:)  = mean(Cproj_temp_across, 4);

        clear Cproj_temp_within Cproj_temp_across

    end
    %Save data
    savfile = sprintf('%s%sAnalysis_output/Shape_decoding_target_searchlight', datapath,fs);
    save(savfile, 'time_orig','Cproj_within', 'Cproj_across');
else
    %Load data
    loadfile = sprintf('%s%sAnalysis_output/Shape_decoding_target_searchlight.mat', datapath,fs);
    load(loadfile)
end

%% Plotting
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('qual', 'Set1', 12);
cmap2 = cbrewer('qual', 'Accent', 12);
cmap3 = cbrewer('div', 'PRGn', 12);
cmap4 = cbrewer('qual', 'Set2', 8);
cmap5 = cbrewer('div', 'BrBG', 12);

%load fieldtrip template file for accurate structure for plotting (we will
%overwrite the actual data)
substrg  = sprintf('S%02d',02);
subpath  = sprintf('%s%s%s',datapath,fs,substrg);
filename = sprintf('%s%s_ftEEGfinal_filtfirst_%s.mat',subpath,fs,substrg);
if ~exist(filename,'file')
    error(sprintf('%s does not exist!',filename));
end
load(filename, 'clean_data');

%% Within task decoding (Fig S2)

submean_m2s = squeeze(mean(Cproj_within(:,:,:,1),1));
submean_m2c = squeeze(mean(Cproj_within(:,:,:,2),1));

topodat = clean_data;
topodat.time = time_orig;
topodat.label = label;
topodat.dimord = 'time_chan';

%Match-to-sample
topodat.trial = submean_m2s;

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,1),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg = [];
cfg.xlim = [0 0.2];
cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlight ='on';
cfg.highlightchannel = label(p<=0.05);
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1,1,1];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2S Stimulus (0-0.2s)')
savname = strcat(figpath, 'm2s_searchlight_stim.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2s_searchlight_stim.svg');
saveas(gcf, savname)

%Early delay
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,1),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.2 0.7];
cfg.zlim = [0 0.0055];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2S Early delay (0.2-0.7s)')
savname = strcat(figpath, 'm2s_searchlight_early.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2s_searchlight_early.svg');
saveas(gcf, savname)

%Late delay
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,1),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.7 1.7];
cfg.zlim = [0 0.0014];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2S Late delay (0.7-1.7s)')
savname = strcat(figpath, 'm2s_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2s_searchlight_late.svg');
saveas(gcf, savname)


%Match-to-category
topodat.trial = submean_m2c;

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,2),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg = [];
cfg.xlim = [0 0.2];
cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlight ='on';
cfg.highlightchannel = label(p<=0.05/27);
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1,1,1];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2C Stimulus (0-0.2s)')
savname = strcat(figpath, 'm2c_searchlight_stim.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2c_searchlight_stim.svg');
saveas(gcf, savname)

%Early delay
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,2),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.2 0.7];
cfg.zlim = [0 0.0055];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2C Early delay (0.2-0.7s)')
savname = strcat(figpath, 'm2c_searchlight_early.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2c_searchlight_early.svg');
saveas(gcf, savname)

%Late delay
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(Cproj_within(:,toi,ichan,2),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05/27);
cfg.xlim = [0.7 1.7];
cfg.zlim = [0 0.0014];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('M2C Late delay (0.7-1.7s)')
savname = strcat(figpath, 'm2c_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2c_searchlight_late.svg');
saveas(gcf, savname)

%% Task difference topography
diff = squeeze(bsxfun(@minus,Cproj_within(:,:,:,2), Cproj_within(:,:,:,1)));
submean_diff = squeeze(mean(diff));

topodat = clean_data;
topodat.time = time_orig;
topodat.trial = submean_diff;
topodat.label = label;
topodat.dimord = 'time_chan';

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(diff(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg = [];
cfg.xlim = [0 0.2];
% cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlightchannel = label(p<=0.05);
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Task difference Stimulus (0-0.2s)')
savname = strcat(figpath, 'diff_searchlight_stim.png');
saveas(gcf, savname)
savname = strcat(figpath, 'diff_searchlight_stim.svg');
saveas(gcf, savname)


%Early delay
cfg = [];
cfg.xlim = [0.2 0.7];
% cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(diff(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg.highlightchannel = label(p<=0.05);
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Task difference Early delay (0.2-0.7s)')
savname = strcat(figpath, 'diff_searchlight_early.png');
saveas(gcf, savname)
savname = strcat(figpath, 'diff_searchlight_early.svg');
saveas(gcf, savname)

%Late delay
cfg = [];
cfg.xlim = [0.7 1.7];
% cfg.zlim = [0 0.0012];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(diff(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg.highlightchannel = label(p<=0.05);
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Task difference Late delay (1.2-1.7s)')
savname = strcat(figpath, 'diff_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'diff_searchlight_late.svg');
saveas(gcf, savname)

%% Within vs across task decoding (Fig S3B)

mean_within = squeeze(mean(Cproj_within,4));
mean_across = squeeze(mean(Cproj_across,4));
within_minus_across = mean_within - mean_across;

submean_within = squeeze(mean(mean_within,1));
submean_across = squeeze(mean(mean_across,1));
submean_diff = squeeze(mean(within_minus_across,1));

topodat = clean_data;
topodat.time = time_orig;
topodat.label = label;
topodat.dimord = 'time_chan';

%Within
topodat.trial = submean_within;

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(mean_within(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg = [];
cfg.xlim = [0 0.2];
cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlight ='on';
cfg.highlightchannel = label(p<=0.05);
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1,1,1];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Within Stimulus (0-0.2s)')

%Early delay
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(mean_within(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.2 0.7];
cfg.zlim = [0 0.0055];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Within Early delay (0.2-0.7s)')

%Late delay
p=[];
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(mean_within(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05/27);
cfg.xlim = [1.2 1.7];
cfg.zlim = [-0.0001 0.0012];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Within Late delay (1.2-1.7s)')
savname = strcat(figpath, 'within_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'within_searchlight_late.svg');
saveas(gcf, savname)


%Across task decoding
topodat.trial = submean_across;

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(mean_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg = [];
cfg.xlim = [0 0.2];
cfg.zlim = [0 0.0055];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlight ='on';
cfg.highlightchannel = label(p<=0.05);
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1,1,1];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Across Stimulus (0-0.2s)')

%Early delay
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(mean_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.2 0.7];
cfg.zlim = [0 0.0055];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Across Early delay (0.2-0.7s)')

%Late delay
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(mean_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05/27);
cfg.xlim = [0.7 1.7];
cfg.zlim = [-0.0001 0.0012];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Across Late delay (1.2-1.7s)')
savname = strcat(figpath, 'across_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'across_searchlight_late.svg');
saveas(gcf, savname)


%Within minus across task decoding
topodat.trial = submean_diff;

%Stimulus
for ichan = 1:length(label)
    toi = time_orig>=0 & time_orig<=0.2;
    time_mean = squeeze(mean(within_minus_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg = [];
cfg.xlim = [0 0.2];
cfg.zlim = [0 0.001];
cfg.layout = 'fukuda_layout.mat';
cfg.colorbar = 'southoutside';
cfg.highlight ='on';
cfg.highlightchannel = label(p<=0.05);
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1,1,1];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Diff Stimulus (0-0.2s)')

%Early delay
for ichan = 1:length(label)
    toi = time_orig>=0.2 & time_orig<=0.7;
    time_mean = squeeze(mean(within_minus_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end

cfg.highlightchannel = label(p<=0.05);
cfg.xlim = [0.2 0.7];
cfg.zlim = [0 0.001];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Diff Early delay (0.2-0.7s)')

%Late delay
for ichan = 1:length(label)
    toi = time_orig>=0.7 & time_orig<=1.7;
    time_mean = squeeze(mean(within_minus_across(:,toi,ichan),2));
    [h, p(ichan)] = ttest(time_mean,0);
end
cfg.highlightchannel = label(p<=0.05/27);
cfg.xlim = [0.7 1.7];
cfg.zlim = [0 0.001];
figure()
ft_topoplotER(cfg, topodat); colormap viridis; 
title('Diff Late delay (1.2-1.7s)')
savname = strcat(figpath, 'withinAcrossDiff_searchlight_late.png');
saveas(gcf, savname)
savname = strcat(figpath, 'withinAcrossDiff_searchlight_late.svg');
saveas(gcf, savname)
