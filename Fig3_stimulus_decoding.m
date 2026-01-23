%% Decoding analysis
% Figure 3 in 'Categorisation dynamically
% reconfigures working memory in a task-dependent manner' by Frida
% Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2024

%Dependencies:
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- Matlab 'Signal Processing Toolbox', 'Image Processing Toolbox', 
%   'Statistics and Machine Learning toolbox' & 'Parallel Computing
%   toolbox'

clear all
close all
clc

%Datapath
fs = filesep;
studypath = ''; %define path where the data folder sits
datapath = [studypath strcat(fs,'Data/EEG/Preprocessed')];  % input data
addpath(genpath('catWM_toolbox')) %required functions
figpath = ''; %define where to save the figures 

run_analysis = true; %run analysis or load saved results

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);
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

        % Folding variable for within-task decoding analysis
        foldvar = new_behav.blockCount;

        %Baseline correction
        onset = 0;
        baseline_window = time_orig>=-0.2 & time_orig<=-0.1;
        dat = bsxfun(@minus,dat,mean(dat(:,:,baseline_window),3));

        %Which channels to use?
        eeg_chans = find(ismember(clean_data.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
            'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'}));

        %% Decoder
        params.nreps = 50; %100; use 1 to debug; how many times to run the decoder to average over results? Use 100
        params.nfold = 8;
        params.chans = eeg_chans; %posterior chans

        %Use this for saving trialwise decoding results. Cannot use in combination with parallel pools...
        % trialwise_dist_temp = nan(ntrials, length(time_orig), nsets, params.nreps);
        % trialwise_Cmap_temp = nan(length(unique(binned_shapes)), ntrials, length(time_orig), nsets, params.nreps);

        fprintf('\nDecoding...')
        parfor irep = 1:params.nreps
            fprintf('%02d,',irep);
            for iset = 1:nsets
                angs = binned_shapes(iset,:); %angles you are trying to decode
                itrl = new_behav.acc==1; %Correct trials
                % itrl = 1:size(dat,1); %All trials
                [dist_cos, dist] = mahal_theta_kfold_b3(dat(itrl,params.chans,:), angs(itrl), params.nfold);
                Cproj_temp(:,iset,irep) = squeeze(mean(dist_cos, 2)); %average over trials
                Cmap_temp(:,:,iset,irep) = squeeze(mean(dist, 2)); %average over trials

                % %Save trialwise results
                % trialwise_dist_temp(itrl,:,iset,irep) = squeeze(dist_cos);
                % trialwise_Cmap_temp(:,itrl,:,iset,irep) = squeeze(dist);

                %Decoding within and across tasks
                for itask=1:2
                    %within task
                    itrl_train = new_behav.task_type==itask & new_behav.acc==1; %Correct trials
                    % itrl_train = new_behav.task_type==itask; %All trials
                    [dist_cos_within, dist_within] = mahal_theta_block_b3(dat(itrl_train,params.chans,:), angs(itrl_train), foldvar(itrl_train));

                    %across task
                    itrl_test = new_behav.task_type==abs(itask-3) & new_behav.acc==1; %Correct trials
                    % itrl_test = new_behav.task_type==abs(itask-3); %All trials
                    [dist_cos_across, dist_across] = mahal_theta_b_ind(dat(itrl_train,params.chans,:), angs(itrl_train), dat(itrl_test,params.chans,:), angs(itrl_test));

                    Cproj_temp_within(:,itask,iset,irep) = squeeze(mean(dist_cos_within, 2)); %average over trials
                    Cmap_temp_within(:,:,itask,iset,irep) = squeeze(mean(dist_within, 2)); %average over trials
                    Cproj_temp_across(:,itask,iset,irep) = squeeze(mean(dist_cos_across, 2)); %average over trials
                    Cmap_temp_across(:,:,itask,iset,irep) = squeeze(mean(dist_across, 2)); %average over trials

                end %for itask=1:2
            end %for iset = 1:nsets
        end %for irep = 1:params.nreps

        %Average decoding accuracy over repetitions of the decoder
        Cproj_temp = mean(Cproj_temp, 3);
        Cmap_temp = mean(Cmap_temp, 4);
        Cproj_temp_within = mean(Cproj_temp_within, 4);
        Cmap_temp_within = mean(Cmap_temp_within, 5);
        Cproj_temp_across = mean(Cproj_temp_across, 4);
        Cmap_temp_across = mean(Cmap_temp_across, 5);
        % trialwise_dist_temp = mean(trialwise_dist_temp, 4);
        % trialwise_Cmap_temp = mean(trialwise_Cmap_temp, 5);

        %Average decoding accuracy over bin-sets and save for each subject
        Cproj(isub,:,:)  = mean(Cproj_temp, 2);
        Cmap(isub,:,:,:) = mean(Cmap_temp, 3);
        Cproj_within(isub,:,:,:)  = mean(Cproj_temp_within, 3);
        Cmap_within(isub,:,:,:,:) = mean(Cmap_temp_within, 4);
        Cproj_across(isub,:,:,:)  = mean(Cproj_temp_across, 3);
        Cmap_across(isub,:,:,:,:) = mean(Cmap_temp_across, 4);
        % trialwise_dist{isub} = mean(trialwise_dist_temp, 3);
        % trialwise_Cmap{isub} = mean(trialwise_Cmap_temp, 4);

        clear Cproj_temp Cmap_temp Cproj_temp_within Cmap_temp_within Cproj_temp_across Cmap_temp_across trialwise_dist_temp trialwise_Cmap_temp

    end

    %% Save data

    savfile = sprintf('%s%sAnalysis_output/Shape_decoding_target_allchans_corrtrials', datapath,fs);
    % trlws_file = sprintf('%sShape_decoding_target_trlws_allchans_alltrials', datapath);

    save(savfile, 'time_orig','Cproj','Cproj_within','Cproj_across', 'Cmap');
    % save(trlws_file, 'trialwise_dist', 'trialwise_Cmap', '-v7.3');
else
    %Load data
    loadfile = sprintf('%s%sAnalysis_output/Shape_decoding_target_allchans_corrtrials.mat', datapath,fs); %All channels, correct trials
    load(loadfile)
end

%% Prepare for plotting
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('qual', 'Set1', 12);
cmap2 = cbrewer('qual', 'Accent', 12);
cmap3 = cbrewer('div', 'PRGn', 12);
cmap4 = cbrewer('qual', 'Set2', 8);
cmap5 = cbrewer('div', 'BrBG', 12);

% Smooth data
ftype   = 'gaussian'; % smooth the data with a gaussian filter
fsample = 250;
fsize   = (0.04*fsample); % 0.02*500 = 10 ms
smooth_Cproj = filtfast(Cproj,2,[],ftype,fsize);
smooth_Cmap  = filtfast(Cmap,3,[],ftype,fsize);
smooth_Cproj_within = filtfast(Cproj_within,2,[],ftype,fsize);
smooth_Cproj_across = filtfast(Cproj_across,2,[],ftype,fsize);

%% Plot pattern similarity over time (Figure 3A)
figure('Color', 'white', 'Position', [100 100 500 250])
imagesc(time_orig, [], squeeze(mean(smooth_Cmap,1))); cb=colorbar; colormap bone
xlim([-0.2 2.2])
yticks(1:18)
yticklabels(-180:30:179)
ylims=ylim;
line([0, 0], [ylims(1), ylims(2)], 'color', [1,1,1])
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', [1,1,1])
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', [1,1,1])
xlabel('Time from stimulus onset (s)', 'FontSize',12)
ylabel('Angular distance from target shape (deg)', 'FontSize',12)
ylabel(cb,'Pattern similarity', 'FontSize',12,'Rotation',270)
savname = strcat(figpath, 'tuning_function_overall.png');
saveas(gcf, savname)
savname = strcat(figpath, 'tuning_function_overall.svg');
saveas(gcf, savname)

%% Plot decoding over time (Figure 3B)

%Run cluster-corrected stats over time-resolved decoding
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1]; %cued and uncued
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = smooth_Cproj(plotsubs,:,plotcond);
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = cmap4(3,:);
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
hplot = plotpatch(smooth_Cproj,time_orig,cmap4(3,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
xlim([-0.2 2.2])
ylim(ylims);
ylabel('Decoding quality','FontSize',12)
xlabel('time from stimulus onset (s)','FontSize',12)
savname = strcat(figpath, 'decoding_quality_overall.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decoding_quality_overall.svg');
saveas(gcf, savname)

% Calculate time-point where decoding plateaus
dy=diff(mean(smooth_Cproj))./diff(time_orig);
figure()
plot(time_orig(2:end),dy); hold on
plot(time_orig(2:end), 0.01*ones(1,length(time_orig(2:end))), 'k--')
plot(time_orig(2:end), -0.01*ones(1,length(time_orig(2:end))), 'k--')

%% Plot within-task decoding (Figure 3CD)

%Run cluster-corrected stats over time-resolved decoding
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2]; %DMS and DMC
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = smooth_Cproj_within(plotsubs,:,plotcond);
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = cmap2(1:2,:);
ylims = [-0.002 0.012];
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

cfg=[];
hplot = plotpatch(smooth_Cproj_within,time_orig,cmap2,cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
ylims = ylim;
xlim([-0.2 2.2])
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
legend(hplot,{'m2sample', 'm2category'})
title('Within task decoding','FontSize',12)
ylabel('Decoding quality','FontSize',12)
xlabel('Time from stimulus onset (s)','FontSize',12)
savname = strcat(figpath, 'decoding_quality_withinTask.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decoding_quality_withinTask.svg');
saveas(gcf, savname)

%Averging over time window of interest
toi = time_orig>=0 & time_orig <=0.7;
within_early = squeeze(mean(smooth_Cproj_within(:,toi,:),2));

N = nsubs;
yMean_early = squeeze(mean(within_early,1));
ySEM_early = squeeze(std(within_early)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_early = bsxfun(@times, ySEM_early, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(within_early(:,1))), within_early(:,1), 30, cmap2(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(within_early(:,2))), within_early(:,2), 30, cmap2(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_early(1,1), yCI95_early(2,1), '.', 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_early(1,2), yCI95_early(2,2), '.', 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [within_early(:,1), within_early(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
title('decoding quality (early 0-0.7 s)')
savname = strcat(figpath, 'decodingquality_task_earlyTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decodingquality_task_earlyTOI.svg');
saveas(gcf, savname)

%Stats
[h p ci stat] = ttest(within_early(:,1),0,'tail','right')
es = meanEffectSize(within_early(:,1),Effect="cohen")

[h p ci stat] = ttest(within_early(:,2),0,'tail','right')
es = meanEffectSize(within_early(:,2),Effect="cohen")

[h p ci stat] = ttest(within_early(:,1),within_early(:,2))

toi = time_orig>=0.7 & time_orig <=1.7;
within_late = squeeze(mean(smooth_Cproj_within(:,toi,:),2));

N = nsubs;
yMean_late = squeeze(mean(within_late,1));
ySEM_late = squeeze(std(within_late)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(within_late(:,1))), within_late(:,1), 30, cmap2(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(within_late(:,2))), within_late(:,2), 30, cmap2(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [within_late(:,1), within_late(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
title('decoding quality (late 0.7-1.7 s)')
savname = strcat(figpath, 'decodingquality_task_lateTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decodingquality_task_lateTOI.svg');
saveas(gcf, savname)

%Stats
[h p ci stat] = ttest(within_late(:,1),0,'tail','right')
es = meanEffectSize(within_late(:,1),Effect="cohen")

[h p ci stat] = ttest(within_late(:,2),0,'tail','right')
es = meanEffectSize(within_late(:,2),Effect="cohen")

[h p ci stat] = ttest(within_late(:,1),within_late(:,2))

%% Plot cross-task decoding (Figure 3FG)

%Cross-task decoding across both tasks
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2]; %cued and uncued
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = cat(3, mean(smooth_Cproj_within,3), mean(smooth_Cproj_across,3));
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = [cmap2(3:4,:);[0,0,0]];
% ylims = [-0.002 0.013];
ylims = [-0.002 0.01];
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

cfg=[];
hplot = plotpatch(mean(smooth_Cproj_within,3),time_orig,cmap2(3,:),cfg); hold on;
gplot = plotpatch(mean(smooth_Cproj_across,3),time_orig,cmap2(4,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
ylim(ylims)
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
legend([hplot,gplot],{'within', 'across'})
xlim([-0.2 2.2])
title('X-task decoding')
ylabel('decoding quality')
xlabel('time from stimulus onset (s)')
savname = strcat(figpath, 'decoding_quality_acrossTask_both.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decoding_quality_acrossTask_both.svg');
saveas(gcf, savname)


statdata = cat(3, mean(Cproj_within,3), mean(Cproj_across,3));


toi = time_orig>=0 & time_orig <=0.7;
within_early_both = squeeze(mean(statdata(:,toi,1),2));
across_early_both = squeeze(mean(statdata(:,toi,2),2));
within_across_early = [within_early_both,across_early_both];

N = nsubs;
yMean_late = squeeze(mean(within_across_early,1));
ySEM_late = squeeze(std(within_across_early)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(within_across_early(:,1))), within_across_early(:,1), 30, cmap2(3,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(within_across_early(:,2))), within_across_early(:,2), 30, cmap2(4,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap2(3,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap2(4,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [within_across_early(:,1), within_across_early(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'within', 'across'})
title('decoding quality (early 0-0.7 s)')

[h p ci stat] = ttest(within_across_early(:,1),0,'tail','right')
es = meanEffectSize(within_across_early(:,1),Effect="cohen")

[h p ci stat] = ttest(within_across_early(:,2),0,'tail','right')
es = meanEffectSize(within_across_early(:,2),Effect="cohen")

[h p ci stat] = ttest(within_across_early(:,1),within_across_early(:,2),'tail','right')


toi = time_orig>=0.7 & time_orig <=1.7;
within_late_both = squeeze(mean(statdata(:,toi,1),2));
across_late_both = squeeze(mean(statdata(:,toi,2),2));
within_across_late = [within_late_both,across_late_both];

N = nsubs;
yMean_late = squeeze(mean(within_across_late,1));
ySEM_late = squeeze(std(within_across_late)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(within_across_late(:,1))), within_across_late(:,1), 30, cmap2(3,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(within_across_late(:,2))), within_across_late(:,2), 30, cmap2(4,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap2(3,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap2(4,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [within_across_late(:,1), within_across_late(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'within', 'across'})
title('decoding quality (late 0.7-1.7 s)')
savname = strcat(figpath, 'decodingquality_xtask_lateTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decodingquality_xtask_lateTOI.svg');
saveas(gcf, savname)

[h p ci stat] = ttest(within_across_late(:,1),0,'tail','right')
es = meanEffectSize(within_across_late(:,1),Effect="cohen")

[h p ci stat] = ttest(within_across_late(:,2),0,'tail','right')
es = meanEffectSize(within_across_late(:,2),Effect="cohen")

[h p ci stat] = ttest(within_across_late(:,1),within_across_late(:,2),'tail','right')
es = meanEffectSize(within_across_late(:,1),within_across_late(:,2), 'Paired', true, Effect="cohen")

%% Plot cross-task decoding separately by task (Figure S3)

%DMS task
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2]; 
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = cat(3, smooth_Cproj_within(:,:,1), smooth_Cproj_across(:,:,2));
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = cmap2(3:4,:);
ylims = [-0.002 0.01];
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

cfg=[];
hplot = plotpatch(smooth_Cproj_within(:,:,1),time_orig,cmap2(3,:),cfg); hold on;
gplot = plotpatch(smooth_Cproj_across(:,:,2),time_orig,cmap2(4,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
% ylim([-4*(10^-3) 18*(10^-3)])
ylims = ylim;
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
legend([hplot,gplot],{'within', 'across'})
xlim([-0.2 2.2])
title('test m2sample')
ylabel('decoding quality')
xlabel('time from stimulus onset (s)')
savname = strcat(figpath, 'decoding_quality_acrossTask_m2s.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decoding_quality_acrossTask_m2s.svg');
saveas(gcf, savname)

toi = time_orig>=0.7 & time_orig <=1.7;
within_late_m2s = squeeze(mean(smooth_Cproj_within(:,toi,1),2));
across_late_m2s = squeeze(mean(smooth_Cproj_across(:,toi,2),2));
within_across_late_m2s = [within_late_m2s,across_late_m2s];
m2s_cost = bsxfun(@minus, within_late_m2s, across_late_m2s);

[h p ci stat] = ttest(within_across_late_m2s(:,1),0,'tail','right')
[h p ci stat] = ttest(within_across_late_m2s(:,2),0,'tail','right')
[h p ci stat] = ttest(within_across_late_m2s(:,1),within_across_late_m2s(:,2),'tail','right')
es = meanEffectSize(within_across_late_m2s(:,1),within_across_late_m2s(:,2), 'Paired', true, Effect="cohen")


%DMC task
plotsubs = 1:length(sublist);
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2];
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = cat(3, smooth_Cproj_within(:,:,2), smooth_Cproj_across(:,:,1));
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);


%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = [cmap2(3:4,:);[0,0,0]];
ylims = [-0.002 0.01];
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

hplot = plotpatch(smooth_Cproj_within(:,:,2),time_orig,cmap2(3,:),cfg); hold on;
gplot = plotpatch(smooth_Cproj_across(:,:,1),time_orig,cmap2(4,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
% ylim([-4*(10^-3) 18*(10^-3)])
ylims = ylim;
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
legend([hplot,gplot],{'within', 'across'})
xlim([-0.2 2.2])
title('test m2category')
ylabel('decoding quality')
xlabel('time from stimulus onset (s)')
savname = strcat(figpath, 'decoding_quality_acrossTask_m2c.png');
saveas(gcf, savname)
savname = strcat(figpath, 'decoding_quality_acrossTask_m2c.svg');
saveas(gcf, savname)


toi = time_orig>=0.7 & time_orig <=1.7;
within_late_m2c = squeeze(mean(smooth_Cproj_within(:,toi,2),2));
across_late_m2c = squeeze(mean(smooth_Cproj_across(:,toi,1),2));
within_across_late_m2c = [within_late_m2c,across_late_m2c];
m2c_cost = bsxfun(@minus, within_late_m2c, across_late_m2c);

[h p ci stat] = ttest(within_across_late_m2c(:,1),0,'tail','right')
[h p ci stat] = ttest(within_across_late_m2c(:,2),0,'tail','right')
[h p ci stat] = ttest(within_across_late_m2c(:,1),within_across_late_m2c(:,2),'tail','right')
es = meanEffectSize(within_across_late_m2c(:,1),within_across_late_m2c(:,2), 'Paired', true, Effect="cohen")

%Test for difference in cost between DMS and DMC tasks
[h p ci stat] = ttest(m2s_cost,m2c_cost)


%% Split cross-decoding by degree of category-dependent representational shift (Figure 6DE)

%Load category-shift file 
load(sprintf('%s%sAnalysis_output/category_shift_m2c_lateTOI.mat', datapath,fs))
median_warp = median(m2c_toi);
strong_warp_subs = m2c_toi>median_warp;
weak_warp_subs = m2c_toi<median_warp;

for i_warp=1:2
    %Cross-task decoding across both tasks
    switch i_warp
        case 1, plotsubs = strong_warp_subs;
        case 2, plotsubs = weak_warp_subs;
    end
    corrwin  = [0 2.2];
    nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
    plotcond      = [1,2]; %cued and uncued
    np            = length(plotcond);
    pthresh       = 0.1; %threshold for plotting significance bar

    plotdata = cat(3, mean(smooth_Cproj_within(plotsubs,:,:),3), mean(smooth_Cproj_across(plotsubs,:,:),3));
    timeseries = time_orig;

    [decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

    %Plot time-resolved decoding quality and cluster-corrected stats
    figure('Color', 'white', 'Position', [100 100 600 400])
    statcols = [cmap2(3:4,:);[0,0,0]];
    % ylims = [-0.002 0.013];
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

    cfg=[];
    hplot = plotpatch(mean(smooth_Cproj_within(plotsubs,:,:),3),time_orig,cmap2(3,:),cfg); hold on;
    gplot = plotpatch(mean(smooth_Cproj_across(plotsubs,:,:),3),time_orig,cmap2(4,:),cfg);
    plot(time_orig, zeros(size(time_orig)), 'k--')
    ylim(ylims)
    line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
    line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
    line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
    line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
    legend([hplot,gplot],{'within', 'across'})
    xlim([-0.2 2.2])
    title('X-task decoding')
    ylabel('decoding quality')
    xlabel('time from stimulus onset (s)')
    switch i
        case 1
            savname = strcat(figpath, 'xtask_decoding_strongWarpSubs.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'xtask_decoding_strongWarpSubs.svg');
            saveas(gcf, savname)
        case 2
            savname = strcat(figpath, 'xtask_decoding_weakWarpSubs.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'xtask_decoding_weakWarpSubs.svg');
            saveas(gcf, savname)
    end


    toi = time_orig>=0.7 & time_orig <=1.7;
    within_late_both = squeeze(mean(plotdata(:,toi,1),2));
    across_late_both = squeeze(mean(plotdata(:,toi,2),2));
    within_across_late = [within_late_both,across_late_both];

    N = nsubs;
    yMean_late = squeeze(mean(within_across_late,1));
    ySEM_late = squeeze(std(within_across_late)./sqrt(N));
    CI95 = tinv([0.025 0.975], N-1);
    yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

    figure('Color', 'white', 'Position', [100 100 300 300])
    ax1=scatter(1.1*ones(1,length(within_across_late(:,1))), within_across_late(:,1), 30, cmap2(3,:), 'filled'); hold on
    ax2=scatter(1.9*ones(1,length(within_across_late(:,2))), within_across_late(:,2), 30, cmap2(4,:), 'filled');
    ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
    errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap2(3,:), 'linewidth', 2, 'MarkerSize', 30)
    errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap2(4,:), 'linewidth', 2, 'MarkerSize', 30)
    plot([1.2,1.8], [within_across_late(:,1), within_across_late(:,2)], 'color', [0.7,0.7,0.7])
    plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
    xlim([0.5, 2.5])
    xticks([1,2])
    xticklabels({'within', 'across'})
    title('decoding quality (late 0.7-1.7 s)')
    switch i_warp
        case 1
            savname = strcat(figpath, 'xtask_decoding_lateTOI_strongWarpSubs.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'xtask_decoding_lateTOI_strongWarpSubs.svg');
            saveas(gcf, savname)
        case 2
            savname = strcat(figpath, 'xtask_decoding_lateTOI_weakWarpSubs.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'xtask_decoding_lateTOI_weakWarpSubs.svg');
            saveas(gcf, savname)
    end

    [h p ci stat] = ttest(within_across_late(:,1),0,'tail','right')
    [h p ci stat] = ttest(within_across_late(:,2),0,'tail','right')
    [h p ci stat] = ttest(within_across_late(:,1),within_across_late(:,2),'tail','right')
    es = meanEffectSize(within_across_late(:,1),within_across_late(:,2), 'Paired', true, Effect="cohen")

    xdecoding_cost(i_warp,:) = within_across_late(:,1) - within_across_late(:,2);
end

[h p ci stats] = ttest2(xdecoding_cost(1,:), xdecoding_cost(2,:), 'tail', 'right')
es = meanEffectSize(xdecoding_cost(1,:),xdecoding_cost(2,:), 'Paired', false, Effect="cohen")
