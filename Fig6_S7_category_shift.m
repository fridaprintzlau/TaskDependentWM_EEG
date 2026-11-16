%% Category-dependent representational shift
% Figure 6 in 'Categorisation dynamically
% reconfigures working memory in a task-dependent manner' by Frida
% Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2024

%Dependencies:
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- Matlab 'Signal Processing Toolbox', 'Image Processing Toolbox' &
%   'Statistics and Machine Learning toolbox'

clear all
close all
clc

%Datapath
fs = filesep;
studypath = ''; %define path where the data folder sits
datapath = [studypath strcat(fs,'Data/EEG/Preprocessed')];  % input data
addpath(genpath('catWM_toolbox')) %required functions
figpath = ''; %define where to save the figures

run_analysis = false; %run analysis or load saved results

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);

%% Subject loop
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

        %Rotate stimulus values according to category
        new_bound = 180; %align everything to boundary=180
        shift = new_bound - new_behav.bound(1);
        new_behav.target_rot = mod((new_behav.target+shift), 360); %rotate and wrap between 0 and 360

        %Get target labels & assign to bins
        target_ang = new_behav.target_rot.*pi/180; %convert to radians
        theta = target_ang;

        %Baseline correction
        onset = 0;
        baseline_window = time_orig>=-0.2 & time_orig<=-0.1;
        dat = bsxfun(@minus,dat,mean(dat(:,:,baseline_window),3));

        %Which channels to use?
        eeg_chans = find(ismember(clean_data.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
            'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'}));
        chans = eeg_chans;
        %% Compute RDM and fit models

        stim = target_ang;
        u_stim = unique(stim);
        task = new_behav.task_type;
        u_task = unique(task);

        condxchanxtime = nan(length(u_stim)*length(u_task), length(chans), length(time_orig));
        resids = nan(ntrials, length(chans), length(time_orig));
        for itask = 1:length(u_task)
            for istim = 1:length(u_stim)
                curr_trls = stim==u_stim(istim) & task==itask;
                curr_cond = ((itask-1)*length(u_stim))+istim;
                resids(curr_trls,:,:) = dat(curr_trls,chans,:)-mean(dat(curr_trls,chans,:),1);
                condxchanxtime(curr_cond,:,:) = mean(dat(curr_trls,chans,:),1); %take within condition average
            end
        end
        timexchanxcond = permute(condxchanxtime, [3 2 1]);

        for t = 1:length(time_orig)
            make2D = squeeze(timexchanxcond(t,:,:)); %chans by conds
            stimDistances = pdist(make2D', 'mahalanobis', covdiag(resids(:,:,t)));

            %compute unzscored version for illustrative RDMs
            stimDistances_nozscore = stimDistances;
            stimdismatrix_nozscore = squareform(stimDistances_nozscore);
            RDM_output_nozscore(isub,:,:,t) = stimdismatrix_nozscore;

            %zscore condition distances to make coding models comparable
            stimDistances = zscore(stimDistances);
            stimdismatrix = squareform(stimDistances);
            RDM_output(isub,:,:,t) = stimdismatrix(:,:,:);

        end
    end

    %% Save data
    savfile = sprintf('%s%sAnalysis_output/RDM_results_allstim', datapath,fs);
    save(savfile, 'RDM_output', 'RDM_output_nozscore', 'time_orig','theta');
else
    %Load data
    load(sprintf('%s%sAnalysis_output/RDM_results_allstim.mat', datapath,fs))
end
%% Prepare for plotting
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('qual', 'Accent', 12);

%% Category-dependent shift (Figure 6BC)

% Calculate stimulus distances relative to each target
u_theta = unique(theta);
u_theta=circ_dist(u_theta,0); % to make sure that orientations are centered around 0
theta_dist=circ_dist2(u_theta',u_theta)';
theta_dist_abs = abs(theta_dist);
catmodel = squareform([ones(18)*0 ones(18)*1;...
    ones(18)*1 ones(18)*0]);
cats = squareform(catmodel);

RDM_m2s = RDM_output_nozscore(:,1:length(u_theta),1:length(u_theta),:);%Match-to-sample
RDM_m2c = RDM_output_nozscore(:,length(u_theta)+1:length(u_theta)*2,length(u_theta)+1:length(u_theta)*2,:);%Match-to-category

RDM_m2s_same=[];
RDM_m2s_diff=[];
RDM_m2c_same=[];
RDM_m2c_diff=[];
diff_score_m2s=[];
diff_score_m2c=[];
dup=[];
for i=1:length(u_theta)
    index_same_cat=[];
    index_diff_cat=[];
    same_cat = cats(i,i)==cats(i,:);
    for j=1:length(u_theta)
        dup = round(theta_dist_abs(i,:),3)==round(theta_dist_abs(i,j),3); %identify repeated distances
        index_same = find(dup & same_cat); %find repeated values from the same category
        index_diff = find(dup & ~same_cat); %find repeated values from the other category
        if isempty(index_same) || isempty(index_diff) %if duplicates are not present both within and across category
            continue
        else
            index_same_cat(j) = index_same;
            index_diff_cat(j) = index_diff;
        end
    end
    index_same_cat = nonzeros(unique(index_same_cat));
    index_diff_cat = nonzeros(unique(index_diff_cat));
    
    RDM_m2s_same = RDM_m2s(:,i,index_same_cat,:);
    RDM_m2s_diff = RDM_m2s(:,i,index_diff_cat,:);
    diff_score_m2s(i,:,:,:) = mean(bsxfun(@minus, RDM_m2s_diff, RDM_m2s_same),3);
    
    RDM_m2c_same = RDM_m2c(:,i,index_same_cat,:);
    RDM_m2c_diff = RDM_m2c(:,i,index_diff_cat,:);
    diff_score_m2c(i,:,:,:) = mean(bsxfun(@minus, RDM_m2c_diff, RDM_m2c_same),3);
end

mean_diff_score_m2s = squeeze(mean(diff_score_m2s,1));
mean_diff_score_m2c = squeeze(mean(diff_score_m2c,1));

% Smooth data for each subject
ftype   = 'gaussian'; % smooth the data with a gaussian filter (for regression analysis, not ERP/TF)
fsample = 250;
fsize   = (0.04*fsample);
smooth_diff_m2s = filtfast(mean_diff_score_m2s,2,[],ftype,fsize);
smooth_diff_m2c = filtfast(mean_diff_score_m2c,2,[],ftype,fsize);

%Plot
plotsubs = 1:nsubs;
corrwin  = [0 2.2];
nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
plotcond      = [1,2]; %DMS and DMC
np            = length(plotcond);
pthresh       = 0.1; %threshold for plotting significance bar

plotdata = squeeze(cat(3,smooth_diff_m2s,smooth_diff_m2c));
timeseries = time_orig;

[decod_b, decod_stats] = computeStats(plotdata, corrwin, timeseries, nperm, pthresh);

%Plot time-resolved decoding quality and cluster-corrected stats
figure('Color', 'white', 'Position', [100 100 600 400])
statcols = [cmap1(1:2,:);[0,0,0]];
ylims = [-0.015 0.025];
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

hold on
cfg = [];
hplot = plotpatch(smooth_diff_m2s,time_orig,cmap1(1,:),cfg); hold on
gplot = plotpatch(smooth_diff_m2c,time_orig,cmap1(2,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
ylim(ylims);
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
xlim([-0.2 2.2])
legend([hplot,gplot], {'M2S', 'M2C'})
ylabel('Neural distances between-within category')
xlabel('Time from stimulus onset (s)')
savname = strcat(figpath, 'distances_within_between.png');
saveas(gcf, savname)
savname = strcat(figpath, 'distances_within_between.svg');
saveas(gcf, savname)

%TOI
early_toi = time_orig>0 & time_orig <=0.7;

m2s_toi = mean(mean_diff_score_m2s(:,early_toi),2);
m2c_toi = mean(mean_diff_score_m2c(:,early_toi),2);
early_dat = [m2s_toi,m2c_toi];

N = nsubs;
yMean_early = squeeze(mean(early_dat,1));
ySEM_early = squeeze(std(early_dat)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_early = bsxfun(@times, ySEM_early, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(early_dat(:,1))), early_dat(:,1), 30, cmap1(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(early_dat(:,2))), early_dat(:,2), 30, cmap1(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_early(1,1), yCI95_early(2,1), '.', 'color', cmap1(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_early(1,2), yCI95_early(2,2), '.', 'color', cmap1(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [early_dat(:,1), early_dat(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
title('decoding quality (early 0-1 s)')
savname = strcat(figpath, 'distances_within_between_earlyTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'distances_within_between_earlyTOI.svg');
saveas(gcf, savname)

late_toi = time_orig>0.7 & time_orig <=1.7;
m2s_toi = mean(mean_diff_score_m2s(:,late_toi),2);
m2c_toi = mean(mean_diff_score_m2c(:,late_toi),2);

late_dat = [m2s_toi,m2c_toi];

N = nsubs;
yMean_late = squeeze(mean(late_dat,1));
ySEM_late = squeeze(std(late_dat)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(late_dat(:,1))), late_dat(:,1), 30, cmap1(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(late_dat(:,2))), late_dat(:,2), 30, cmap1(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap1(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap1(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [late_dat(:,1), late_dat(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:2.5], zeros(size([0.5:2.5])), 'k--')
plot([0.5:2.5], 0.0089*ones(size([0.5:2.5])), 'b--')
xlim([0.5, 2.5])
ylim([-0.05 0.06])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
title('category warping (late 0.7-1.7 s)')
savname = strcat(figpath, 'distances_within_between_lateTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'distances_within_between_lateTOI.svg');
saveas(gcf, savname)

[h,p,stats,ci]=ttest(m2s_toi,0,'tail','right')
es = meanEffectSize(m2s_toi,Effect="cohen")

[h,p,stats,ci]=ttest(m2c_toi,0,'tail','right')
es = meanEffectSize(m2c_toi,Effect="cohen")

[h,p,stats,ci]=ttest(m2s_toi,m2c_toi,'tail','left')
es = meanEffectSize(m2s_toi,m2c_toi, 'Paired', true, Effect="cohen")

% fname = sprintf('%s%sAnalysis_output/category_shift_m2c_lateTOI', datapath,fs);
% save(fname, 'm2c_toi') %Save for median split analysis

%% Figure S7 representational shift as a function of absolute distance to category boundary

near = [1,2,3,16,17,18,19,20,21,34,35,36];
mid = [4,5,6,13,14,15,22,23,24,31,32,33];
far = [7,8,9,10,11,12,25,26,27,28,29,30];
near_diff_score_m2s = squeeze(mean(diff_score_m2s(near,:,:),1));
near_diff_score_m2c = squeeze(mean(diff_score_m2c(near,:,:),1));
mid_diff_score_m2s = squeeze(mean(diff_score_m2s(mid,:,:),1));
mid_diff_score_m2c = squeeze(mean(diff_score_m2c(mid,:,:),1));
far_diff_score_m2s = squeeze(mean(diff_score_m2s(far,:,:),1));
far_diff_score_m2c = squeeze(mean(diff_score_m2c(far,:,:),1));

%Plot match-to-category only for each distance
% Smooth data for each subject
smooth_near = filtfast(near_diff_score_m2c,2,[],ftype,fsize);
smooth_middle = filtfast(mid_diff_score_m2c,2,[],ftype,fsize);
smooth_far = filtfast(far_diff_score_m2c,2,[],ftype,fsize);

figure('Color', 'white', 'Position', [100 100 600 400])
cfg = [];
jplot = plotpatch(smooth_far,time_orig,cmap8(12,:),cfg); hold on
gplot = plotpatch(smooth_middle,time_orig,cmap8(9,:),cfg);
hplot = plotpatch(smooth_near,time_orig,cmap8(4,:),cfg);
plot(time_orig, zeros(size(time_orig)), 'k--')
ylims = [-0.015 0.035];
ylim(ylims);
line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
xlim([-0.2 2.2])
legend([hplot,gplot,jplot], {'Near', 'Middle', 'Far'})
ylabel('Neural distances between-within category DMC task','FontSize',12)
xlabel('Time from stimulus onset (s)','FontSize',12)

late_toi = time_orig>0.7 & time_orig <=1.7;
m2c_near = mean(near_diff_score_m2c(:,late_toi),2);
m2c_mid = mean(mid_diff_score_m2c(:,late_toi),2);
m2c_far = mean(far_diff_score_m2c(:,late_toi),2);

late_dat = [m2c_near,m2c_mid,m2c_far];

N = nsubs;
yMean_late = squeeze(mean(late_dat,1));
ySEM_late = squeeze(std(late_dat)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_late = bsxfun(@times, ySEM_late, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(late_dat(:,1))), late_dat(:,1), 30, cmap8(4,:), 'filled'); hold on
ax2=scatter(2.1*ones(1,length(late_dat(:,2))), late_dat(:,2), 30, cmap8(9,:), 'filled');
ax3=scatter(3.1*ones(1,length(late_dat(:,3))), late_dat(:,3), 30, cmap8(12,:), 'filled');
ax1.MarkerFaceAlpha=0.3; ax2.MarkerFaceAlpha=0.3; ax3.MarkerFaceAlpha=0.3;
errorbar(1, yMean_late(1,1), yCI95_late(2,1), '.', 'color', cmap8(4,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_late(1,2), yCI95_late(2,2), '.', 'color', cmap8(9,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(3, yMean_late(1,3), yCI95_late(2,3), '.', 'color', cmap8(12,:), 'linewidth', 2, 'MarkerSize', 30)
% plot([1.2,1.8], [late_dat(:,1), late_dat(:,2)], 'color', [0.7,0.7,0.7])
plot([0.5:3.5], zeros(size([0.5:3.5])), 'k--')
xlim([0.5, 3.5])
% ylim([-0.15 0.15])
xticks([1,2,3])
xticklabels({'near', 'middle', 'far'})
title('category shift (late 0.7-1.7 s)','FontSize',12)

[h, p, ci, stat] = ttest(m2c_near, 0, 'tail', 'right')
[h, p, ci, stat] = ttest(m2c_mid, 0, 'tail', 'right')
[h, p, ci, stat] = ttest(m2c_far, 0, 'tail', 'right')

%one-way anova
subjects = (1:size(late_dat, 1))';
tbl_shift = array2table(late_dat, 'VariableNames', {'near', 'mid', 'far'});
tbl_shift.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_shift, 'near-far ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);