%% Behavioural analysis
% Figure 1 & 2 in 'Categorisation dynamically
% reconfigures working memory in a task-dependent manner' by Frida
% Printzlau, Olya Bulatova, Michael Mack & Keisuke Fukuda

% - Frida Printzlau, October 2024

%Dependencies:
%-- 'catWM_toolbox' contains decoding, statistics and plotting functions
%   required for running the script
%-- Matlab 'Statistics and Machine Learning toolbox'


clear all
close all
clc

% Define paths
studypath = ''; %define path where the data folder sits
datapath = [studypath, '/Data/Behaviour/']; %Input data
figpath = ''; %add path for saving figures
colorpath = 'catWM_toolbox/cbrewer';
addpath(genpath(colorpath));

%Define subjects
sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs   = length(sublist);

%Initialise arrays
ntrials = 864;
WM_stim_allsub = NaN(nsubs,ntrials);
WM_stim_rot_allsub = NaN(nsubs,ntrials);
WM_acc_dist=[];
WM_rt_dist=[];

%% Loop through subjects
for isub = 1:nsubs
    %Load memory data
    substrg  = sprintf('s%02d',sublist(isub));
    filestr = fullfile(datapath,  strcat(substrg,'*__catWM_EEG_postLearn.mat'));
    fstruct = dir(filestr);
    if ~isempty(fstruct)
        filename = fstruct(1).name;
        WM = load(strcat(datapath,filename));
    else, error(sprintf('%s does not exist!',substrg));
    end
    
    %Load learning data
    filestr    = fullfile(datapath,  strcat(substrg,'*__catWM_learning.mat'));
    fstruct = dir(filestr);
    if ~isempty(filestr)
        filename = fstruct(1).name;
        learn = load(strcat(datapath,filename));
    else, error(sprintf('%s does not exist!',substrg));
    end

    %Which category condition
    counterbalance_cond = learn.dat.SubjDetails.Condition;
    switch counterbalance_cond
        case 1
            bound = [170 350];
            cat_mean = [80 260];
        case 2
            bound = [200 20];
            cat_mean = [290 110];
        case 3
            bound = [230 50];
            cat_mean = [320 140];
        case 4
            bound = [260 80];
            cat_mean = [350 170];
        case 5
            bound = [290 110];
            cat_mean = [20 200];
        case 6
            bound = [320 140];
            cat_mean = [50 230];
    end

    %% Demographics info
    age(isub) = str2num(WM.dat.SubjDetails.Age);
    switch WM.dat.SubjDetails.Sex
        case {'F','f','female'}
            sex(isub) = 1;
        case {'M','m','male'}
            sex(isub) = 2;
    end
    switch WM.dat.SubjDetails.Hand
        case {'R','r','right'}
            handedness(isub) = 1;
        case {'L','l','left'}
            handedness(isub) = 2;
        case {'both','L/R'}
            handedness(isub) = 3;
    end
    
    %% Learning
    learn_acc  = learn.dat.rsp.accuracy;
    learn_rt   = learn.dat.rsp.RT;
    learn_stim = learn.dat.stms.stim;
    block      = learn.dat.condition.blockCount;
    
    % Learning blocks
    unique_blocks = unique(block);
    for i=1:length(unique_blocks)
        block_acc(isub,i) = mean(learn_acc(block==i));
        block_rt(isub,i)  = nanmedian(learn_rt(block==i & learn_acc==1));
    end
    
    %Distance to boundary
    learn_dist_to_bound = abs(mod((learn_stim-bound(1))+90, 90*2)-90); % wrap between -90 and 90, so determining distance to either boundary
    unique_dist = unique(learn_dist_to_bound);
    for i=1:length(unique_dist)
        i_dist = unique_dist(i);
        dist_acc(isub,i) = mean(learn_acc(learn_dist_to_bound==i_dist));
        dist_rt(isub,i)  = median(learn_rt(learn_dist_to_bound==i_dist & learn_acc==1));
    end
    
    %% Working memory
    
    WM_trialcounter = WM.dat.condition.trialCount;
    WM_rt           = WM.dat.rsp.RT;
    WM_acc          = WM.dat.rsp.accuracy;
    WM_stim         = WM.dat.stms.target;
    WM_task         = WM.dat.condition.block_type; %1=match-to-sample, 2=match-to-category
    WM_match        = WM.dat.condition.trial_type; %1=match, 2=non-match
    WM_probe        = WM.dat.stms.probe;

    %Performance as a function of task type
    for i = 1:2
        task_acc(i) = mean(WM_acc(WM_task==i));
        task_rt(i) = nanmedian(WM_rt(WM_task==i & WM_acc==1));
    end

    task_acc_allsub(isub,:,:) = task_acc;
    task_rt_allsub(isub,:,:) = task_rt;

    %Performance as a function of target-probe distance
    tp_dist = abs(mod((WM_stim-WM_probe)+180, 360)-180);
    tp_dist(WM_stim==WM_probe) = 0;

    u_tp_dist = unique(tp_dist);
    for itask = 1:2
        for idist = 1:length(u_tp_dist)
            WM_acc_dist(isub,idist,itask) = mean(WM_acc(tp_dist==u_tp_dist(idist) & WM_task==itask));
            WM_rt_dist(isub,idist,itask) = median(WM_rt(tp_dist==u_tp_dist(idist) & WM_acc==1 & WM_task==itask));
        end
    end

    %Align all subjects to same boundary
    WM_stim_allsub(isub,:) = WM_stim;
    new_bound = 180; %align everything to boundary=180
    shift = new_bound - bound(1);
    WM_stim_rot = mod((WM_stim+shift), 360); %rotate and wrap between 0 and 360
    WM_stim_rot_allsub(isub,:) = WM_stim_rot;

    %Performance as a function of target-boundary distance
    tbound_dist_abs = abs(mod((WM_stim-bound(1))+90, 180)-90);
    tbound_dist_abs(WM_stim==bound(1)) = 0;
    tbound_dist_abs(WM_stim==bound(2)) = 0;

    u_bound_dist = unique(tbound_dist_abs);
    for itask = 1:2
        for idist = 1:length(u_bound_dist)
            WM_acc_bound(isub,idist,itask) = mean(WM_acc(tbound_dist_abs==u_bound_dist(idist) & WM_task==itask));
            WM_rt_bound(isub,idist,itask) = median(WM_rt(tbound_dist_abs==u_bound_dist(idist) & WM_acc==1 & WM_task==itask));
        end
    end

    %Follow-up analysis of surprising effect of target-boundary distance on
    %DMS performance
    exact_match = tp_dist==0;
    small_tp_dist = tp_dist>0 & tp_dist<=90;
    large_tp_dist = tp_dist>90;

    tbound_dist_dir = mod((WM_stim-bound(1))+180, 360)-180;
    pbound_dist_dir = mod((WM_probe-bound(1))+180, 360)-180;
    pbound_dist_abs = abs(mod((WM_probe-bound(1))+90, 180)-90);
    toward_trial = zeros(1,length(WM_stim));
    away_trial = zeros(1,length(WM_stim));
    for i_trial = 1:length(WM_stim)
        if small_tp_dist(i_trial)
            %Restrict analysis to trials where target and probe are within
            %the same quadrant
            if ((tbound_dist_dir(i_trial) >0 && tbound_dist_dir(i_trial)<90)...
                    && (pbound_dist_dir(i_trial) >0 && pbound_dist_dir(i_trial)<90))...
                    || (tbound_dist_dir(i_trial) <-90 && pbound_dist_dir(i_trial) <-90)
                toward_trial(i_trial) = (pbound_dist_dir(i_trial) > tbound_dist_dir(i_trial)) & (pbound_dist_abs(i_trial)>tbound_dist_abs(i_trial));
                away_trial(i_trial) = pbound_dist_dir(i_trial) < tbound_dist_dir(i_trial) & (pbound_dist_abs(i_trial)<tbound_dist_abs(i_trial));
            elseif (tbound_dist_dir(i_trial) >90 && pbound_dist_dir(i_trial) >90)...
                    || ((tbound_dist_dir(i_trial)<0 && tbound_dist_dir(i_trial)>-90)...
                    && (pbound_dist_dir(i_trial)<0 && pbound_dist_dir(i_trial)>-90))
                toward_trial(i_trial) = pbound_dist_dir(i_trial) < tbound_dist_dir(i_trial) & (pbound_dist_abs(i_trial)>tbound_dist_abs(i_trial));
                away_trial(i_trial) = pbound_dist_dir(i_trial) > tbound_dist_dir(i_trial) & (pbound_dist_abs(i_trial)<tbound_dist_abs(i_trial));
            end
        end
    end

    m2s_acc_toward(isub) = mean(WM_acc(WM_task==1 & small_tp_dist & toward_trial));
    m2s_acc_away(isub) = mean(WM_acc(WM_task==1 & small_tp_dist & away_trial));

    for itb_dist = 1:length(u_bound_dist)
        m2s_acc_bound_match(isub,itb_dist) = mean(WM_acc(tbound_dist_abs==u_bound_dist(itb_dist) & WM_task==1 & tp_dist==0));
        m2s_acc_bound_small(isub,itb_dist) = mean(WM_acc(tbound_dist_abs==u_bound_dist(itb_dist) & WM_task==1 & small_tp_dist));
        m2s_acc_bound_large(isub,itb_dist) = mean(WM_acc(tbound_dist_abs==u_bound_dist(itb_dist) & WM_task==1 & large_tp_dist));
    end

end

%% Plotting
cmap1 = cbrewer('qual', 'Set1', 12);
cmap2 = cbrewer('qual', 'Accent', 12);
cmap3 = cbrewer('div', 'PRGn', 12);
cmap4 = cbrewer('qual', 'Set2', 8);

%% Learning results (Figure 1)
%Learning across blocks
N = nsubs;
yMean_acc = squeeze(mean(block_acc,1));
ySEM_acc = squeeze(std(block_acc)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc = bsxfun(@times, ySEM_acc, CI95(:));

figure()
plot(unique_blocks, block_acc, 'k--'); hold on
errorbar(yMean_acc, yCI95_acc(2,:), 'color', 'k', 'linewidth', 2, 'MarkerSize', 20)
xticks(unique_blocks)
xlabel('Block number')
ylabel('Accuracy')
title('Accuracy over blocks')
savname = strcat(figpath, 'learning_acc_blocks.png');
saveas(gcf, savname)
savname = strcat(figpath, 'learning_acc_blocks.svg');
saveas(gcf, savname)

yMean_rt = squeeze(nanmean(block_rt,1));
ySEM_rt = squeeze(nanstd(block_rt)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt = bsxfun(@times, ySEM_rt, CI95(:));

figure()
plot(unique_blocks, block_rt, 'k--'); hold on
errorbar(yMean_rt, yCI95_rt(2,:), 'color', 'k', 'linewidth', 2, 'MarkerSize', 20)
xticks(unique_blocks)
xlabel('Block number')
ylabel('RT')
title('RT over blocks')
savname = strcat(figpath, 'learning_rt_blocks.png');
saveas(gcf, savname)
savname = strcat(figpath, 'learning_rt_blocks.svg');
saveas(gcf, savname)

% Stats
%accuracy
subjects = (1:size(block_acc, 1))';
tbl_acc = array2table(block_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6','Cond7','Cond8'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond8 ~ 1', 'WithinDesign', [1 2 3 4 5 6 7 8]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
tbl_rt = array2table(block_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6','Cond7','Cond8'});
tbl_rt.Subject = subjects;

% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond8 ~ 1', 'WithinDesign', [1 2 3 4 5 6 7 8]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%Learning relative to boundary 
yMean_acc = squeeze(mean(dist_acc,1));
ySEM_acc = squeeze(std(dist_acc)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc = bsxfun(@times, ySEM_acc, CI95(:));

figure()
plot(unique_dist, dist_acc, 'k--'); hold on
errorbar(unique_dist,yMean_acc, yCI95_acc(2,:), 'color', 'k', 'linewidth', 2, 'MarkerSize', 20)
xticks(unique_dist)
xlabel('distance to boundary')
ylabel('accuracy')
title('Accuracy - distance to boundary')
savname = strcat(figpath, 'learning_acc_boundary.png');
saveas(gcf, savname)
savname = strcat(figpath, 'learning_acc_boundary.svg');
saveas(gcf, savname)

yMean_rt = squeeze(mean(dist_rt,1));
ySEM_rt = squeeze(std(dist_rt)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt = bsxfun(@times, ySEM_rt, CI95(:));

figure()
plot(unique_dist, dist_rt, 'k--'); hold on
errorbar(unique_dist,yMean_rt, yCI95_rt(2,:), 'color', 'k', 'linewidth', 2, 'MarkerSize', 20)
xticks(unique_dist)
xlabel('distance to boundary')
ylabel('RT')
title('RT - distance to boundary')
savname = strcat(figpath, 'learning_rt_boundary.png');
saveas(gcf, savname)
savname = strcat(figpath, 'learning_rt_boundary.svg');
saveas(gcf, savname)

%Stats
%accuracy
subjects = (1:size(dist_acc, 1))';
tbl_acc = array2table(dist_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
tbl_rt = array2table(dist_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_rt.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%Pairwise comparisons
[h,p,ci,stat] = ttest(dist_acc(:,1), dist_acc(:,2))
es = meanEffectSize(dist_acc(:,1),dist_acc(:,2),Paired=true,Effect="cohen")

[h,p,ci,stat] = ttest(dist_acc(:,2), dist_acc(:,3))
es = meanEffectSize(dist_acc(:,2),dist_acc(:,3),Paired=true,Effect="cohen")

[h,p,ci,stat] = ttest(dist_rt(:,1), dist_rt(:,2))
es = meanEffectSize(dist_rt(:,1),dist_rt(:,2),Paired=true,Effect="cohen")

[h,p,ci,stat] = ttest(dist_rt(:,2), dist_rt(:,3))
es = meanEffectSize(dist_rt(:,2),dist_rt(:,3),Paired=true,Effect="cohen")

%% Overall WM performance (Figure 2AB)
acc_allsub = squeeze(task_acc_allsub);
rt_allsub = squeeze(task_rt_allsub);

N = nsubs;
yMean_acc = squeeze(mean(acc_allsub,1));
ySEM_acc = squeeze(std(acc_allsub)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc = bsxfun(@times, ySEM_acc, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(acc_allsub)), acc_allsub(:,1), 30, cmap2(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(acc_allsub)), acc_allsub(:,2), 30, cmap2(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_acc(1,1), yCI95_acc(2,1), '.', 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_acc(1,2), yCI95_acc(2,2), '.', 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [acc_allsub(:,1), acc_allsub(:,2)], 'color', [0.7,0.7,0.7])
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
ylabel('Accuracy','FontSize',20)
ylim([0.7 1])
title('Accuracy', 'FontSize', 20)
savname = strcat(figpath, 'behavior_overall_accuracy.png');
saveas(gcf, savname)
savname = strcat(figpath, 'behavior_overall_accuracy.svg');
saveas(gcf, savname)

N = nsubs;
yMean_rt = squeeze(mean(rt_allsub,1));
ySEM_rt = squeeze(std(rt_allsub)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt = bsxfun(@times, ySEM_rt, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(rt_allsub)), rt_allsub(:,1), 30, cmap2(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(rt_allsub)), rt_allsub(:,2), 30, cmap2(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_rt(1,1), yCI95_rt(2,1), '.', 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_rt(1,2), yCI95_rt(2,2), '.', 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [rt_allsub(:,1), rt_allsub(:,2)], 'color', [0.7,0.7,0.7])
xlim([0.5, 2.5])
ylim([0.4 1.2])
xticks([1,2])
xticklabels({'match-to-sample', 'match-to-category'})
ylabel('Reaction times (s)', 'FontSize',20)
title('Reaction times','FontSize',20)
savname = strcat(figpath, 'behavior_overall_rt.png');
saveas(gcf, savname)
savname = strcat(figpath, 'behavior_overall_rt.svg');
saveas(gcf, savname)

%Stats
[h p ci stat] = ttest(acc_allsub(:,1),acc_allsub(:,2))
es = meanEffectSize(acc_allsub(:,1),acc_allsub(:,2),Paired=true,Effect="cohen")

[h p ci stat] = ttest(rt_allsub(:,1),rt_allsub(:,2))
es = meanEffectSize(rt_allsub(:,1),rt_allsub(:,2),Paired=true,Effect="cohen")

%% WM performance: target-probe distance (Figure 2CD)
N = nsubs;
% Bin data for plotting and stats
WM_acc_dist_new = [];
% WM_acc_dist_new(:,1,:) = WM_acc_dist(:,1,:);
WM_acc_dist_new(:,1,:) = mean(WM_acc_dist(:,2:4,:),2);
WM_acc_dist_new(:,2,:) = mean(WM_acc_dist(:,5:7,:),2);
WM_acc_dist_new(:,3,:) = mean(WM_acc_dist(:,8:10,:),2);
WM_acc_dist_new(:,4,:) = mean(WM_acc_dist(:,11:13,:),2);
WM_acc_dist_new(:,5,:) = mean(WM_acc_dist(:,14:16,:),2);
WM_acc_dist_new(:,6,:) = mean(WM_acc_dist(:,17:19,:),2);

WM_rt_dist_new = [];
% WM_rt_dist_new(:,1,:) = WM_rt_dist(:,1,:);
WM_rt_dist_new(:,1,:) = nanmean(WM_rt_dist(:,2:4,:),2);
WM_rt_dist_new(:,2,:) = nanmean(WM_rt_dist(:,5:7,:),2);
WM_rt_dist_new(:,3,:) = nanmean(WM_rt_dist(:,8:10,:),2);
WM_rt_dist_new(:,4,:) = nanmean(WM_rt_dist(:,11:13,:),2);
WM_rt_dist_new(:,5,:) = nanmean(WM_rt_dist(:,14:16,:),2);
WM_rt_dist_new(:,6,:) = nanmean(WM_rt_dist(:,17:19,:),2);

u_tp_dist_new = [];
% u_tp_dist_new(1) = u_tp_dist(1);
u_tp_dist_new(1) = mean(u_tp_dist(2:4));
u_tp_dist_new(2) = mean(u_tp_dist(5:7));
u_tp_dist_new(3) = mean(u_tp_dist(8:10));
u_tp_dist_new(4) = mean(u_tp_dist(11:13));
u_tp_dist_new(5) = mean(u_tp_dist(14:16));
u_tp_dist_new(6) = mean(u_tp_dist(17:19));

%Accuracy
yMean_acc_m2s = squeeze(mean(WM_acc_dist_new(:,:,1),1));
ySEM_acc_m2s = squeeze(std(WM_acc_dist_new(:,:,1))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_m2s = bsxfun(@times, ySEM_acc_m2s, CI95(:));
ySEM_acc_m2s_match = squeeze(std(WM_acc_dist(:,1,1))./sqrt(N));
yCI95_acc_m2s_match = bsxfun(@times, ySEM_acc_m2s_match, CI95(:));

yMean_acc_m2c = squeeze(mean(WM_acc_dist_new(:,:,2),1));
ySEM_acc_m2c = squeeze(std(WM_acc_dist_new(:,:,2))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_m2c = bsxfun(@times, ySEM_acc_m2c, CI95(:));
ySEM_acc_m2c_match = squeeze(std(WM_acc_dist(:,1,2))./sqrt(N));
yCI95_acc_m2c_match = bsxfun(@times, ySEM_acc_m2c_match, CI95(:));

figure()
a = errorbar(u_tp_dist_new, yMean_acc_m2s, yCI95_acc_m2s(2,:), 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 20); hold on
b = errorbar(u_tp_dist_new, yMean_acc_m2c, yCI95_acc_m2c(2,:), 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 20);
errorbar(u_tp_dist(1), squeeze(mean(WM_acc_dist(:,1,1))), yCI95_acc_m2s_match(2,:), 'color', cmap2(1,:), "LineStyle","none", 'linewidth', 2, 'MarkerSize', 20); %plot match trials separately
errorbar(u_tp_dist(1), squeeze(mean(WM_acc_dist(:,1,2))), yCI95_acc_m2c_match(2,:), 'color', cmap2(2,:), "LineStyle","none", 'linewidth', 2, 'MarkerSize', 20); %plot match trials separately
legend([a,b],{'DMS', 'DMC'}, 'Location', 'best')
ylabel('Accuracy')
xlabel('Target-probe distance')
ylim([0 1])
xlim([-10 180])
savname = strcat(figpath, 'accuracy_targetProbeDist.png');
saveas(gcf, savname)
savname = strcat(figpath, 'accuracy_targetProbeDist.svg');
saveas(gcf, savname)

%RT
yMean_rt_m2s = squeeze(nanmean(WM_rt_dist_new(:,:,1),1));
ySEM_rt_m2s = squeeze(nanstd(WM_rt_dist_new(:,:,1))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt_m2s = bsxfun(@times, ySEM_rt_m2s, CI95(:));
ySEM_rt_m2s_match = squeeze(std(WM_rt_dist(:,1,1))./sqrt(N));
yCI95_rt_m2s_match = bsxfun(@times, ySEM_rt_m2s_match, CI95(:));

yMean_rt_m2c = squeeze(nanmean(WM_rt_dist_new(:,:,2),1));
ySEM_rt_m2c = squeeze(nanstd(WM_rt_dist_new(:,:,2))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt_m2c = bsxfun(@times, ySEM_rt_m2c, CI95(:));
ySEM_rt_m2c_match = squeeze(std(WM_rt_dist(:,1,2))./sqrt(N));
yCI95_rt_m2c_match = bsxfun(@times, ySEM_rt_m2c_match, CI95(:));

figure()
a=errorbar(u_tp_dist_new, yMean_rt_m2s, yCI95_rt_m2s(2,:), 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 20); hold on
b=errorbar(u_tp_dist_new, yMean_rt_m2c, yCI95_rt_m2c(2,:), 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 20);
errorbar(u_tp_dist(1), squeeze(mean(WM_rt_dist(:,1,1))), yCI95_rt_m2s_match(2,:), '.', 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 20); %plot match trials separately
errorbar(u_tp_dist(1), squeeze(mean(WM_rt_dist(:,1,2))), yCI95_rt_m2c_match(2,:), '.','color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 20); %plot match trials separately
legend([a,b],{'DMS', 'DMC'})
ylabel('Reaction times (s)')
xlabel('Target-probe distance')
xlim([-10 180])
savname = strcat(figpath, 'rt_targetProbeDist.png');
saveas(gcf, savname)
savname = strcat(figpath, 'rt_targetProbeDist.svg');
saveas(gcf, savname)

%Stats
%Match-to-sample
%Acc
m2s_dist_acc = squeeze(WM_acc_dist_new(:,:,1));
subjects = (1:size(m2s_dist_acc, 1))';
tbl_acc = array2table(m2s_dist_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond6 ~ 1', 'WithinDesign', [1 2 3 4 5 6]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
m2s_dist_rt = squeeze(WM_rt_dist_new(:,:,1));
subjects = (1:size(m2s_dist_rt, 1))';
tbl_rt = array2table(m2s_dist_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6'});
tbl_rt.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond6 ~ 1', 'WithinDesign', [1 2 3 4 5 6]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))


%Match-to-category
%Acc
m2c_dist_acc = squeeze(WM_acc_dist_new(:,:,2));
subjects = (1:size(m2c_dist_acc, 1))';
tbl_acc = array2table(m2c_dist_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond6 ~ 1', 'WithinDesign', [1 2 3 4 5 6]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
m2c_dist_rt = squeeze(WM_rt_dist_new(:,:,2));
subjects = (1:size(m2c_dist_rt, 1))';
tbl_rt = array2table(m2c_dist_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5', 'Cond6'});
tbl_rt.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond6 ~ 1', 'WithinDesign', [1 2 3 4 5 6]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%% WM performance: target-boundary distance (Figure 2EF)
N = nsubs;
% bin data for plotting and stats
WM_acc_bound_new = [];
WM_acc_bound_new(:,1,:) = mean(WM_acc_bound(:,1:3,:),2);
WM_acc_bound_new(:,2,:) = mean(WM_acc_bound(:,4:6,:),2);
WM_acc_bound_new(:,3,:) = mean(WM_acc_bound(:,7:9,:),2);

WM_rt_bound_new = [];
WM_rt_bound_new(:,1,:) = mean(WM_rt_bound(:,1:3,:),2);
WM_rt_bound_new(:,2,:) = mean(WM_rt_bound(:,4:6,:),2);
WM_rt_bound_new(:,3,:) = mean(WM_rt_bound(:,7:9,:),2);

u_bound_dist_new = [];
u_bound_dist_new(1) = mean(u_bound_dist(1:3));
u_bound_dist_new(2) = mean(u_bound_dist(4:6));
u_bound_dist_new(3) = mean(u_bound_dist(7:9));

%Accuracy
N = nsubs;
yMean_acc_m2s = squeeze(mean(WM_acc_bound_new(:,:,1),1));
ySEM_acc_m2s = squeeze(std(WM_acc_bound_new(:,:,1))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_m2s = bsxfun(@times, ySEM_acc_m2s, CI95(:));

yMean_acc_m2c = squeeze(mean(WM_acc_bound_new(:,:,2),1));
ySEM_acc_m2c = squeeze(std(WM_acc_bound_new(:,:,2))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_m2c = bsxfun(@times, ySEM_acc_m2c, CI95(:));

figure()
a = errorbar(u_bound_dist_new, yMean_acc_m2s, yCI95_acc_m2s(2,:), 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 20); hold on
b = errorbar(u_bound_dist_new, yMean_acc_m2c, yCI95_acc_m2c(2,:), 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 20);
legend([a,b],{'DMS', 'DMC'},'Location','best')
ylabel('Accuracy')
xlabel('Target-boundary distance')
ylim([0.7 1])
savname = strcat(figpath, 'accuracy_targetBoundDist.png');
saveas(gcf, savname)
savname = strcat(figpath, 'accuracy_targetBoundDist.svg');
saveas(gcf, savname)

%RT
yMean_rt_m2s = squeeze(nanmean(WM_rt_bound_new(:,:,1),1));
ySEM_rt_m2s = squeeze(nanstd(WM_rt_bound_new(:,:,1))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt_m2s = bsxfun(@times, ySEM_rt_m2s, CI95(:));

yMean_rt_m2c = squeeze(nanmean(WM_rt_bound_new(:,:,2),1));
ySEM_rt_m2c = squeeze(nanstd(WM_rt_bound_new(:,:,2))./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_rt_m2c = bsxfun(@times, ySEM_rt_m2c, CI95(:));

figure()
a=errorbar(u_bound_dist_new, yMean_rt_m2s, yCI95_rt_m2s(2,:), 'color', cmap2(1,:), 'linewidth', 2, 'MarkerSize', 20); hold on
b=errorbar(u_bound_dist_new, yMean_rt_m2c, yCI95_rt_m2c(2,:), 'color', cmap2(2,:), 'linewidth', 2, 'MarkerSize', 20);
legend([a,b],{'m2sample', 'm2category'})
ylabel('Reaction times (s)')
xlabel('Target-boundary distance')
savname = strcat(figpath, 'rt_targetBoundDist.png');
saveas(gcf, savname)
savname = strcat(figpath, 'rt_targetBoundDist.svg');
saveas(gcf, savname)

%Stats
%Match-to-sample
%Acc
m2s_bound_acc = squeeze(WM_acc_bound_new(:,:,1));
subjects = (1:size(m2s_bound_acc, 1))';
tbl_acc = array2table(m2s_bound_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
m2s_bound_rt = squeeze(WM_rt_bound_new(:,:,1));
subjects = (1:size(m2s_bound_rt, 1))';
tbl_rt = array2table(m2s_bound_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_rt.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))


%Match-to-category
%Acc
m2c_bound_acc = squeeze(WM_acc_bound_new(:,:,2));
subjects = (1:size(m2c_bound_acc, 1))';
tbl_acc = array2table(m2c_bound_acc, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%RT
m2c_bound_rt = squeeze(WM_rt_bound_new(:,:,2));
subjects = (1:size(m2c_bound_rt, 1))';
tbl_rt = array2table(m2c_bound_rt, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_rt.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_rt, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%% Follow-up analysis of DMS category effect (Figure 2GH)

m2s_acc_bound_match_new = [];
m2s_acc_bound_match_new(:,1,:) = mean(m2s_acc_bound_match(:,1:3),2);
m2s_acc_bound_match_new(:,2,:) = mean(m2s_acc_bound_match(:,4:6),2);
m2s_acc_bound_match_new(:,3,:) = mean(m2s_acc_bound_match(:,7:9),2);

m2s_acc_bound_small_new = [];
m2s_acc_bound_small_new(:,1,:) = nanmean(m2s_acc_bound_small(:,1:3),2);
m2s_acc_bound_small_new(:,2,:) = nanmean(m2s_acc_bound_small(:,4:6),2);
m2s_acc_bound_small_new(:,3,:) = nanmean(m2s_acc_bound_small(:,7:9),2);

m2s_acc_bound_large_new = [];
m2s_acc_bound_large_new(:,1,:) = mean(m2s_acc_bound_large(:,1:3),2);
m2s_acc_bound_large_new(:,2,:) = mean(m2s_acc_bound_large(:,4:6),2);
m2s_acc_bound_large_new(:,3,:) = mean(m2s_acc_bound_large(:,7:9),2);

u_bound_dist_new = [];
u_bound_dist_new(1) = mean(u_bound_dist(1:3));
u_bound_dist_new(2) = mean(u_bound_dist(4:6));
u_bound_dist_new(3) = mean(u_bound_dist(7:9));

%Accuracy
N = nsubs;
yMean_acc_match = squeeze(mean(m2s_acc_bound_match_new,1));
ySEM_acc_match = squeeze(std(m2s_acc_bound_match_new)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_match = bsxfun(@times, ySEM_acc_match, CI95(:));

yMean_acc_small = squeeze(mean(m2s_acc_bound_small_new,1));
ySEM_acc_small = squeeze(std(m2s_acc_bound_small_new)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_small = bsxfun(@times, ySEM_acc_small, CI95(:));

yMean_acc_large = squeeze(mean(m2s_acc_bound_large_new,1));
ySEM_acc_large = squeeze(std(m2s_acc_bound_large_new)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_large = bsxfun(@times, ySEM_acc_large, CI95(:));

new_map = cbrewer('seq', 'YlGn', 12);
figure()
a = errorbar(u_bound_dist_new, yMean_acc_match, yCI95_acc_match(2,:), 'color', new_map(4,:), 'linewidth', 2, 'MarkerSize', 20); hold on
b = errorbar(u_bound_dist_new, yMean_acc_small, yCI95_acc_small(2,:), 'color', new_map(8,:), 'linewidth', 2, 'MarkerSize', 20);
c = errorbar(u_bound_dist_new, yMean_acc_large, yCI95_acc_large(2,:), 'color', new_map(12,:), 'linewidth', 2, 'MarkerSize', 20);
legend([a,b,c],{'match', 'small', 'large'},'Location','best')
ylabel('Accuracy')
xlabel('Target-boundary distance')
savname = strcat(figpath, 'm2s_acc_tpdist_bounddist.png');
saveas(gcf, savname)
savname = strcat(figpath, 'm2s_acc_tpdist_bounddist.svg');
saveas(gcf, savname)

%Match-to-sample
y_dist = m2s_acc_bound_small_new;
% y_dist   = m2s_acc_bound_match_new;
% y_dist   = m2s_acc_bound_large_new;
subjects = (1:size(y_dist, 1))';
tbl_acc = array2table(y_dist, 'VariableNames', {'Cond1', 'Cond2', 'Cond3'});
tbl_acc.Subject = subjects;
% Fit a repeated measures model
rm = fitrm(tbl_acc, 'Cond1-Cond3 ~ 1', 'WithinDesign', [1 2 3]);

% Perform repeated measures ANOVA
ranovatbl = ranova(rm);
mauchlytbl = mauchly(rm);

% Display ANOVA results
disp(ranovatbl);
disp(mauchlytbl);
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(end))

%plot
yMean_acc_toward = squeeze(mean(m2s_acc_toward',1));
ySEM_acc_toward = squeeze(std(m2s_acc_toward)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_toward = bsxfun(@times, ySEM_acc_toward, CI95(:));

yMean_acc_away = squeeze(mean(m2s_acc_away',1));
ySEM_acc_away = squeeze(std(m2s_acc_away)./sqrt(N));
CI95 = tinv([0.025 0.975], N-1);
yCI95_acc_away = bsxfun(@times, ySEM_acc_away, CI95(:));

figure('Color', 'white', 'Position', [100 100 300 300])
ax1=scatter(1.1*ones(1,length(m2s_acc_toward)), m2s_acc_toward, 30, cmap1(1,:), 'filled'); hold on
ax2=scatter(1.9*ones(1,length(m2s_acc_away)), m2s_acc_away, 30, cmap1(2,:), 'filled');
ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5;
errorbar(1, yMean_acc_toward, yCI95_acc_toward(2), '.', 'color', cmap1(1,:), 'linewidth', 2, 'MarkerSize', 30)
errorbar(2, yMean_acc_away, yCI95_acc_away(2), '.', 'color', cmap1(2,:), 'linewidth', 2, 'MarkerSize', 30)
plot([1.2,1.8], [m2s_acc_toward; m2s_acc_away], 'color', [0.7,0.7,0.7])
xlim([0.5, 2.5])
xticks([1,2])
xticklabels({'toward', 'away'})
title('all target-boundary distances')

%Stats
[h p ci stat] = ttest(m2s_acc_toward, m2s_acc_away)
es = meanEffectSize(m2s_acc_toward,m2s_acc_away,Paired=true,Effect="cohen")
