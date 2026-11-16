%% RSA analysis
% Figure 4 & 5 in 'Categorisation dynamically
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
figpath = ''; %define where to save the figures ;

run_analysis = false; %run analysis or load saved results

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs = length(sublist);
%% Create model RDMs

% Fit models separately for the two tasks
%Category coding model
catmodel = squareform([ones(6)*0 ones(6)*1;...
    ones(6)*1 ones(6)*0]);

%Distance to boundary model
dist = [0 0.5 1 1 0.5 0; 0.5 0 0.5 0.5 0 0.5; 1 0.5 0 0 0.5 1; 1 0.5 0 0 0.5 1; 0.5 0 0.5 0.5 0 0.5; 0 0.5 1 1 0.5 0];
disttobound_model = squareform([dist dist; dist dist]);

%Stimulus model
m1 = [0 0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667;... 
    0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333;...
    0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667 0.5;...
    0.5 0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667;...
    0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333;...
    0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333 1;...
    1 0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667 0.8333;...
    0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333 0.5 0.6667;...
    0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333 0.5;...
    0.5 0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667 0.3333;...
    0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667 0 0.1667;...
    0.1667 0.3333 0.5 0.6667 0.8333 1 0.8333 0.6667 0.5 0.3333 0.1667 0];
stimmodel = squareform(m1);

X = [catmodel' stimmodel' disttobound_model'];
X = zscore(X);

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

        %Color/Shape bins
        bin_size = 30; % 30 for 12 bins between 1-360 deg.
        bin_angles_deg = 0:bin_size:360;
        nangs = length(bin_angles_deg) - 1;
        bin_angles = bin_angles_deg;
        nsets=1;
        bin_angles = [zeros(nsets,1), bin_angles];
        bin_angles_rad = bin_angles.*pi/180;

        %Rotate stimulus values according to category
        new_bound = 180; %align everything to boundary=180
        shift = new_bound - new_behav.bound(1);
        new_behav.target_rot = mod((new_behav.target+shift), 360); %rotate and wrap between 0 and 360

        %Get target labels & assign to bins
        target_ang = new_behav.target_rot.*pi/180; %convert to radians
        [N, edges, bins] = histcounts(target_ang, bin_angles_rad);
        bins=bins-1;

        labels = bin_angles_deg.*pi/180;
        binned_shapes = NaN(nsets, ntrials);
        for i_bin = 1:length(labels)
            binned_shapes(bins == i_bin) = labels(i_bin);
        end

        %Baseline correction
        onset = 0;
        baseline_window = time_orig>=-0.2 & time_orig<=-0.1;
        dat = bsxfun(@minus,dat,mean(dat(:,:,baseline_window),3));

        %Which channels to use?
        eeg_chans = find(ismember(clean_data.label,{'AF3', 'F3', 'Fz', 'FC5', 'FC1' 'C3' 'T7', 'P7', 'P3',...
            'Pz', 'POz', 'PO3', 'PO7', 'O1', 'O2', 'PO4', 'PO8', 'P4', 'P8', 'T8', 'C4', 'Cz', 'FC2', 'FC6', 'F4', 'AF4', 'AFz'}));
        chans = eeg_chans;
        %% Compute RDM and fit models

        stim = binned_shapes;
        u_stim = unique(stim);
        task = new_behav.task_type;
        u_task = unique(task);

        condxchanxtime = nan(length(u_stim)*length(u_task), length(chans), length(time_orig));
        resids = nan(ntrials, length(chans), length(time_orig));
        for itask = 1:length(u_task)
            for istim = 1:length(u_stim)
                curr_trls = stim==u_stim(istim) & task==itask;
                curr_cond = ((itask-1)*length(u_stim))+istim;
                resids(curr_trls,:,:) = dat(curr_trls,chans,:)-mean(dat(curr_trls,chans,:),1); %subtract condition mean
                condxchanxtime(curr_cond,:,:) = mean(dat(curr_trls,chans,:),1); %take within condition average
            end
        end
        timexchanxcond = permute(condxchanxtime, [3 2 1]);

        for t = 1:length(time_orig) %For each time point
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

            %Fit models separately for each task
            for itask = 1:length(u_task)
                distforfit = squareform(stimDistances);
                distforfit = squareform(distforfit(((itask-1)*length(u_stim))+1:((itask-1)*length(u_stim))+12,...
                    ((itask-1)*length(u_stim))+1:((itask-1)*length(u_stim))+12)); %do fit separately for each task
                b = glmfit(X,distforfit','normal','constant','on');
                RDM_fit_glm(isub,itask,t,:) = b;
            end
        end
    end

    %% Save data
    savfile = sprintf('%s%sAnalysis_output/RDM_results', datapath, fs);
    save(savfile, 'RDM_output', 'RDM_output_nozscore', 'RDM_fit_glm', 'time_orig');
else
    %Load data
    load(sprintf('%s%sAnalysis_output/RDM_results.mat', datapath,fs))
end
%% Prepare for plotting
colorpath = '/Volumes/FridaPrintz/IncidentalWM_EEG_Exp1/scripts/toolbox/cbrewer';
addpath(genpath(colorpath));
cmap1 = cbrewer('div', 'PuOr', 12); %
cmap2 = cbrewer('seq', 'Greens', 12);
cmap3 = cbrewer('qual', 'Set1', 12);
cmap4 = cbrewer('qual', 'Accent', 12);
cmap5 = cbrewer('div', 'PiYG', 12);
cmap6 = cbrewer('qual', 'Set2', 8);
cmap7 = cbrewer('div', 'BrBG', 8);

%% Plot RDMs (Figure 4A)

%Stimulus onset
toi = time_orig>=0 & time_orig<=0.2; 
plotRDM = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));

figure()
colormap('bone');
imagesc(plotRDM,[0.9 1.35]);
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Data RDM (0-0.2 s)');
xticks(1:24)
yticks(1:24)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'dataRDM_stimTOI.png');
saveas(gcf, savname)
savname = strcat(figpath, 'dataRDM_stimTOI.svg');
saveas(gcf, savname)

%Early delay
toi = time_orig>=0.2 & time_orig<=0.7; 
plotRDM = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));

figure()
colormap(bone);
imagesc(plotRDM,[0.9 1.35]);
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Data RDM (0.2-0.7 s)');
xticks(1:24)
yticks(1:24)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'dataRDM_earlyDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'dataRDM_earlyDelay.svg');
saveas(gcf, savname)

%Late delay
toi = time_orig>=0.7 & time_orig<=1.7; 
plotRDM = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));

figure()
colormap(bone);
imagesc(plotRDM,[0.9 1.35]);
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Data RDM (0.7-1.7 s)');
xticks(1:24)
yticks(1:24)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B','1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'dataRDM_lateDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'dataRDM_lateDelay.svg');
saveas(gcf, savname)

%% MDS plots (Figure 4B)

%Encoding 
toi = time_orig>=0 & time_orig<=0.2; 
RDM_for_MDS = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));
s = [30 90 160 160 90 30 30 90 160 160 90 30];

%2D
[Y,e] = cmdscale(RDM_for_MDS,2);
figure('Color', 'white', 'Position', [100 100 450 400])
p1 = plot(Y(1:12,1),Y(1:12,2), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter(Y(1:12,1),Y(1:12,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot(Y(13:24,1),Y(13:24,2), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter(Y(13:24,1),Y(13:24,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Stim onset (0-0.2)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
savname = strcat(figpath, 'MDS_2D_stimOnset.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_2D_stimOnset.svg');
saveas(gcf, savname)

%3D
[Y,e] = cmdscale(RDM_for_MDS,3);
figure('Color', 'white', 'Position', [100 100 450 400],'Renderer','Painters')
p1 = plot3(Y(1:12,1),Y(1:12,2),Y(1:12,3), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter3(Y(1:12,1),Y(1:12,2),Y(1:12,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot3(Y(13:24,1),Y(13:24,2),Y(13:24,3), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter3(Y(13:24,1),Y(13:24,2),Y(13:24,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Stim onset (0-0.2)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
zlabel('Dimension 3','FontSize',12)
savname = strcat(figpath, 'MDS_3D_stimOnset.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_3D_stimOnset.svg');
saveas(gcf, savname)

%Early Delay
toi = time_orig>=0.2 & time_orig<=0.7; 
RDM_for_MDS = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));
s = [30 90 160 160 90 30 30 90 160 160 90 30];

%2D
[Y,e] = cmdscale(RDM_for_MDS,2);
figure('Color', 'white', 'Position', [100 100 450 400])
p1 = plot(Y(1:12,1),Y(1:12,2), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter(Y(1:12,1),Y(1:12,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot(Y(13:24,1),Y(13:24,2), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter(Y(13:24,1),Y(13:24,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Delay (0.2-0.7)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
savname = strcat(figpath, 'MDS_2D_earlyDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_2D_earlyDelay.svg');
saveas(gcf, savname)

%3D
[Y,e] = cmdscale(RDM_for_MDS,3);
figure('Color', 'white', 'Position', [100 100 450 400],'Renderer','Painters')
p1 = plot3(Y(1:12,1),Y(1:12,2),Y(1:12,3), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter3(Y(1:12,1),Y(1:12,2),Y(1:12,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot3(Y(13:24,1),Y(13:24,2),Y(13:24,3), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter3(Y(13:24,1),Y(13:24,2),Y(13:24,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Delay (0.2-0.7)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
zlabel('Dimension 3','FontSize',12)
savname = strcat(figpath, 'MDS_3D_earlyDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_3D_earlyDelay.svg');
saveas(gcf, savname)

%Late delay
toi = time_orig>=0.7 & time_orig<=1.7; 
RDM_for_MDS = squeeze(mean(mean(RDM_output_nozscore(:,:,:,toi),4),1));
s = [30 90 160 160 90 30 30 90 160 160 90 30];

%2D (dimension 1 & 2)
[Y,e] = cmdscale(RDM_for_MDS,2);
figure('Color', 'white', 'Position', [100 100 450 400])
p1 = plot(Y(1:12,1),Y(1:12,2), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter(Y(1:12,1),Y(1:12,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot(Y(13:24,1),Y(13:24,2), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter(Y(13:24,1),Y(13:24,2), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Late delay (0.7-1.7 s)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
savname = strcat(figpath, 'MDS_2D_lateDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_2D_lateDelay.svg');
saveas(gcf, savname)

% Dimension 1 & 3
[Y,e] = cmdscale(RDM_for_MDS,3);
figure('Color', 'white', 'Position', [100 100 450 400])
p1 = plot(Y(1:12,1),Y(1:12,3), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter(Y(1:12,1),Y(1:12,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot(Y(13:24,1),Y(13:24,3), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter(Y(13:24,1),Y(13:24,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Later delay (0.7-1.7 s)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 3','FontSize',12)
savname = strcat(figpath, 'MDS_2D_D1_3_lateDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_2D_D1_3_lateDelay.svg');
saveas(gcf, savname)

% Dimension 2 & 3
figure('Color', 'white', 'Position', [100 100 450 400])
p1 = plot(Y(1:12,2),Y(1:12,3), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter(Y(1:12,2),Y(1:12,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot(Y(13:24,2),Y(13:24,3), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter(Y(13:24,2),Y(13:24,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Later delay (0.7-1.7 s)')
xlabel('Dimension 2','FontSize',12)
ylabel('Dimension 3','FontSize',12)
savname = strcat(figpath, 'MDS_2D_D2_3_lateDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_2D_D2_3_lateDelay.svg');
saveas(gcf, savname)

%3D
[Y,e] = cmdscale(RDM_for_MDS,3);
figure('Color', 'white', 'Position', [100 100 450 400],'Renderer','Painters')
p1 = plot3(Y(1:12,1),Y(1:12,2),Y(1:12,3), 'LineStyle', '-', 'marker', 'None', 'color', 'k'); hold on
scatter3(Y(1:12,1),Y(1:12,2),Y(1:12,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
p2 = plot3(Y(13:24,1),Y(13:24,2),Y(13:24,3), 'LineStyle', '--', 'marker', 'None', 'color', 'k');
scatter3(Y(13:24,1),Y(13:24,2),Y(13:24,3), s, cmap1, 'filled', 'MarkerEdgeColor', 'k')
title('Late delay (0.7-1.7 s)')
xlabel('Dimension 1','FontSize',12)
ylabel('Dimension 2','FontSize',12)
zlabel('Dimension 3','FontSize',12)
savname = strcat(figpath, 'MDS_3D_lateDelay.png');
saveas(gcf, savname)
savname = strcat(figpath, 'MDS_3D_lateDelay.svg');
saveas(gcf, savname)

%% Plot model fits (Figure 5)

%Plot models
cmap7 = cbrewer('div', 'BrBG', 12);
figure('Color', 'white', 'Position', [100 100 500 400])
colormap gray
imagesc(squareform(stimmodel)); 
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Stimulus model');
xticks(1:12)
yticks(1:12)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'modelRDM_stimModel.png');
saveas(gcf, savname)
savname = strcat(figpath, 'modelRDM_stimModel.svg');
saveas(gcf, savname)

figure('Color', 'white', 'Position', [100 100 500 400])
colormap gray
imagesc(squareform(catmodel)); 
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Category model');
xticks(1:12)
yticks(1:12)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'modelRDM_catModel.png');
saveas(gcf, savname)
savname = strcat(figpath, 'modelRDM_catModel.svg');
saveas(gcf, savname)

figure('Color', 'white', 'Position', [100 100 500 400])
colormap gray
imagesc(squareform(disttobound_model)); 
ax = gca;
ax.TickDir = 'in';
cbh=colorbar('v');
title('Distance-to-boundary model');
xticks(1:12)
yticks(1:12)
xticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
yticklabels({'1A','2A','3A','4A','5A','6A','7B','8B','9B','10B','11B','12B'})
xlabel('stimulus','FontSize',12)
ylabel('stimulus','FontSize',12)
savname = strcat(figpath, 'modelRDM_disttoboundModel.png');
saveas(gcf, savname)
savname = strcat(figpath, 'modelRDM_disttoboundModel.svg');
saveas(gcf, savname)


cmap7 = cbrewer('div', 'BrBG', 8);
% Plot model fits
for itask = 1:2
    toi = time_orig>=-0.5 & time_orig<=4;
    plotcondition = 2:4;
    plottime = time_orig(toi);
    plotdata = squeeze(RDM_fit_glm(:,itask,toi,plotcondition));

    %smooth data
    fsample = round(median(1./diff(time_orig)));
    ftype = 'gaussian';
    fsize = 0.04*fsample;
    plotdata_smooth = filtfast(plotdata,2,[],ftype,fsize);
    plotcol = [cmap7(1,:); cmap7(end,:); cmap5(1,:); [0,0,0]];

    plotsubs = 1:nsubs;
    corrwin  = [0 2.2];
    nperm         = 10000; %number of permutations (set to 0 for uncorrected pvalues)
    pthresh       = 0.15; %threshold for plotting significance bar
    timeseries = time_orig;

    [decod_b, decod_stats] = computeStats(plotdata_smooth, corrwin, timeseries, nperm, pthresh);

    %Plot time-resolved decoding quality and cluster-corrected stats
    figure('Color', 'white', 'Position', [100 100 600 400])
    statcols = plotcol;
    ylims = [-0.1 0.25];
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
    hplot = plotpatch(plotdata_smooth,plottime,plotcol,cfg);
    plot(time_orig, zeros(size(time_orig)), 'k--')
    line([0, 0], [ylims(1), ylims(2)], 'color', 'k')
    line([0.15, 0.15], [ylims(1), ylims(2)], 'color', 'k')
    line([1.7, 1.7], [ylims(1), ylims(2)], 'color', 'k')
    line([2.2, 2.2], [ylims(1), ylims(2)], 'color', 'k')
    legend(hplot, {'category coding', 'stimulus coding', 'dist-to-bound coding'})
    ylim(ylims)
    xlim([-0.2 2.2])
    switch itask
        case 1
            savname = strcat(figpath, 'RSA_modelfit_m2s.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'RSA_modelfit_m2s.svg');
            saveas(gcf, savname)
        case 2
            savname = strcat(figpath, 'RSA_modelfit_m2c.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'RSA_modelfit_m2c.svg');
            saveas(gcf, savname)
    end

    toi = time_orig>=0.7 & time_orig<=1.7;
    task_toi = squeeze(mean(RDM_fit_glm(:,itask,toi,plotcondition),3));

    N = nsubs;
    yMean_toi = squeeze(mean(task_toi,1));
    ySEM_toi = squeeze(std(task_toi)./sqrt(N));
    CI95 = tinv([0.025 0.975], N-1);
    yCI95_toi = bsxfun(@times, ySEM_toi, CI95(:));

    figure('Color', 'white', 'Position', [100 100 300 300])
    ax1=scatter(1.1*ones(1,length(task_toi(:,1))), task_toi(:,1), 30, cmap7(1,:), 'filled'); hold on
    ax2=scatter(2.1*ones(1,length(task_toi(:,2))), task_toi(:,2), 30, cmap7(end,:), 'filled');
    ax3=scatter(3.1*ones(1,length(task_toi(:,3))), task_toi(:,3), 30, cmap5(1,:), 'filled');
    ax1.MarkerFaceAlpha=0.5; ax2.MarkerFaceAlpha=0.5; ax3.MarkerFaceAlpha=0.5;
    errorbar(1, yMean_toi(1,1), yCI95_toi(2,1), '.', 'color', cmap7(1,:), 'linewidth', 2, 'MarkerSize', 30)
    errorbar(2, yMean_toi(1,2), yCI95_toi(2,2), '.', 'color', cmap7(end,:), 'linewidth', 2, 'MarkerSize', 30)
    errorbar(3, yMean_toi(1,3), yCI95_toi(2,3), '.', 'color', cmap5(1,:), 'linewidth', 2, 'MarkerSize', 30)
    plot([0.5:3.5], zeros(size([0.5:3.5])), 'k--')
    xlim([0.5, 3.5])
    ylim([-0.2 0.5])
    xticks([1,2,3])
    xticklabels({'Category coding', 'Stimulus coding', 'Dist-to-bound coding'})
    title('RSA model fits (late 0.7-1.7 s)')
    switch itask
        case 1
            savname = strcat(figpath, 'RSA_modelfit_lateTOI_m2s.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'RSA_modelfit_lateTOI_m2s.svg');
            saveas(gcf, savname)
        case 2
            savname = strcat(figpath, 'RSA_modelfit_lateTOI_m2c.png');
            saveas(gcf, savname)
            savname = strcat(figpath, 'RSA_modelfit_lateTOI_m2c.svg');
            saveas(gcf, savname)
    end

    switch itask
        case 1, [h,p,stats,ci]=ttest(task_toi(:,1),0) %do two tailed for m2s
        case 2, [h,p,stats,ci]=ttest(task_toi(:,1),0, 'tail', 'right')
    end
    es = meanEffectSize(task_toi(:,1),Effect="cohen")
    [h,p,stats,ci]=ttest(task_toi(:,2),0,'tail','right')
    es = meanEffectSize(task_toi(:,2),Effect="cohen")
    [h,p,stats,ci]=ttest(task_toi(:,3),0,'tail','right')
    es = meanEffectSize(task_toi(:,3),Effect="cohen")
end

%Compare tasks
toi = time_orig>=0.7 & time_orig<=1.7;
cat_code_toi = squeeze(mean(RDM_fit_glm(:,:,toi,2),3));
stim_code_toi = squeeze(mean(RDM_fit_glm(:,:,toi,3),3));
dist_code_toi = squeeze(mean(RDM_fit_glm(:,:,toi,4),3));

[h,p,stats,ci]=ttest(cat_code_toi(:,1),cat_code_toi(:,2))
es = meanEffectSize(cat_code_toi(:,1),cat_code_toi(:,2), 'Paired', true, Effect="cohen")
[h,p,stats,ci]=ttest(stim_code_toi(:,1),stim_code_toi(:,2))
es = meanEffectSize(stim_code_toi(:,1),stim_code_toi(:,2), 'Paired', true, Effect="cohen")
[h,p,stats,ci]=ttest(dist_code_toi(:,1),dist_code_toi(:,2))
es = meanEffectSize(dist_code_toi(:,1),dist_code_toi(:,2), 'Paired', true, Effect="cohen")
