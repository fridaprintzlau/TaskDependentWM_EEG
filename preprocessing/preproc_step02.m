%% EEG preprocessing - step 2: identify artefacts
% 'Categorisation dynamically reconfigures working memory in a 
% task-dependent manner' by Frida Printzlau, Olya Bulatova, Michael Mack 
% & Keisuke Fukuda

% - Frida Printzlau, October 2024

%Dependencies:
%-- Fieldtrip toolbox available for free download here: https://www.fieldtriptoolbox.org/
%   Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs
%   Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, 
%   EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, 2011; 2011:156869.
%-- Matlab 'Symbolic Math Toolbox', 'Signal Processing Toolbox', 'Image Processing Toolbox', 
%   'Statistics and Machine Learning toolbox' & 'Parallel Computing toolbox'

clear all;
close all;
clc;

studypath = ''; %define path where the data folder sits
fieldtrippath = ''; %add fieldtrip path here
addpath(fieldtrippath);
ft_defaults;
fs = filesep;

dat_dir = [studypath strcat(fs,'Data/EEG/Preprocessed')];

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs   = length(sublist);

%%
window_size = 50; %number of samples 50= 200ms
step_size = 25;
eMove_threshold = 30;%in microvolts
for isub = 1:nsubs
    %Load filtered, epoched data
    substrg  = sprintf('S%02d%s',sublist(isub));
    subpath  = sprintf('%s%s%s',dat_dir,fs, substrg);    
    filename = sprintf('%s%s%s_ftbeforeAR.mat',subpath,fs,substrg);
    if ~exist(filename,'file')
        error(sprintf('%s does not exist!',filename));
    end
    load(filename);

    % Identify eye artefacts
    cfg = [];
    cfg.artfctdef.eog.channel = {'VEOG'};
    cfg.artfctdef.eog.bpfilter   = 'yes';
    cfg.artfctdef.eog.bpfilttype = 'but';
    cfg.artfctdef.eog.bpfreq     = [1 15];
    cfg.artfctdef.eog.bpfiltord  = 4;
    cfg.artfctdef.eog.hilbert    = 'yes';
    % cfg.artfctdef.eog.interactive = 'yes';
    cfg = ft_artifact_eog(cfg, data);
    eog.artfctdef = cfg.artfctdef;

    %Identify other major artefacts based on variance
    cfg = [];
    cfg.method = 'summary';
    cfg.channel = 'eeg';
    cfg.keepchannel = 'yes';
    cfg.latency = [-0.3 2.2];
    cfg.keeptrial   = 'nan'; %keep the trials order/numbering, just fill with nan (important for matching to behaviour)
    cfg.keepchannel = 'nan';
    data_summary = ft_rejectvisual(cfg, data);

    %Visual inspection of marked artefacts
    cfg = [];
    % cfg.channel = 'EOG'
    cfg.ylim = 'maxmin';
    cfg.viewmode = 'vertical';
    cfg.allowoverlap = 'yes';
    cfg.artfctdef.eog = eog.artfctdef.eog; % reuse the automatic artifacts from above
    cfg.artfctdef.summary   = data_summary.cfg.artfctdef.summary;
    artf = ft_databrowser(cfg, data);

    dat = [];
    ntrials = size(data.trial,2);
    for itrial = 1:ntrials
        dat(itrial,:,:) = data.trial{itrial};
    end

    %Saccade detection
    disp('detecting eMove...');
    eMove= zeros(ntrials, size(dat,3));
    %first identify rHEOG and lHEOG
    for chan = 1:length(data.label)
        curr_chan_name = char(data.label(chan,:));
        if isequal(curr_chan_name,'L_HEOG')==1
            lHEOG_channel_index = chan;
        elseif isequal(curr_chan_name,'R_HEOG')==1
            rHEOG_channel_index = chan;
        end
    end

    for itrial=1:ntrials
        eMove_data = squeeze(dat(itrial,lHEOG_channel_index,:)-dat(itrial,rHEOG_channel_index,:));%lHEOG-rHEOG
        t_init = 1;
        while 1
            t_end = t_init+window_size;
            current_data = eMove_data(t_init:t_end);
            window_center = round(window_size/2);
            pre_amp = mean(current_data(1:window_center-1));
            post_amp = mean(current_data(window_center:length(current_data)));

            if abs(pre_amp-post_amp)>eMove_threshold
                eMove(itrial,t_init:t_end) = 1;
            end
            t_init = t_init + step_size;
            if t_init + window_size > length(eMove_data)
                break
            end
        end
    end

    %Inspect and exclude saccade trials
    sacc_toi = data.time{1}>=-0.2 & data.time{1}<=2.2; %full trial
    sacc_trials = find(sum(eMove(:,sacc_toi),2));
    if ~isempty(sacc_trials)
        cfg = [];
        cfg.trials = sacc_trials;
        sacc_data = ft_selectdata(cfg, data);
        cfg = [];
        cfg.ylim = 'maxmin';
        cfg.viewmode = 'vertical';
        cfg.allowoverlap = 'yes';
        cfg.channel = {'R_HEOG', 'L_HEOG'};
        sacc_data = ft_databrowser(cfg, sacc_data);
    else, disp('no saccades')
    end
    rejsaccade = zeros(ntrials,1);
    rejsaccade(sacc_trials) = 1;

    % Exclude blink trials
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.eog.artifact = artf.artfctdef.eog.artifact;
    clean_data = ft_rejectartifact(cfg, data);

    stim_toi = [-0.2 0.2]; %stimulus + mask
    ping_toi = [1.7 1.8]; %ping
    toi_indx = clean_data.time{1}>=stim_toi(1) & clean_data.time{1}<=stim_toi(2) | clean_data.time{1}>=ping_toi(1) & clean_data.time{1}<=ping_toi(2);
    dat_full = []; dat_toi = [];
    for itrial = 1:length(clean_data.trial)
        dat_toi(itrial,:,:) = cat(3,clean_data.trial{itrial}(:,toi_indx));
    end
    rejblink = any(isnan(reshape(dat_toi,[size(dat_toi,1) size(dat_toi,2)*size(dat_toi,3)])),2); %reject any trial with a marked artifact in critical periods

    fprintf('percent blink trials: %f', sum(rejblink)/ntrials*100);
    fprintf('\npercent saccade trials total: %f', length(sacc_trials)/ntrials*100);

    %Exclude summary trials
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.summary.artifact   = artf.artfctdef.summary.artifact;
    clean_data = ft_rejectartifact(cfg, data);

    full_toi = [-0.2 2.2]; %full trial
    time_indx = clean_data.time{1}>=full_toi(1) & clean_data.time{1}<=full_toi(2);
    dat = [];
    for itrial = 1:length(clean_data.trial)
        dat(itrial,:,:) = cat(3,clean_data.trial{itrial}(:,time_indx));
    end
    rejother = any(isnan(reshape(dat,[size(dat,1) size(dat,2)*size(dat,3)])),2); %reject any trial with other marked artifact

    artfct_def_eye   = rejblink | rejsaccade;
    artfct_def_other = rejother;
    
    %Save marked artefact trials for each subject
    rejfile     = sprintf('%s%s%s_artefacts.mat',subpath,fs,substrg);
    save(rejfile,'artfct_def_eye','artfct_def_other');
    
    clear dat data data_bp
end
