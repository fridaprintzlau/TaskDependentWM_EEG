%% EEG preprocessing - step 3: Reject artefacts and merge EEG and behavioural data
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

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs   = length(sublist);

%%
for sub = 1:nsubs

    % Load behavioural data
    substrg  = sprintf('s%02d',sublist(sub));
    behavpath = [studypath strcat(fs,'Data/Behaviour')];
    filestr = fullfile(behavpath,  strcat(substrg,'*__catWM_EEG_postLearn.mat'));
    fstruct = dir(filestr);
    if ~isempty(fstruct)
        filename = fstruct(1).name;
        behav = load(strcat(behavpath,fs,filename));
    else, error(sprintf('%s does not exist!',substrg));
    end
    behav = behav.dat;
      
    % Load EEG data
    EEGpath   = [studypath strcat(fs,'Data/EEG/Preprocessed')];
    substrg  = sprintf('S%02d',sublist(sub));
    subpath  = sprintf('%s%s%s',EEGpath,fs,substrg);
    filestr = sprintf('%s%s%s_ftbeforeAR.mat',subpath,fs,substrg);
    if exist(filestr), load(filestr)
    else, error(sprintf('%s does not exist!',filestr))
    end
    data = data;
    
    % Load artefact data
    rejfile = sprintf('%s%s%s_artefacts.mat',subpath,fs,substrg);
    if exist(rejfile,'file')
        load(rejfile);
    else
        artfct_def_eye = [];
        artfct_def_other = [];
    end
    
    %Select variables of interest for EEG analysis
    counterbalance_cond = behav.SubjDetails.Condition;
    switch counterbalance_cond %category boundary
        case 1, bound = [170 350];
        case 2, bound = [200 20];
        case 3, bound = [230 50];
        case 4, bound = [260 80];
        case 5, bound = [290 110];
        case 6, bound = [320 140];
    end
    temp.trialCount  = behav.condition.trialCount;
    temp.blockCount  = behav.condition.blockCount;
    temp.task_type   = behav.condition.block_type;
    temp.trial_type  = behav.condition.trial_type;
    temp.target      = behav.stms.target;
    temp.RT          = behav.rsp.RT;
    temp.acc         = behav.rsp.accuracy;
    stim_to_bound = mod((temp.target-bound(1))+180, 180*2)-180; % wrap between -90 and 90, so determining distance to either boundary
    category = zeros(1,length(temp.target));
    category(stim_to_bound>0) = 1;
    category(stim_to_bound<0) = 2;
    temp.category = category;
    
    behav = temp;
    ntrial = size(behav.trialCount,2);
    artefacts = artfct_def_eye | artfct_def_other;

    %Deal with issues
    if sublist(sub) == 17 %first three trials are from a false start
        artefacts(1:3) = 1;
    end

    nartefacts(sub,1) = nnz(artefacts);
    fprintf('\nRejecting %03d trials for %s.\n\n',nartefacts(sub,1),substrg);
    good_trials = ((ntrial-nartefacts(sub,1))/ntrial).*100;
    fprintf('\n %d percent good trials remaining for %s.\n\n',round(good_trials),substrg);

    %Excluding trials from EEG data
    keep_trials = 1:size(data.trial,2);
    cfg = [];
    keep_trials = keep_trials(~artefacts);

    cfg.trials = keep_trials;
    clean_data = ft_preprocessing(cfg, data);

    %Deal with issues
    if sublist(sub) == 17 %first three trials are from a false start
        artefacts(1:3) = [];
    elseif sublist(sub) ==37
        artefacts = [artefacts(1:528); true; artefacts(529:end)]; %trial 529 missing from EEG data
    elseif sublist(sub) ==42
        artefacts = [artefacts(1:792); true(9,1); artefacts(793:end)]; %data overflow in block 22, so first 9 trials in block 23 missing
    end

    %Excluding trials from behavioural data
    new_behav = behav;
    behavfields = fieldnames(behav);
    nfields     = length(behavfields);
    for ifield = 1:nfields
        if ndims(behav.(behavfields{ifield})) == 1
            new_behav.(behavfields{ifield})(artefacts) = []; %exclude artefacts
        elseif ndims(behav.(behavfields{ifield})) == 2
            new_behav.(behavfields{ifield})(:, artefacts) = []; %exclude artefacts
        elseif ndims(behav.(behavfields{ifield})) == 3
            new_behav.(behavfields{ifield})(:, :, artefacts) = []; %excludeartefacts
        end
    end

    new_behav.ntrial = size(new_behav.trialCount,2);
    new_behav.bound  = bound;

    %Save final data
    savfile = sprintf('%s%s%s_ftEEGfinal.mat',subpath,fs,substrg);
    save(savfile,'clean_data','new_behav');

    %
    clear clean_data behav new_behav
    fprintf('\n%s saved.\n',substrg);

end