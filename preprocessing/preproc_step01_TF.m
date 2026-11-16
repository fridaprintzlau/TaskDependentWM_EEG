%% EEG preprocessing - step 1: filter and epoch for time-frequency analysis
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

clear all
close all
clc

studypath = ''; %define path where the data folder sits
fieldtrippath = ''; %add fieldtrip path here
addpath(fieldtrippath);
ft_defaults;

fs = filesep;
input_dir = [studypath strcat(fs,'Data/EEG/Raw')];  % input directory 
output_dir = [studypath strcat(fs,'Data/EEG/Preprocessed')]; %output directory

sublist = [02,06,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,36,37,38,39,40,41,42];
nsubs   = length(sublist);

%%
for isub = 1:nsubs
    %% Subject details

    sub = sublist(isub);
    substrg = sprintf('%02d%s',sub);
    subDir  = sprintf('%s%sS%s',output_dir,fs,substrg);
    if ~exist(subDir); mkdir(subDir); end
    
    data_name = strcat(input_dir, sprintf('/S%s%ss%s_WM.eeg', substrg, fs, substrg));% define the data path and its name
    
    cfg = [];
    cfg.dataset = data_name;
    
    %% Read in raw data

    cfg                    = [];
    cfg.dataset            = data_name;
    cfg.channel            = 'all';             % define channel type
    data                   = ft_preprocessing(cfg); % read raw data

    %% HP and LP filtering

    cfg                = [];
    cfg.bpfilter       = 'yes';    % enable band-pass filtering
    cfg.bpfreq         = [0.1 45]; % set up the frequency for high-pass and low-pass filter
    cfg.dftfilter      = 'no';     % enable notch filtering to eliminate power line noise
    data_bp            = ft_preprocessing(cfg,data);
    
    %% Epoch data

    % Read events
    cfg                     = [];
    cfg.trialdef.prestim    = 1;                   % in seconds
    cfg.trialdef.poststim   = 4;                   % in seconds
    cfg.trialdef.eventtype  = 'Stimulus';            % get a list of the available types
    cfg.trialdef.eventvalue = {'S 12', 'S 22'};      % trigger values of interest (onset of WM stimulus)
    cfg.dataset             = data_name;             % set the name of the dataset
    cfg_tr_def              = ft_definetrial(cfg);   % read the list of the specific stimulus

    data = ft_redefinetrial(cfg_tr_def, data_bp);

    % Append data in case of multiple sessions
    if sub == 31 || sub ==33 || sub==41 || sub==42
        data_old = data;
        data_name_new = strcat(input_dir, sprintf('/S%s%ss%sb_WM.eeg', substrg, fs, substrg));% define the data path and its name
    

        cfg                     = [];
        cfg.trialdef.prestim    = 0.5;                   % in seconds
        cfg.trialdef.poststim   = 3.2;                   % in seconds
        cfg.trialdef.eventtype  = 'Stimulus';            % get a list of the available types
        cfg.trialdef.eventvalue = {'S 12','S 22'};       % trigger values of interest
        cfg.dataset             = data_name_new;         % set the name of the dataset
        cfg_tr_def              = ft_definetrial(cfg);   % read the list of the specific stimulus
        % Read raw data
        cfg                    = [];
        cfg.dataset            = data_name_new;
        cfg.channel            = {'all'};
        data_temp              = ft_preprocessing(cfg); % read raw data
        %Epoch data
        data_new = ft_redefinetrial(cfg_tr_def, data_temp);

        cfg = [];
        cfg.keepsampleinfo = 'no';
        data_append = ft_appenddata(cfg, data_old, data_new);
        data = data_append;
    end

    %Save data before artefact rejection
    save(strcat(subDir, fs, sprintf('S%s_ftbeforeAR_TF',substrg)), 'data')
    
end