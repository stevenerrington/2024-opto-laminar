%%
%{
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab opto script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%}

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
optoLog = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'opto'));

session_list = [316, 318, 319, 320];
neuron_frontal_label_list = {'DSP26b','DSP26b','DSP26b','DSP26b'};
neuron_auditory_label_list = {'DSP15c','DSP15b','DSP15b','DSP15b'};

for session_list_i = 1:length(session_list)

    session_i = session_list(session_list_i);

    monkey = optoLog.monkey{session_i}; % Monkey name [troy, chief]

    % Experimental parameters -------------------------------------------
    n_channels = 32; % Number of channels recorded in session

    % Key setup variables
    exp_filename = optoLog.data_folder{session_i}; % Experimental raw data
    task = optoLog.task{session_i}; % Experiment type [agl, opto]
    session_n = optoLog.file_n{session_i}; % Experimental file tag

    % Define experimental/data directories -------------------------------
    outfile_name = optoLog.session{session_i}; % Processed file name

    dirs.raw_data = optoLog.data_dir{session_i};

    % Behavioral data -------------------------------------------------------
    % Read in events
    ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
    clear event_table_raw event_table
    event_table_raw = get_event_table(ops);
    opto_event = get_opto_trials(event_table_raw);
    aligntime = opto_event.laserOnset_ms;


    % Local field potential data -------------------------------------------------------
    filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n '.ncs'],32);
    lfp = ft_read_neuralynx_interp(filelabels_lfp);
    lfp = lfp.trial{1};

    ops.timewin = [-1000:5000];
    ops.freq = [1 60];
    ops.ch_extract = [1:32];
    lfp = patch_fault_ch(lfp,23);
    [~, lfp_array{session_list_i,1}] = get_lfp_aligned(lfp,aligntime,ops);

end

%%
count = 0;
for session_list_i = 1:4

    nan_trials = []; nan_trials_idx = []; valid_trials_idx = [];
    nan_trials = isnan(lfp_array{session_list_i,1});
    nan_trials_idx = squeeze(nan_trials(1,:,:));
    valid_trials_idx = find(nan_trials_idx(1,:) == 0);

    baseline_lfp_activity{session_list_i,1} = reshape(lfp_array{session_list_i,1}(:,1000+[-999:0],valid_trials_idx), n_channels, length([-999:0]) * size(valid_trials_idx,2));
    stimulation_lfp_activity{session_list_i,1} = reshape(lfp_array{session_list_i,1}(:,1000+[0:999],valid_trials_idx), n_channels, length([0:999]) * size(valid_trials_idx,2));

    % Parameters for pwelch
    window = 500; % Length of each segment
    noverlap = 250; % Number of overlapping samples
    nfft = 5000; % Number of FFT points

    for channel_i = 1:32
        [power_baseline{session_list_i,1}(channel_i,:), f] = pwelch(baseline_lfp_activity{session_list_i,1}(channel_i,:), window, noverlap, nfft, 1000, 'power');
        [power_stimulation{session_list_i,1}(channel_i,:), ~] = pwelch(stimulation_lfp_activity{session_list_i,1}(channel_i,:), window, noverlap, nfft, 1000, 'power');
    end

    % Find indices for the desired frequency range
    gamma_idx = find((f >= 39) & (f <= 41));
    theta_idx = find((f >= 3) & (f <= 9));

    for channel_i = 1:32
        count = count + 1;
        gamma_power_baseline(count,1) = sum(10*log10(power_baseline{session_list_i,1}(channel_i, gamma_idx)));
        theta_power_baseline(count,1) = sum(10*log10(power_baseline{session_list_i,1}(channel_i, theta_idx)));

        gamma_power_stimulation(count,1) = sum(10*log10(power_stimulation{session_list_i,1}(channel_i, gamma_idx)));
        theta_power_stimulation(count,1) = sum(10*log10(power_stimulation{session_list_i,1}(channel_i, theta_idx)));
    end

end

cond_label = [];

for cond_i = 1:4
    cond_label = [cond_label; repmat({[optoLog.laser_color{session_list(cond_i)} '-' optoLog.laser_freq{session_list(cond_i)}]},32,1)];
end

cond_label = repmat(cond_label,4,1);
area_label = repmat([repmat({'auditory'},16,1); repmat({'frontal'},16,1)],16,1);

epoch_label = [repmat({'baseline'}, 256,1); repmat({'stimulation'}, 256,1)];


power_data = [gamma_power_baseline; theta_power_baseline; gamma_power_stimulation; theta_power_stimulation];


clear power_bar_fig

power_bar_fig(1,1) = gramm('x',cond_label,'y',power_data,'color', epoch_label, 'column',area_label);
power_bar_fig(1,1).stat_summary('geom',{'bar','black_errorbar'},'type','sem');
power_bar_fig(1,1).axe_property('YLim',[180 500])
figure('Position',[100 100 800 550]);
power_bar_fig.draw();


electrode_contact = repmat([1:32]',16,1);

clear power_depth_fig
power_depth_fig(1,1) = gramm('x',electrode_contact,'y',power_data, 'color', cond_label,...
    'subset', strcmp(area_label,'auditory') & (strcmp(cond_label,'blue-40 Hz') | strcmp(cond_label,'red-40 Hz')),...
    'column', epoch_label);
power_depth_fig(1,1).stat_summary('geom',{'line','point'});
power_depth_fig(1,1).axe_property('YLim',[320 500])
figure('Position',[100 100 800 550]);
power_depth_fig.draw();



clear power_depth_fig
power_depth_fig(1,1) = gramm('x',electrode_contact,'y',power_data, 'color', epoch_label,...
    'subset', strcmp(area_label,'auditory') & (strcmp(cond_label,'blue-40 Hz')));
power_depth_fig(1,1).stat_summary('geom',{'line','point'});
power_depth_fig(1,1).axe_property('YLim',[320 420])
figure('Position',[100 100 800 550]);
power_depth_fig.draw();

