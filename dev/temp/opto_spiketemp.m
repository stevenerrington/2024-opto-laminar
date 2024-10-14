
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


    % % Local field potential data -------------------------------------------------------
    % filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n '.ncs'],32);
    % lfp = ft_read_neuralynx_interp(filelabels_lfp);
    % lfp = lfp.trial{1};
    %
    % ops.timewin = [-1000:5000];
    % ops.freq = [1 60];
    % ops.ch_extract = [1:32];
    % lfp = patch_fault_ch(lfp,23);
    % [~, lfp_array] = get_lfp_aligned(lfp,aligntime,ops);
    % trial_average_lfp = nanmean(lfp_array,3);


    % Spiking data
    ops = struct();
    ops.rootZ = fullfile(dirs.kilosort,outfile_name);
    ops.bin_file = [dirs.bin_data outfile_name '.dat'];
    ops.nCh = 32;
    ops.fs = 32000;

    [spikes] = phy2mat(ops);
    [spk_info] = phyinfo2mat(ops);

    ops.aligntime = aligntime;
    ops.timewin = -1000:5000;
    ops.sdf_filter = 'Gauss';
    [sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

    frontal_sdf{session_list_i,1} = sdf.(neuron_frontal_label_list{session_list_i});
    auditory_sdf{session_list_i,1} = sdf.(neuron_auditory_label_list{session_list_i});

    frontal_raster{session_list_i,1} = raster.(neuron_frontal_label_list{session_list_i});
    auditory_raster{session_list_i,1} = raster.(neuron_auditory_label_list{session_list_i});

end


frontal_sdf_plotdata = []; auditory_sdf_plotdata = [];
frontal_raster_plotdata = []; auditory_raster_plotdata = [];
cond_label = [];

for cond_i = 1:4
    frontal_sdf_plotdata = [ frontal_sdf_plotdata ; frontal_sdf{cond_i,1} ];
    auditory_sdf_plotdata = [ auditory_sdf_plotdata ; auditory_sdf{cond_i,1} ];

    frontal_raster_plotdata = [ frontal_raster_plotdata ; frontal_raster{cond_i,1} ];
    auditory_raster_plotdata = [ auditory_raster_plotdata ; auditory_raster{cond_i,1} ];

    cond_label = [cond_label; repmat({[optoLog.laser_color{session_list(cond_i)} '-' optoLog.laser_freq{session_list(cond_i)}]},75,1)];
end

%%

clear opto_singleunit_fig

opto_singleunit_fig(1,1)=gramm('x',auditory_raster_plotdata,'color',cond_label);
opto_singleunit_fig(1,1).geom_raster('geom',{'line'});
opto_singleunit_fig(1,1).axe_property('XLim',[-100 1000]);

opto_singleunit_fig(2,1)=gramm('x',ops.timewin,'y',auditory_sdf_plotdata,'color',cond_label);
opto_singleunit_fig(2,1).stat_summary();
opto_singleunit_fig(2,1).axe_property('XLim',[-100 1000],'YLim',[0 40]);

opto_singleunit_fig(1,2)=gramm('x',frontal_raster_plotdata,'color',cond_label);
opto_singleunit_fig(1,2).geom_raster('geom',{'line'});
opto_singleunit_fig(1,2).axe_property('XLim',[-100 1000]);

opto_singleunit_fig(2,2)=gramm('x',ops.timewin,'y',frontal_sdf_plotdata,'color',cond_label);
opto_singleunit_fig(2,2).stat_summary();
opto_singleunit_fig(2,2).axe_property('XLim',[-100 1000],'YLim',[0 40]);


opto_singleunit_fig(1,1).set_layout_options...
    ('Position',[0.1 0.8 0.35 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

opto_singleunit_fig(2,1).set_layout_options...
    ('Position',[0.1 0.1 0.35 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

opto_singleunit_fig(1,2).set_layout_options...
    ('Position',[0.5 0.8 0.35 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

opto_singleunit_fig(2,2).set_layout_options...
    ('Position',[0.5 0.1 0.35 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
opto_singleunit_fig.draw();