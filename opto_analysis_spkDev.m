
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


session_name = {'troy-opto-2021-11-22a', 'troy-opto-2021-11-22b', 'troy-opto-2021-11-22c'};
cond_name = {'1_Blue 05 Hz', '2_Blue 40 Hz', '3_Red 40 Hz'};


for session_loop_i = 1:4
    session_idx = find(strcmp(optoLog.session,session_name{session_loop_i}));
    session_i = session_idx;
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

    % Spiking data
    ops = struct();
    ops.rootZ = fullfile(dirs.kilosort,outfile_name);
    ops.bin_file = [dirs.bin_data outfile_name '.dat'];
    ops.nCh = 32;
    ops.fs = 32000;

    [spikes] = phy2mat(ops);
    %[spk_info] = phyinfo2mat(ops);

    ops.aligntime = aligntime;
    ops.timewin = -1000:5000;
    ops.sdf_filter = 'PSP';
    [sdf{session_loop_i}, raster{session_loop_i}] = get_spikes_aligned(spikes,aligntime,ops);

    waveform_out{session_loop_i} = spikes.waveform;

end

%% Spike density function

unit = 'DSP26b';
xlim_val = [-500 1500];
clear opto_singleunit_fig

cond = 1;

opto_singleunit_fig(1,1)=gramm('x',[raster{cond}.(unit)],...
    'color',[repmat(cond_name(1),length(raster{cond}.(unit)),1)]);
opto_singleunit_fig(1,1).geom_raster('geom',{'line'});
opto_singleunit_fig(1,1).axe_property('XLim',xlim_val);

opto_singleunit_fig(2,1)=gramm('x',ops.timewin,'y',[sdf{cond}.(unit)],...
    'color',[repmat(cond_name(1),size(sdf{cond}.(unit),1),1)]);
opto_singleunit_fig(2,1).stat_summary();
opto_singleunit_fig(2,1).axe_property('XLim',xlim_val,'YLim',[0 40]);



opto_singleunit_fig(1,1).set_layout_options...
    ('Position',[0.1 0.6 0.35 0.25],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

opto_singleunit_fig(2,1).set_layout_options...
    ('Position',[0.1 0.1 0.35 0.4],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
opto_singleunit_fig.draw();




%% Waveform



opto_waveform_fig(1,1)=gramm('x',1:82,'y',spikes.waveform.WAV15a);
opto_waveform_fig(1,1).stat_summary();
opto_waveform_fig(1,1).axe_property('XLim',[20 60],'YLim',[-15 20]);

figure('Renderer', 'painters', 'Position', [100 100 200 200]);
opto_waveform_fig.draw();



%% Autocorr
clear sdf_trial acorr

[acorr_blue_40Hz, ~] = xcorr(nanmean(sdf{1}.(unit)(:,1000+[0:1000])), 'coeff');
[acorr_red_40Hz, ~] = xcorr(nanmean(sdf{2}.(unit)(:,1000+[0:1000])), 'coeff');
[acorr_blue_05Hz, ~] = xcorr(nanmean(sdf{4}.(unit)(:,1000+[0:1000])), 'coeff');
[acorr_red_05Hz, lags] = xcorr(nanmean(sdf{3}.(unit)(:,1000+[0:1000])), 'coeff');






figuren('Renderer', 'painters', 'Position', [100 100 300 500]);
subplot(2,1,1); hold on
plot(lags,acorr_blue_05Hz,'b')
plot(lags,acorr_red_05Hz,'r')
xlim([-500 500]); ylim([0 1])

subplot(2,1,2); hold on
plot(lags,acorr_blue_40Hz,'b')
plot(lags,acorr_red_40Hz,'r')
xlim([-500 500]); ylim([0 1])

