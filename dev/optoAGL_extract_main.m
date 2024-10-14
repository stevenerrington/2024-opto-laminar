
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

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session
exp_filename = '2021-10-18_09-57-28';
session_n = '_0010'; 
outfile_name = 'troy-optoAGL-2021-10-18';

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array
% Convert .ncs files into a binary file for use in Kilosort --------/

if ~exist([dirs.bin_data outfile_name '.dat']) == 1
    for ch_n = 1:n_channels
        clear filepart_name NCSpath spk_ncs_in

        filepart_name = ['CSC' int2str(ch_n) session_n];
        NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

        spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
        spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
    end

    % Create a binary file and export the restructure broadband data
    clear bin_out_file
    bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
    fwrite(bin_out_file,spk_ncs_out,'int16');
    fclose(bin_out_file);
    mkdir(fullfile(dirs.kilosort,outfile_name));


else
    ops = struct();
    ops.rootZ = fullfile(dirs.kilosort,outfile_name);
    ops.bin_file = [dirs.bin_data outfile_name '.dat'];
    ops.nCh = 32;
    ops.fs = 32000;

    if exist(fullfile(ops.rootZ,'params.py')) == 2
        try
            [spikes] = phy2mat(ops);
            [spk_info] = phyinfo2mat(ops);
            fprintf(['- phy import successful!: ' outfile_name ' \n'])
        catch
            fprintf(['- error importing phy curation: ' outfile_name ' \n'])
        end
    else
        fprintf(['- no kilosort output detected: ' outfile_name ' \n'])
        spikes = [];
        spk_info = [];
    end

end

% Local field potential data -------------------------------------------------------
filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n '.ncs'],32);
lfp = ft_read_neuralynx_interp(filelabels_lfp);
lfp = lfp.trial{1};

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
clear event_table_raw event_table
event_table_raw = get_event_table(ops);
ops.event_port = 2;
event_table = get_agl_t_trials_nlx(event_table_raw,ops);

ops.aligntime = event_table.laserOnset_ms;

[lfp_aligned, lfp_array_aligned] = get_lfp_aligned(lfp,ops.aligntime,ops);



stim_trials = find(~isnan(ops.aligntime));


stim_average = nanmean(lfp_array_aligned(:,:,stim_trials),3);


figuren('Renderer', 'painters', 'Position', [100 100 1200 400]);

channels = 1:16;
subplot(1,1,1); hold on
for ch_i = 1:length(channels)
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  stim_average(ch_i,:)+10*(ch_i-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-50 (length(channels)*10)+50]); yticks([10*([1:32]-1)]); yticklabels(num2cell([1:32]))
xlim([-250 1500])




%% Spikes
ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(spikes,ops.aligntime,ops);

spike_labels = fieldnames(sdf);

for spike_i = 1:length(spike_labels)
    figuren;
    subplot(1,2,1); hold on
    plot(ops.timewin,smooth(nanmean(sdf.(spike_labels{spike_i})(isnan(event_table.laser_gamma_ms),:)),100))
    plot(ops.timewin,smooth(nanmean(sdf.(spike_labels{spike_i})(~isnan(event_table.laser_gamma_ms),:)),100))
    title('gamma')

    subplot(1,2,2); hold on
    plot(ops.timewin,smooth(nanmean(sdf.(spike_labels{spike_i})(isnan(event_table.laser_theta_ms),:)),100))
    plot(ops.timewin,smooth(nanmean(sdf.(spike_labels{spike_i})(~isnan(event_table.laser_theta_ms),:)),100))
    title('theta')
end