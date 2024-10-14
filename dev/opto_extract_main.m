
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
session_i = 175;
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
opto_event = get_opto_trials(event_table_raw);

% %% Output data
% save(fullfile(dirs.mat_data,[outfile_name '.mat']),'event_table','spikes','spk_info','lfp','-v7.3')
% fprintf('Extracted data successfully saved to %s    \n', fullfile(dirs.mat_data,[outfile_name '.mat']))
% fprintf(' - Events  ✓   \n')
% fprintf(' - Spikes  ✓   \n')
% fprintf(' - LFP     ✓   \n')