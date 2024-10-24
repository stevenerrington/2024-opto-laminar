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

%% LFP Conditions
clear session_idx
session_idx.blue_05hz = 184;
session_idx.blue_40hz = 185;
session_idx.red_05hz = 187;
session_idx.red_40hz = 186;
cond_labels = fieldnames(session_idx);
ops.timewin = [-1000:5000];

%% ERSP for each condition
for cond_i = 1:4
    clear data_in stim_trials

    data_in = load(fullfile(dirs.mat_data,optoLog.session{session_idx.(cond_labels{cond_i})}));
    ops.aligntime = data_in.opto_event.laserOnset_ms;
    stim_trials = find(~isnan(ops.aligntime));

    clear lfp* data  nan_trials nan_trials_idx valid_trials_idx

    % Run alignment algorithms
    ops.timewin = -1000:5000;
    ops.freq = [2 200];
    [~, signal_out] = get_lfp_aligned(data_in.lfp,ops.aligntime,ops);


    nan_trials = []; nan_trials_idx = []; valid_trials_idx = [];
    nan_trials = isnan(signal_out);
    nan_trials_idx = squeeze(nan_trials(1,:,:));
    valid_trials_idx = find(nan_trials_idx(1,:) == 0);

    data = signal_out(1:16,:,valid_trials_idx);

    % EEGlab analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters and configuration -----------------------------------------
    srate =  1000; % Sampling rate (Hz)
    epochmin = ops.timewin(1)/1000; % Epoch start (sec)
    epochmax = ops.timewin(end)/1000; % Epoch end (sec)
    basemin = -500; % Baseline window start (ms)
    basemax = 0; % Baseline window end (ms)

    clear nTr in
    nTr=size(data,3); % Number of trials
    in(:,1)=[1:nTr]'; % Dummy variable for trial n
    in(:,2)=abs(epochmin)*ones(nTr,1); % Dummy variable for trial n

    freq_range=[2.5 100]; % Frequency range for ERSP analysis
    maxfreq = max(freq_range); % Max frequency
    padratio = 2; % Pad ratio
    alpha_val = 0.01; % Alpha
    maxersp = 6;


    % Setup data ------------------------------------------------------------
    clear EEG

    % - Structure data for EEGlab -------------------------------------------
    % (adapted from YK earlier matlab code)
    EEG = pop_importdata('dataformat', 'array', 'data', 'data', 'srate',srate, 'nbchan',32);
    EEG = eeg_checkset(EEG);

    % - Define epochs (although data is already aligned)
    EEG = pop_importepoch(EEG, in, { 'Epoch', 'stim'}, 'latencyfields',{ 'stim'}, 'timeunit',1, 'headerlines',0);
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'stim'  }, [epochmin         epochmax], 'newname', 'Level epochs', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );

    % - Perform baseline correction
    EEG = pop_rmbase( EEG, [basemin    0]);
    EEG = eeg_checkset( EEG );

    chan_i = 10;

    clear ersp itc powbase times freqs erspboot itcboot alltfX

    fprintf('Running analysis on channel %i \n', chan_i)
    [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
        1, chan_i, [EEG.xmin EEG.xmax]*srate, [3 0.7], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
        'plotphase', 'off', 'alpha', alpha_val, 'naccu', 200, 'baseboot',1,'rmerp','off', ...
        'erspmax', maxersp, 'plotersp','off', 'plotitc','off','baseline',[basemin basemax],'marktimes',0);

    clear psd_freqs psd_values
    [psd_values, psd_freqs] = pop_spectopo(EEG, 1, [0 1] * srate, 'EEG', ...
                                      'freqrange', freq_range, 'plotchan', chan_i, 'electrodes', 'off','plot','off');

    ersp_out{cond_i} = ersp;
    times_out{cond_i} = times;
    freqs_out{cond_i} = freqs;
    psd_out{cond_i} = psd_values;
    psd_freqs_out{cond_i} = psd_freqs;
end

% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorscale_blue = abs(cbrewer('seq','Blues',500));
colorscale_red = abs(cbrewer('seq','Reds',500));

figuren('Renderer', 'painters', 'Position', [100 100 900 400]);
ax1 = subplot(2,2,1);
imagesc('XData',times_out{1},'YData',freqs_out{1},'CData',ersp_out{1})
xlim([-200 1500]); ylim([min(freqs_out{1}) 50])
vline([0 1000], 'k-');
colorbar; colormap(ax1, colorscale_blue); clim([0 5])

ax2 = subplot(2,2,3);
imagesc('XData',times_out{3},'YData',freqs_out{3},'CData',ersp_out{3})
xlim([-200 1500]); ylim([min(freqs_out{3}) 50])
vline([0 1000], 'k-');
colorbar; colormap(ax2, colorscale_red); clim([0 5])

ax3 = subplot(2,2,2);
imagesc('XData',times_out{2},'YData',freqs_out{2},'CData',ersp_out{2})
xlim([-200 1500]); ylim([min(freqs_out{2}) 50])
vline([0 1000], 'k-');
colorbar; colormap(ax3, colorscale_blue); clim([0 10])

ax4 = subplot(2,2,4);
imagesc('XData',times_out{4},'YData',freqs_out{4},'CData',ersp_out{4})
xlim([-200 1500]); ylim([min(freqs_out{4}) 50])
vline([0 1000], 'k-');
colorbar; colormap(ax4, colorscale_red); clim([0 10])


%% Extract example LFP's for each condition
figuren('Renderer', 'painters', 'Position', [100 100 900 100]);

for cond_i = 1:length(cond_labels)
    clear data_in stim_trials
    data_in = load(fullfile(dirs.mat_data,optoLog.session{session_idx.(cond_labels{cond_i})}));
    ops.aligntime = data_in.opto_event.laserOnset_ms;
    stim_trials = find(~isnan(ops.aligntime));

    % Behavioral data -------------------------------------------------------
    % Read in events
    clear lfp*
    [lfp_aligned, lfp_array_aligned] = get_lfp_aligned(data_in.lfp,ops.aligntime,ops);

    % Set line properties
    if cond_i < 3
        line_color = 'blue';
    else
        line_color = 'red';
    end

    if cond_i == 1 | cond_i == 3
        subplot_i = 1;
        ylim_i = [-80 80];
    else 
        subplot_i = 2;
        ylim_i = [-40 40];
    end

    subplot(1,2,subplot_i); hold on
    plot(ops.timewin, nanmean(lfp_array_aligned(10,:,stim_trials),3),'LineWidth',1,'Color',line_color)
    xlim([-250 1500]); ylim(ylim_i); yticks([-80:40:80]); box off
    vline(0,'k')

    if cond_i == 1 | cond_i == 2
        set(gca,'XColor', 'None')
    end

end

%% 
figuren('Renderer', 'painters', 'Position', [100 100 600 450]); hold on
subplot(2,2,[1 3]); hold on
plot(psd_freqs_out{1},psd_out{1}(10,:),'color',colorscale_blue(size(colorscale_blue,1),:),'LineWidth',1)
plot(psd_freqs_out{2},psd_out{2}(10,:),'color',colorscale_blue(size(colorscale_blue,1)/2,:),'LineWidth',1)
plot(psd_freqs_out{1},psd_out{3}(10,:),'color',colorscale_red(size(colorscale_red,1),:),'LineWidth',1)
plot(psd_freqs_out{2},psd_out{4}(10,:),'color',colorscale_red(size(colorscale_red,1)/2,:),'LineWidth',1)
xlim([0 45])

subplot(2,2,[2]); hold on
plot(psd_freqs_out{1},psd_out{1}(10,:),'color',colorscale_blue(size(colorscale_blue,1),:),'LineWidth',1)
plot(psd_freqs_out{2},psd_out{2}(10,:),'color',colorscale_blue(size(colorscale_blue,1)/2,:),'LineWidth',1)
plot(psd_freqs_out{1},psd_out{3}(10,:),'color',colorscale_red(size(colorscale_red,1),:),'LineWidth',1)
plot(psd_freqs_out{2},psd_out{4}(10,:),'color',colorscale_red(size(colorscale_red,1)/2,:),'LineWidth',1)
xlim([0 10]); ylim([20 35])

%% Extract laminar LFP's
clear data_in stim_trials
cond_i = 4;
data_in = load(fullfile(dirs.mat_data,optoLog.session{session_idx.(cond_labels{cond_i})}));
ops.aligntime = data_in.opto_event.laserOnset_ms;
stim_trials = find(~isnan(ops.aligntime));

% Behavioral data -------------------------------------------------------
% Read in events
clear lfp*
[lfp_aligned, lfp_array_aligned] = get_lfp_aligned(data_in.lfp,ops.aligntime,ops);

nan_trials = []; nan_trials_idx = []; valid_trials_idx = [];
nan_trials = isnan(lfp_array_aligned);
nan_trials_idx = squeeze(nan_trials(1,:,:));
valid_trials_idx = find(nan_trials_idx(1,:) == 0);

baseline_lfp_activity = reshape(lfp_array_aligned(:,1000+[-999:0],valid_trials_idx), 32, length([-999:0]) * size(valid_trials_idx,2));
stimulation_lfp_activity = reshape(lfp_array_aligned(:,1000+[0:999],valid_trials_idx), 32, length([0:999]) * size(valid_trials_idx,2));

% Parameters for pwelch
window = 500; % Length of each segment
noverlap = 250; % Number of overlapping samples
nfft = 5000; % Number of FFT points

for channel_i = 1:32
    [power_baseline(channel_i,:), f] = pwelch(baseline_lfp_activity(channel_i,:), window, noverlap, nfft, 1000, 'power');
    [power_stimulation(channel_i,:), ~] = pwelch(stimulation_lfp_activity(channel_i,:), window, noverlap, nfft, 1000, 'power');
end

% Find indices for the desired frequency range
gamma_idx = find((f >= 39) & (f <= 41));
theta_idx = find((f >= 4) & (f <= 6));

for channel_i = 1:32
    gamma_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, gamma_idx)));
    theta_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, theta_idx)));

    gamma_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, gamma_idx)));
    theta_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, theta_idx)));
end


% Plot LFP and power x channel
figuren('Renderer', 'painters', 'Position', [100 100 600 400]);

subplot(1,3,[1 2]); hold on
trial_average_lfp = nanmean(lfp_array_aligned,3);

for ch_i = 1:16
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+30*(ch_i-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-30 (16*30)]); yticks([30*([1:32]-1)]); yticklabels(num2cell([1:16]))
xlim([-500 1000]); vline(0,'k')
set(gca,'Ycolor', 'None')

subplot(1,3,3); hold on
delta_power = ((theta_power_stimulation(1:16)./theta_power_baseline(1:16))*100)-100;
plot(delta_power,1:16,'Color',color_line)
set(gca,'ydir', 'reverse')
set(gca,'Ycolor', 'None')
xlim([-5 30]); ylim([0 17]); yticks([-5:5:30])


%% Single unit
cond_i = 1;

outfile_name = optoLog.session{session_idx.(cond_labels{cond_i})};
data_in = load(fullfile(dirs.mat_data,optoLog.session{session_idx.(cond_labels{cond_i})}));

ops = struct();
ops.rootZ = fullfile(dirs.kilosort,outfile_name);
ops.bin_file = [dirs.bin_data outfile_name '.dat'];
ops.nCh = 32;
ops.fs = 32000;

[spikes] = phy2mat(ops);
fprintf(['- phy import successful!: ' outfile_name ' \n']);


ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';
ops.aligntime = data_in.opto_event.laserOnset_ms;

[sdf, raster] = get_spikes_aligned(spikes,ops.aligntime,ops);

single_unit_rastersdf_figure(spikes,ops)