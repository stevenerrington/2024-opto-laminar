
for trial_i = 1 : size(opto_event,1)
    test(:,trial_i) = opto_event.laser_on{trial_i}-opto_event.laserOnset_ms(trial_i);
    test2(:,trial_i) = opto_event.laser_off{trial_i}-opto_event.laserOnset_ms(trial_i);
end


%% 

aligntime = opto_event.laserOnset_ms;

ops.timewin = [-1000:5000];
ops.freq = [1 60];
ops.ch_extract = [1:32];
lfp = patch_fault_ch(lfp,23);
[~, lfp_array] = get_lfp_aligned(lfp,aligntime,ops);
trial_average_lfp = nanmean(lfp_array,3);

figuren('Renderer', 'painters', 'Position', [100 100 400 1200]);

subplot(1,1,1); hold on
for ch_i = 1:n_channels
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+10*(ch_i-1),'color',color_line)
end

% vline(test,'k-')
% vline(test2,'r-')
set(gca,'ydir', 'reverse')
ylim([-15 (length(1:n_channels)*10)]); yticks([10*([1:32]-1)]); yticklabels(num2cell([1:32]))
xlim([-250 1000])



ops.aligntime = aligntime;
ops.plot_ch = fieldnames(spikes.time);
single_unit_rastersdf_figure(spikes,ops)

%%




baseline_lfp_activity = reshape(lfp_array(:,1000+[-999:0],:), n_channels, length([-999:0]) * size(lfp_array,3));
stimulation_lfp_activity = reshape(lfp_array(:,1000+[0:999],:), n_channels, length([-999:0]) * size(lfp_array,3));


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
theta_idx = find((f >= 3) & (f <= 9));

for channel_i = 1:32
    gamma_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, gamma_idx)));
    theta_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, theta_idx)));

    gamma_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, gamma_idx)));
    theta_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, theta_idx)));
end


power_data = [gamma_power_baseline; theta_power_baseline; gamma_power_stimulation; theta_power_stimulation];
epoch_label = [repmat({'baseline'}, 64, 1); repmat({'stimulation'}, 64, 1)];
freq_label = [repmat({'gamma'}, 32, 1); repmat({'theta'}, 32, 1); repmat({'gamma'}, 32, 1); repmat({'theta'}, 32, 1)];

power_bar_fig(1,1) = gramm('x',freq_label,'y',power_data,'color', epoch_label, 'subset', strcmp(freq_label, 'gamma'));
power_bar_fig(1,1).stat_summary('geom',{'bar','black_errorbar'});
figure('Position',[100 100 800 550]);
power_bar_fig.draw();


