
blue_05hz_optosession = find(strcmp(optoLog.laser_freq,'40 Hz') & strcmp(optoLog.laser_color,'blue') & (optoLog.laser_probe == 1) & (optoLog.extract_bin == 1));

%% Extract laminar LFP's
for session_i = 1:length(blue_05hz_optosession)
    try
        fprintf('session %i of %i \n', session_i, length(blue_05hz_optosession))
        clear data_in stim_trials

        data_in = load(fullfile(dirs.mat_data,optoLog.session{blue_05hz_optosession(session_i)}));
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

        clear baseline_lfp_activity stimulation_lfp_activity
        baseline_lfp_activity = reshape(lfp_array_aligned(:,1000+[-999:0],valid_trials_idx), 32, length([-999:0]) * size(valid_trials_idx,2));
        stimulation_lfp_activity = reshape(lfp_array_aligned(:,1000+[0:999],valid_trials_idx), 32, length([0:999]) * size(valid_trials_idx,2));

        % Parameters for pwelch
        window = 500; % Length of each segment
        noverlap = 250; % Number of overlapping samples
        nfft = 5000; % Number of FFT points

        for channel_i = 1:16
            [power_baseline(channel_i,:), f] = pwelch(baseline_lfp_activity(channel_i,:), window, noverlap, nfft, 1000, 'power');
            [power_stimulation(channel_i,:), ~] = pwelch(stimulation_lfp_activity(channel_i,:), window, noverlap, nfft, 1000, 'power');
        end

        % Find indices for the desired frequency range
        gamma_idx = find((f >= 39) & (f <= 41));
        theta_idx = find((f >= 4) & (f <= 6));

        clear gamma_power_* theta_power_*
        for channel_i = 1:16
            gamma_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, gamma_idx)));
            theta_power_baseline(channel_i,1) = sum(10*log10(power_baseline(channel_i, theta_idx)));

            gamma_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, gamma_idx)));
            theta_power_stimulation(channel_i,1) = sum(10*log10(power_stimulation(channel_i, theta_idx)));
        end




        delta_power(session_i,:) = [((gamma_power_stimulation./gamma_power_baseline)*100)-100]';

    catch
        delta_power(session_i,:) = nan(1,16);

    end



end


for session_i = 1:length(blue_05hz_optosession)

        clear halfMax index1 index2 fwhm
        % Find the half max value.
        halfMax = (min(delta_power(session_i,:)) + max(delta_power(session_i,:))) / 2;
        % Find where the data first drops below half the max.
        index1 = find(delta_power(session_i,:) >= halfMax, 1, 'first');
        % Find where the data last rises above half the max.
        index2 = find(delta_power(session_i,:) >= halfMax, 1, 'last');
        fwhm = index2-index1; % FWHM in indexes.

        try
            fwhm_out(session_i) = fwhm;
        catch
            fwhm_out(session_i) = NaN;
        end
end


figuren('Renderer', 'painters', 'Position', [100 100 600 300]);
histogram(fwhm_out,0:1:16,'LineStyle','none')
xlim([0 16]); xticks(0.5:1:16.5); xticklabels(1:16)
box off

