function opto_event = get_opto_trials(event_table)

clear event_code

event_code.trialStart = 129;
event_code.laser_seq_start = 137;
event_code.laser_on = 1;
event_code.laser_off = 0;

code = event_table.value;

trial_start_codes = find(code == event_code.trialStart);


% Trial information
clear trial_n cond_value cond_label
for i = 1:length(trial_start_codes)

    % Get trial number
    trial_n(i,1) = i;
   
end


event_labels = fieldnames(event_code);
event_labels = event_labels(1:4);

% Timestamps
clear trial_codeblock event_time_ms
for i = 1:length(trial_start_codes)

    if i < length(trial_start_codes)
        trial_codeblock{i,1} = code(trial_start_codes(i):trial_start_codes(i+1)-1);
    else
        trial_codeblock{i,1} = code(trial_start_codes(i):end);
    end


    for event_i = 1:2
        try
            event_time_ms.(event_labels{event_i})(i,1) =...
                event_table.timestamp_ms...
                (trial_start_codes(i)+find(trial_codeblock{i,1}...
                == event_code.(event_labels{event_i}))-1);
        catch
            event_time_ms.(event_labels{event_i})(i,1) = NaN;
        end

    end

    for event_i = 3:4
        try
            event_time_ms.(event_labels{event_i}){i,1} =...
                event_table.timestamp_ms...
                (trial_start_codes(i)+find(trial_codeblock{i,1}...
                == event_code.(event_labels{event_i}))-1);
        catch
            event_time_ms.(event_labels{event_i}){i,1} = NaN;
        end

    end
end



opto_event = table(trial_n, trial_codeblock,...
    event_time_ms.trialStart,event_time_ms.laser_seq_start,...
    'VariableNames',{'trial_n','trial_codeblock','trialStart_ms','laserOnset_ms'});

