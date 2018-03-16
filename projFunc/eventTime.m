function [evt_time] = eventTime(block, evt_Name)
%go through all trials (segments) and return the event time
%   Detailed explanation goes here

trialN = numel(block.annotation_all_trial_ids);
for t = 1:trialN
    EvtStr = block.segments{1,t}.events{1,1}.annotation_trial_event_labels;
    for e = 1:size(EvtStr,1)
        if strfind(EvtStr(e,:),evt_Name)
            EvtInd = e;
            evt_time(t,1) = block.segments{1,t}.events{1,1}.times(EvtInd);
            break
        end
    end
end


end

