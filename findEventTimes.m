function event_times = findEventTimes(message, datastruct)
%pulls times of all events matching input message in a subject's datastruct.
message_cell = repmat({message}, length(datastruct.message.string), 1);
comp = cellfun(@strcmp, message_cell, datastruct.message.string); %, 'UniformOutput', true
times = str2double(datastruct.message.time);
event_times = times(comp);