function times = findEventTimes(event, datastruct)
%pulls times of all events matching ONE input event string in a subject's
%loaded _eye.mat file
ind = strcmp(datastruct.message.string, event);
times = datastruct.message.time(ind);
times = cellfun(@str2num, times);