function times = findEventTimes(event, datastruct)
%pulls times of all events matching ONE input event string in a subject's 'message' struct.
ind = strcmp(datastruct.string, event);
times = datastruct.time(ind);
