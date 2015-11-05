function outputdata = getalleyedata(datadir)
    %Runs importeyedata on all *_eye.tsv in specified folder.
    currdir = pwd;
    %Include '/' at the end of the path!
    if datadir(length(datadir)) ~= '/'
        datadir = [datadir '/'];
    end
    % pull all data filenames and create empty output struct
    filenames = dir([datadir '*_eye.tsv']);
    fnames = {filenames.name}';
    %output = cellfun(@importeyedata, fnames);
    output = struct('subject', cell(1, length(fnames)), 'data', cell(1, length(fnames)), 'message', cell(1, length(fnames)));
    cd(datadir)
    %run importeyedata.m on each file and separate data and message
    for i = 1:length(fnames)
        try
            raw = importeyedata(fnames{i});
            output(i).subject = fnames{i}; %re_matches{1};
            output(i).data = raw.data;
            output(i).message = raw.message;
        catch 
            disp(['An error occurred while importing data from ' fnames{i}]);
            disp('Skipping this file.')
            i = i+1;
        end
    end
    outputdata = output;
    cd(currdir)
end
