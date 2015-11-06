function eye_data = getalleyedata(datadir)
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
            disp(['Importing ' fnames{i}])
            raw = importeyedata(fnames{i});
            output(i).subject = fnames{i};
            output(i).data = raw.data;
            output(i).message = raw.message;
        catch 
            disp(['An error occurred while importing data from ' fnames{i}]);
            output(i).subject = fnames{i};
            disp('Skipping this file.')
            continue
        end
    end
    disp(['Saving data as eye_data.m in ' datadir])
    eye_data = output;
    save('eye_data', 'output');
    cd(currdir)
end
