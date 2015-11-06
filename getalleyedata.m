function eye_data = getalleyedata(datadir)
    %Saves eye data as .mat for all *_eye.tsv files in indicated dir
    currdir = pwd;
    %add '/' at the end of the path if missing
    if datadir(length(datadir)) ~= '/'
        datadir = [datadir '/'];
    end
    % pull all data filenames
    filenames = dir([datadir '*_eye.tsv']);
    fnames = {filenames.name}';
    %choose the pattern for saving files
    %chop off '.tsv' at the end of the data filename
    savenames = cellfun(@(x) x(1:end-4), fnames, 'UniformOutput', false);
    %run importeyedata.m on each file and separate data and message
    cd(datadir)
    for i = 1:length(fnames)
        try
            disp(['Importing ' fnames{i}])
            raw = importeyedata([datadir fnames{i}]);
            eye_data = raw.data;
            message = raw.message;
            %save as eye_data and message to subject's filename as .mat
            disp(['Saving data as ' savenames{i} '.mat'])
            save([savenames{i}], 'eye_data', 'message');
            clear eye_data message
        catch 
            disp(['An error occurred while importing data from ' fnames{i}]);
            disp('Skipping this file.')
            continue
        end
    end
    cd(currdir)
end
