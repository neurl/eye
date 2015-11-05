 function result = importeyedata(filename)
    %filename = 'sampledata.tsv';
    fid = fopen(filename);
    X = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    Y = strvcat(X{1});
    ind_mess = Y(:,1) == 'M';
    M = {X{1}{ind_mess}};
    Z = {X{1}{~ind_mess}};
    eye.data.var_names = strsplit(strvcat(Z(1)))';
    Z = {Z(2:end)};
    A = strvcat(Z{:});
    B = textscan(A','%s%s%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
    eye.data.date = B{1};
    for i = 3:size(B,2)
        ts = eye.data.var_names(i-1);
        eye.data.(ts{1}) = B{i};
    end
    E = strvcat(B{2});
    sc = str2num(E(:,end-5:end));
    mn = str2num(E(:,end-8:end-7));
    hr = str2num(E(:,end-11:end-10));
    tm = hr*3600+mn*60+sc;
    eye.data.hour = hr;
    eye.data.minute = mn;
    eye.data.second = sc;
    eye.data.times = tm;
    C = strvcat(M{:});
    C(:,size(C,2)+1) = repmat(' ',size(C,1),1);
    D = textscan(C','%s%s%s%s%s');
    E = strvcat(D{3});
    sc = str2num(E(:,end-5:end));
    mn = str2num(E(:,end-8:end-7));
    hr = str2num(E(:,end-11:end-10));
    tm = hr*3600+mn*60+sc;
    eye.message.date = D{2};
    eye.message.second = sc;
    eye.message.minute = mn;
    eye.message.hour = hr;
    eye.message.times = tm;
    eye.message.time = D{4};
    eye.message.string = D{5};
    result = eye;
 end