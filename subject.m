classdef subject < handle
    
    properties
        
        datadir
        dataname
        
        data
        
        eye_name
        eye_data
        eye_mess
        
    end
    
    methods
        
        function obj = subject(datadir, dataname)
            obj.datadir = datadir;
            obj.dataname = dataname;
            obj.eye_name = [obj.dataname(1:end-4) '_eye.tsv'];
        end
        
        function load(obj)
            
            % load behavior
            load([obj.datadir obj.dataname]);
            for i = 1:length(data)
                data(i).trialNum = i;
            end
            obj.data = data;
            
        end
        
        function load_eye(obj)
            
            % load eye data
            fid = fopen([obj.datadir obj.eye_name]);
            X = textscan(fid, '%s', 'delimiter', '\n');
            fclose(fid);
            
            T = textscan(X{1}{1}, '%s', 'delimiter', '\t');
            Y = strvcat(X{1}{:});
            ind_mess = Y(:,1) == 'M';
            M = {X{1}{ind_mess}};
            Z = {X{1}{~ind_mess}};
            Z = {Z{2:end}};
            A = strvcat(Z{:});
            B = textscan(A','%s%s%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
            
            E = strvcat(B{2});
            sc = str2num(E(:,end-5:end));
            mn = str2num(E(:,end-8:end-7))*60;
            hr = str2num(E(:,end-11:end-10))*60*60;
            tm = hr+mn+sc;
            
            eye_data.tm = tm;
            eye_data.pd_raw = [B{15} B{22}];
            
            eye_data.var_names = {T{1}{4:end}};
            
            
            % messages
            C = strvcat(M{:});
            C = C(2:end,:);
            CC = C(:,4:end);
            D = textscan(CC','%s%s%s%s');
            E = strvcat(D{2});
            sc = str2num(E(:,end-5:end));
            mn = str2num(E(:,end-8:end-7))*60;
            hr = str2num(E(:,end-11:end-10))*60*60;
            tm = hr+mn+sc;
            
            % remove events after pupil loss
            
            %if ~isempty(min(find(diff(tm)==0)))
            %   idx = 1:min(find(diff(tm)==0));
            %   tm = tm(idx);
            %   str = {D{4}{idx}};
            %else
                str = D{4};
            %end
            mess.tm = tm;
            mess.str = str;
            
            
            obj.eye_data = eye_data;
            obj.eye_mess = mess;
            
        end
        
        function clean_pupil(obj)
            obj.eye_data.pd = clean_pupil(obj.eye_data.pd_raw);
        end
        
    end
    
end