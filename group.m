classdef group < handle
    
    properties
        
        datadir
        
        sub
        
        
    end
    
    methods
        
        function obj = group(datadir)
            
            obj.datadir = datadir;
        end
        
        function load(obj)
            
            d = dir([obj.datadir '*.mat']);
            
            for sn = 1:length(d)
                disp(['subject ' num2str(sn) ' of ' num2str(length(d))])
                disp([d(sn).name])
                sub(sn) = subject(obj.datadir, d(sn).name);
                sub(sn).load;
                sub(sn).load_eye;
            end
            
            obj.sub = sub;
            
        end
        
        function clean_pupil(obj)
            for sn = 1:length(obj.sub)
                obj.sub(sn).clean_pupil;
            end
        end
        
    end
    
end