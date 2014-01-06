classdef Job < handle
    
    properties
        jobName = 'Default';
    end
    
    methods (Static)
        function exec
            disp('here');
        end
    end
    
    methods
        
        function obj = Job( Name )
            if( nargin < 1 )
                Name = 'Default';
            end
            obj.jobName = Name;
        end
        
        function bool = execute(obj)
            bool = false;
            obj.exec;
            bool = true;
        end
    end
    
end