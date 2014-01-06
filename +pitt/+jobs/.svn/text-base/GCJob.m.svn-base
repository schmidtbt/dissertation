classdef GCJob < pitt.Job
    
    properties
        intValue = 3;
        modelOrder = 5;
        x1;
        x2;
        ret;
        idx1;
        idx2;
    end
    
    methods
        function obj = GCJob( Name, Idx1, Idx2, Chan1, Chan2, ModelOrder )
            % Example overriding parent Job class
            if( nargin < 1 )
                Name = 'DefaultGCJob';
            end
            newName = ['GCJob | ', Name];
            obj@pitt.Job( newName );
            
            obj.modelOrder = ModelOrder;
            obj.x1 = Chan1;
            obj.x2 = Chan2;
            obj.idx1 = Idx1;
            obj.idx2 = Idx2;
            
            
        end
        
        
        function bool = execute( obj )
            retI = cca_granger_regress( [obj.x1; obj.x2], obj.modelOrder, 1 );
            obj.ret = retI;
            %disp( obj )
            bool = true;
        end
    end
    
end