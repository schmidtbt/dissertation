classdef Pipe < handle
    
    properties
        functions;
        arguments;
        data;
    end
    
    methods
        
        function add( obj, fcn, varargin )
            obj.functions{ length(obj.functions)+1 } = fcn;
            obj.arguments{ length(obj.arguments)+1 } = varargin;
        end
        
        function initial( obj, input )
            obj.data = input;
        end
        
        function output = execute( obj )
            for i = 1:length( obj.functions )
                args = obj.arguments{i};
                for j = 1:length(args)
                    if( args{j} == 'X' )
                        args{j} = obj.data;
                    end
                end
                obj.data = obj.functions{i}( args{:} );
            end
            output = obj.data;
        end
        
        function report( obj )
            
            for i = 1:length( obj.functions )
                disp('---------------------');
                disp( obj.functions{i});
                
                for j = 1:length( obj.arguments{i} )
                    disp( obj.arguments{i}{j})
                end
                
            end
        end
        
    end
    
end