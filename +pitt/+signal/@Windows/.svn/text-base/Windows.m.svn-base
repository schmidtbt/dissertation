classdef Windows < handle
    
    properties
        chanData = [];
        numscales = 1;
    end
    
    methods
        
        function obj = Windows( Channel, NumScales )
            obj.chanData = Channel;
            obj.numscales = NumScales;
        end
        
        function [levelMatrix, maxIdx, stepSize, numSteps] = getLevelData( obj, level_idx )
            
            if( level_idx > obj.numscales )
                error('Index is greater than scales specified');
            end
            
            [maxIdx, stepSize, numSteps] = pitt.signal.Windows.getResolutionLevel( length( obj.chanData ), level_idx );
            
            levelMatrix = [];
            for i = 0:numSteps-1
                levelMatrix = [levelMatrix; obj.chanData( (i*stepSize)+1:(i*stepSize)+stepSize )];
            end
            
        end
        
    end
    
    
    methods (Static)
        function [maxIdx, stepSize, numSteps] = getResolutionLevel( lengthData, Scale )
            
            winSize = 2.^Scale;
            stepSize = floor( lengthData ./ winSize );
            %stepSize = 2.^Scale;
            if( stepSize > lengthData )
                error('This scale is too large for this data');
            end
            remain = mod( lengthData, stepSize );
            if( remain )
                % Not evenly divisible
                maxIdx = lengthData - remain;
            else
                % Evenly divisible
                maxIdx = lengthData;
            end
            
            numSteps = maxIdx ./ stepSize;
        end
    end
    
end