classdef Scale < handle
    
    properties
        data;
        pipe;
        output;
    end
    
    
    methods
        
        function obj = Scale( data, pipeObj )
            obj.data = data;
            obj.pipe = pipeObj;
            obj.output = cell(0);
        end
        
        function output = runByChannel( obj )
            for i = 1:size( obj.data, 1 )
                    obj.pipe.initial( obj.data(i,:) );
                    obj.output{i} = obj.pipe.execute();
            end
            
            output = obj.output;
        end
        
        %output [Channel_idx][scale_idx][window_idx]
        function output = runByWindow( obj, numScales, REPORTING )
            
            if( nargin == 2 )
                REPORTING = 0;
            end
            
            maxvertex = size( obj.data, 1 );
            tic;
            for i = 1:size( obj.data, 1 )
                    if( REPORTING ); 
                        fprintf('Vertex: %i / %i\n', i, maxvertex); 
                        t = toc;
                        fprintf(' Cycle time:(in s): %i\n', toc); 
                        fprintf(' (Est)Completion (in s): %i\n', pitt.jobs.Scale.computeElapsedFinish( toc, i, maxvertex )); 
                    end;
                    win = pitt.signal.Windows( obj.data(i,:), numScales );
                    for j = 1:numScales
                        if( REPORTING ); fprintf('\tscale: %i\n', j); end;
                        windat = win.getLevelData(j-1);
                        for k = 1:size(windat,1)
                            obj.pipe.initial( windat(k,:) );
                            obj.output{i}{j}{k} = obj.pipe.execute();
                        end
                    end
            end
            
            output = obj.output;
        end
        
        function output = runByPairedWindow( obj, numScales )
            
            %split the data
            if( mod( length(obj.data), 2 ) ~= 0 )
                error('Could not divide in half the input data' );
            end
            
            bidata = [obj.data(1:end/2); obj.data((end/2)+1:end)];
            
            for i = 1%:size( obj.data, 1 )
                
                win1 = pitt.signal.Windows( bidata(1,:), numScales );
                win2 = pitt.signal.Windows( bidata(2,:), numScales );
                for j = 1:numScales
                    windat1 = win1.getLevelData(j-1);
                    windat2 = win2.getLevelData(j-1);
                    for k = 1:size(windat1,1)
                        obj.pipe.initial( [windat1(k,:); windat2(k,:)] );
                        obj.output{i}{j}{k} = obj.pipe.execute();
                    end
                end
            end
            
            output = obj.output;
        end
        
        
    end
    
    methods( Static )
        
        function est = computeFinishTime( tocCycle, indexCurrent, indexOfTotal )
            est = (indexOfTotal - indexCurrent )*tocCycle;
        end
        
        function est = computeElapsedFinish( tocCycle, indexCurrent, indexOfTotal )
            est = (tocCycle/indexCurrent)*(indexOfTotal - indexCurrent );
        end
        
    end
    
    
end