classdef GCScales < handle
    
    methods (Static)
        
        function output = runGC( data, numScales, numlags )
            
            kernel = pitt.exp.GCScales.genGCPipe( numlags );
            kernel.initial( data );
            output = kernel.execute();
            return;
            
            s = pitt.jobs.Scale( data, pitt.exp.GCScales.genGCPipe( numlags ) );
            %warning off all;
            output = s.runByWindow( numScales );
            %warning on all;
        end
        
        function pipe = genGCPipe( numlags )
            pipe = pitt.jobs.Pipe;
            %pipe.add( @pitt.exp.GCScales.splitGC, 'X' );
            pipe.add( @cca_detrend, 'X' );
            pipe.add( @cca_rm_temporalmean, 'X' );
            pipe.add( @cca_granger_regress, 'X', numlags );
        end
        
        function X = splitGC( input )
            if( mod( length(input), 2 ) ~= 0 )
                error('Could not divide in half the input data' );
            end
            X = [input(1:end/2); input((end/2)+1:end)];
        end
        
    end
    
end
