classdef Processing < handle
   
    methods( Static )
       
        % Modified kernel from pitt.exp.GTesting.gregressKernel allowing a
        % parameterized number of lags and also extracting the prb and gc
        % values of the cca_granger_regress output for easier loading and
        % smaller storage overhead.
        function kernel = gregressKernel( lags )
            
            if( nargin == 0 )
                lags = 10;
            end
            
            kernel = legion.Kernel;
            kernel.add( @pitt.signal.Util.varnorm, 'X' );
            kernel.add( @cca_granger_regress, 'X', lags );
            kernel.add( @pitt.exp.simu.Processing.extractPrbGCFromRet, 'X' );
            
            kernel.addEnviro( @pitt.Depend.CCAadd );
        end
        
        % Given the output of cca_granger_regress, extract a structure
        % contianing only the prb and gc values
        function output = extractPrbGCFromRet( ret )
            output.prb = ret.prb;
            output.gc  = ret.gc;
        end
        
        % Get all the label data
        % stackes label1 ; label2
        function [label_data, label_vertices] = extractLabelDataAllAll( sim )
            label_data        = [sim.label1_data; sim.label2_data ];
            label_vertices    = [sim.dipole_loc_label1, sim.dipole_loc_label2 ];
        end
        
        % Get data from random 'number' of vertices from the data set
        function [rand_data, rand_vertices] = getRandomVertices( sim, number )
            max_vertices    = size( sim.lhdata.data, 1 );
            rand_vertices   = randi( max_vertices, 1, number );
            rand_data       = sim.lhdata.data( rand_vertices, : );
        end
        
        % Get random 'number' of vertices from the data set with subset of
        % timepoints
        % Where tpt_idxs is an array of time-point indexes (eg [1000:1100]
        % );
        function [rand_data, rand_vertices] = getRandomVerticesSubTpts( sim, number, tpt_idxs )
            max_vertices    = size( sim.lhdata.data, 1 );
            rand_vertices   = randi( max_vertices, 1, number );
            rand_data       = sim.lhdata.data( rand_vertices, tpt_idxs );
        end
        
        
        
        
    end
    
end