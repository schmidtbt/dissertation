classdef GraphAnalysis < handle

    
    
    methods (Static)
        
        
        %{
        Load vertices and face data
        
        typical
        [v,f] = pitt.exp.simu.GraphAnalysis.loadFaceVertexData( '/synapse/labresources/sMRI/Opt065/surf/lh.pial' );
        
        We auto-increment f from zero-based to 1-based
        %}
        
        function [v,f] = loadFaceVertexData( surface )
            
            pitt.Depend.FreeSurferAdd();
            [v,f]=read_surf( surface );
            f = f+1;
            
        end
        
        function h = displayBrain( v, f )
           
            pitt.Depend.SurfaceCodeAdd();
            h = plot_mesh( v,f );
            colormap( jet );
            
        end
        
        % displayBrain but oriented for LH view
        function h = displayLHBrain( v,f )
            h = pitt.exp.simu.GraphAnalysis.displayBrain( v,f );
            view( -90, 0 );
            %{
            subplot( 3,1,[1 2] );
            h = pitt.exp.simu.GraphAnalysis.displayBrain( v,f );
            view( -90, 0 );
            subplot( 3,1,2 );
            h = pitt.exp.simu.GraphAnalysis.displayBrain( v,f );
            view( -90, 90 );
            %}
        end
        
        function newVertexValues( h, new_values )
            
            cur_data = get(h,'FaceVertexCData');
            
            if( length(cur_data) ~= length(new_values) )
                error( sprintf('Wrong dimensions specified, should be: %i x 1, given: %i x %i', length(cur_data), size(new_values,1), size(new_values,2)) );
            end
            
            set(h,'FaceVertexCData', new_values );
            
        end
        
        function overlayPLVData( h, plv_data )
            
            if( size(plv_data,1) < size(plv_data,2) )
                plv_data = plv_data';
            end
            
            lst = find( plv_data <= 0 );
            plv_data(lst) = 0;
            
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, plv_data );
            
        end
        
        function displayVertexLocation( h, vertex_idx )
            
            cur_data = get(h,'FaceVertexCData');
            
            vd = zeros( length( cur_data ),1 );
            vd(vertex_idx) = 100;
            
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, vd );
            
        end
        
        
        function vsize = getVertexDataSize( h )
            vsize = length(get(h,'FaceVertexCData'));
        end
        
        %{
        Display the positions of sim vertex data onto brain.
        Returns the new vertice value data
        %}
        function vd = displaySimVertices( h, sim )
           
            vd = zeros( pitt.exp.simu.GraphAnalysis.getVertexDataSize(h),1);
            vd(sim.lhdata.vertices) = 1;
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, vd );
            
        end
        
        function displaySimLabel1( h, sim )
            vd = pitt.exp.simu.GraphAnalysis.displaySimVertices( h, sim );
            vd(sim.lhdata.vertices(sim.dipole_loc_label1)) = 150;
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, vd );
        end
        
        
        function displaySimLabel2( h, sim )
            vd = pitt.exp.simu.GraphAnalysis.displaySimVertices( h, sim );
            vd(sim.lhdata.vertices(sim.dipole_loc_label2)) = 150;
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, vd );
        end
        
        
        function displayVertexDataOnSim( h, sim, new_data )
            vd = pitt.exp.simu.GraphAnalysis.displaySimVertices( h, sim );
            vd( sim.lhdata.vertices ) = ceil(new_data)';
            pitt.exp.simu.GraphAnalysis.newVertexValues( h, vd );
        end
        
        function bounded_plv = boundPLV( plv_data, lower, upper )
            
            lst = find( plv_data > lower & plv_data < upper );
            z = zeros( 1, length(plv_data) );
            z(lst) = plv_data(lst);
            
            bounded_plv = z;
            
        end
        
        function bounded_plv_2d = boundPLVMap( plv_data_2d, lower, upper )
            
            lst = find( plv_data_2d > lower & plv_data_2d < upper );
            bounded_plv_2d = zeros(size(plv_data_2d,1), size(plv_data_2d,2));
            bounded_plv_2d(lst) = plv_data_2d(lst);
            
        end
        
        %{
        
        Want to extract the XYZ coordiantes of every LH vertex in sim
        
        Return indexed from 1:4098 with [X Y Z] coords
        
        %}
        function coords = getXYZCoordsOfVertices( v, sim )
            coords = v( sim.lhdata.vertices, : );
        end
        
        
        %{
        Extract the L2 norm (Euclidean) distance from every coordinate
        to every other coordinate
        %}
        function output = getL2NormPairWise( XYZcoords )
            
            [rows, cols] = size( XYZcoords );
            
            if( cols ~= 3 )
                error( 'XYZ requires 3 points. Not enough values provided' );
            end
            
            % initalize an output variable
            output = zeros( rows, rows );
            
            for i = 1:rows
                disp(i);
                for j =1:rows
                
                    if( i > j ); continue; end;
                    
                    
                    output(i,j) = pdist2(XYZcoords(i,:),XYZcoords(j,:),'euclidean');
                    
                end
                
            end
            
            % Complete flip of matrix:
            output = triu(output)+triu(output,1)';
            
        end
        
        %{
        Get the idxs,idxs points in a 2D matrix.
        %}
        function subset = getSubsetVerticesOfMatrix( idxs, matrix )
            
            subset = matrix( idxs,  : );
            subset = subset( :, idxs );
            
        end
        
        function edgeList = convertIncidenceMatrixToEdgeList( matrix )
            
            [row, col] = size( matrix );
            
            edgeList = [];
            
            for i = 1:row
                disp(i);
                for j = 1:col
                    
                    if( matrix(i,j) > 0 )
                        edgeList = [edgeList; [i j matrix(i,j)] ];
                    end
                    
                end
                
            end
            
        end
        
    end
    
end