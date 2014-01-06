classdef Brain < handle
    %{
    sim = pitt.Simu( 'subj', '~/simu/avniel_clone', 'LS1-lh.label', 'RS1-rh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
    sim.post_simulate();
    coords_lh = pitt.exp.Brain.prepare_lh( sim ); close all;
    coords_rh = pitt.exp.Brain.prepare_rh( sim ); close all;
    figure; pitt.exp.Brain.disp_all( coords_lh, coords_rh );
    
    L1 = sim.dipole_loc_label1;  
    L2 = sim.dipole_loc_label2 + 4098; 
    
    freq15 = load( 'freq_015' ); freq15 = freq15.output;
    Afreq15 = freq15 + tril( freq15, 1 )';
    
    %}
    methods ( Static )
        
        function coords_lh = prepare_lh( sim )
           
            %[v,f]       = pitt.exp.simu.GraphAnalysis.loadFaceVertexData( '/synapse/labresources/sMRI/Opt065/surf/lh.pial' ); 
            [v,f]       = pitt.exp.simu.GraphAnalysis.loadFaceVertexData( '/synapse/labresources/sMRI/Opt065/surf/lh.inflated' ); 
            h           = pitt.exp.simu.GraphAnalysis.displayBrain( v, f );
            all_coords  = get( h, 'Vertices' );
            
            coords_lh   = all_coords( sim.lhdata.vertices,: );
            
        end
        
        function coords_rh = prepare_rh( sim )
           
            [v,f]       = pitt.exp.simu.GraphAnalysis.loadFaceVertexData( '/synapse/labresources/sMRI/Opt065/surf/rh.inflated' ); 
            h           = pitt.exp.simu.GraphAnalysis.displayBrain( v, f );
            all_coords  = get( h, 'Vertices' );
            
            coords_rh   = all_coords( sim.rhdata.vertices,: );
            
        end
        
        function disp_lh( coords, labelsize, labelcolor )
            
            numpoints = size( coords,1 );
            
            if( nargin == 1 )
                scatter3( coords(:,1), coords(:,2), coords(:,3) );
            elseif( nargin == 2 )
                if( isempty( labelsize ) )
                    labelsize = ones( numpoints, 1 );
                end
                scatter3( coords(:,1), coords(:,2), coords(:,3), labelsize );
            elseif( nargin == 3 )
                scatter3( coords(:,1), coords(:,2), coords(:,3), labelsize, labelcolor );
            end
            view( -90, 0 );
            
        end
        
        function disp_rh( coords, labelsize, labelcolor )
           
            
            numpoints = size( coords,1 );
            
            if( nargin == 1 )
                scatter3( coords(:,1), coords(:,2), coords(:,3) );
            elseif( nargin == 2 )
                if( isempty( labelsize ) )
                    labelsize = ones( numpoints, 1 );
                end
                scatter3( coords(:,1), coords(:,2), coords(:,3), labelsize );
            elseif( nargin == 3 )
                scatter3( coords(:,1), coords(:,2), coords(:,3), labelsize, labelcolor );
            end
            
            view( 90, 0 );
            
        end
        
        function disp_all( coords_lh, coords_rh )
            
            subplot( 1,2,1 );
            pitt.exp.Brain.disp_lh( coords_lh );
            
            subplot( 1,2,2 );
            pitt.exp.Brain.disp_rh( coords_rh );
        end
        
        function disp_all_with_labels( coords_lh, coords_rh, L1, L2 )
            
            subplot( 1,2,1 );
            pitt.exp.Brain.disp_lh( coords_lh, ones( size(coords_lh,1),1 ) );
            hold on;
            scatter3( coords_lh(L1,1),coords_lh(L1,2),coords_lh(L1,3),'*');
            hold off;
            
            subplot( 1,2,2 );
            pitt.exp.Brain.disp_rh( coords_rh, ones( size(coords_rh,1),1 ) );
            hold on;
            scatter3( coords_lh(L2-4098,1),coords_lh(L2-4098,2),coords_lh(L2-4098,3),'*');
            hold off;
            
        end
        
        function overlay_values( coords_lh, coords_rh, values, seed )
            
            if( nargin == 4 )
                
                subplot( 3,6,[1 2 3 7 8 9] );
                pitt.exp.Brain.disp_lh( coords_lh, values(1:4098,:), values(1:4098,:) );
                hold on; 
                scatter3( seed(:,1), seed(:,2), seed(:,3), 100, 'filled' );
                hold off;
                
                subplot( 3,6, [13 14 15] )
                pitt.exp.Brain.disp_lh( coords_lh, values(1:4098,:), values(1:4098,:) ); view( -90, 90 );
                hold on; 
                scatter3( seed(:,1), seed(:,2), seed(:,3), 100, 'filled' );
                hold off;
                
                subplot( 3,6,[4 5 6 10 11 12] );
                pitt.exp.Brain.disp_rh( coords_rh, values(4099:end,:), values(4099:end,:) );
                
                subplot( 3,6,[16 17 18] );
                pitt.exp.Brain.disp_rh( coords_rh, values(4099:end,:), values(4099:end,:) ); view( -90, 90 );
            else
                
                subplot( 1,2,1 );
                pitt.exp.Brain.disp_lh( coords_lh, values(1:4098,:), values(1:4098,:) );

                subplot( 1,2,2 );
                pitt.exp.Brain.disp_rh( coords_rh, values(4099:end,:), values(4099:end,:) );
            
            end
            
            
        end
        
        %{ plot vertex indexes on the LH %}
        function identify_vertices_lh( coords_lh, vertices_lh )
            subplot( 2,3,[ 1 2 3] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3) );
            hold on;
            scatter3( coords_lh(vertices_lh,1), coords_lh(vertices_lh,2), coords_lh(vertices_lh,3), 300.*ones(size(vertices_lh,2),1), 'filled' );
            hold off;
            view( -90, 0 );
            
            subplot( 2,3,[4 5 6] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3) );
            hold on;
            scatter3( coords_lh(vertices_lh,1), coords_lh(vertices_lh,2), coords_lh(vertices_lh,3), 300.*ones(size(vertices_lh,2),1), 'filled' );
            hold off;
            view( -90, 90 );
        end
        
        
        function disp_movie( coords_lh, coords_rh, values, jump_points, timeout )
            
            if( nargin == 4 )
                timeout = 2;
            end
            
            for i=1:length(jump_points); 
                figure(1); 
                pitt.exp.Brain.overlay_values( coords_lh, coords_rh, values( :, jump_points(i) ), coords_lh(jump_points(i), : ) ); 
                pause(timeout); 
            end;
        end
        
        %{
        
        Followed by: 
        movie2avi( M, 'myavi.avi' );
        
        or use:
        movie2avi( M, 'myavi3.avi', 'FPS', 1 ) 
        to adjust framerate
        
        %}
        function M = make_movie( coords_lh, coords_rh, values, jump_points )
            
            for i=1:length(jump_points); 
                figure(1); 
                pitt.exp.Brain.overlay_values( coords_lh, coords_rh, values( :, jump_points(i) ), coords_lh(jump_points(i), : ) ); 
                M(i) = getframe(1);
            end;
            
        end
        
        
    end
    
    
end