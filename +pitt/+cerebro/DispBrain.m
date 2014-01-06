classdef DispBrain < handle
    
    methods (Static)
        
        % Disp LH/RH brain with no color bar
        function [hL, hR] = dispFullBrain()
            
            left_vf     = load( '~/opt001_v_f.mat' );
            right_vf    = load( '~/opt001_v_f_right.mat' );
            
            figure;
            subplot( 1,2,1 );
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( left_vf.vertices, left_vf.faces ); view( -90, 0 );
            subplot( 1,2,2 );
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( right_vf.vertices, right_vf.faces ); view( 90, 0 );
            
        end
        
        % Disp LH/RH brain with colorbar
        function [hL, hR] = dispFullBrainWithColorbar()
            
            left_vf     = load( '~/opt001_v_f.mat' );
            right_vf    = load( '~/opt001_v_f_right.mat' );
            
            figure;
            
            
            subplot( 3,6,[1 2 3 7 8 9] );
            colorbar;
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( left_vf.vertices, left_vf.faces ); view( -90, 0 );
            subplot( 3,6,[4 5 6 10 11 12] );
            colorbar;
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( right_vf.vertices, right_vf.faces ); view( 90, 0 );
            
        end
        
        function dispData( hL, hR, data )
            
            if( size(data,1) > size(data,2) )
                data = data';
            end
            
            if( size(data,2) ~= 8196 )
                error('size must be 1 x 8196' )
            end
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, data(1:4098) )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, data(4099:end) )
            
        end
        
        function dispDataCNorm( hL, hR, data, lh_caxis, rh_caxis )
            
            if( size(data,1) > size(data,2) )
                data = data';
            end
            
            if( size(data,2) ~= 8196 )
                error('size must be 1 x 8196' )
            end
            
            subplot( 3,6,[1 2 3 7 8 9] );
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, data(1:4098) )
            caxis(lh_caxis);
            subplot( 3,6,[4 5 6 10 11 12] );
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, data(4099:end) )
            caxis(rh_caxis);
            
        end
        
        function dispVertices( hL, hR, vertices )
            
            lstL = find( vertices <= 4098 );
            lstR = find( vertices > 4098 );
            
            pitt.exp.simu.GraphAnalysis.displayVertexLocation( hL, vertices(lstL)' );
            pitt.exp.simu.GraphAnalysis.displayVertexLocation( hR, vertices(lstR) - 4098 );
            
        end
        
        
        % --------------------------------
        % Methods for Flood Fill Brains
        % High-reso mesh with lower vertex sampling (NIH dataset)
        % --------------------------------
        
        % Extract info from fwd or inv *.fif file
        function [LH_lookup_idx, LH_v, LH_f, RH_lookup_idx, RH_v, RH_f, sizeLsubmesh, sizeRsubmesh, LH_v_used, RH_v_used ] = generateMeshLookupTable( fif_fwd_file )
            
            s = mne_read_source_spaces( fif_fwd_file );
            
            used            = s(1).rr( find( s(1).inuse ), : );
            fprintf( 'generating LH lookup table\n' );
            LH_lookup_idx   = dsearchn( used, s(1).rr );
            LH_v            = s(1).rr;
            LH_f            = s(1).tris;
            sizeLsubmesh    = length(find( s(1).inuse ));
            LH_v_used       = used;
            
            used            = s(2).rr( find( s(2).inuse ), : );
            fprintf( 'generating RH lookup table\n' );
            RH_lookup_idx   = dsearchn( used, s(2).rr );
            RH_v            = s(2).rr;
            RH_f            = s(2).tris;
            sizeRsubmesh    = length(find( s(2).inuse ));
            RH_v_used       = used;
            
        end
        
        % Display blank brains and return figure handles for further displaying
        % Displays both hemispheres
        function [hL, hR] = dispFloodBrain( LH_v, LH_f, RH_v, RH_f )
            
            figure;
            subplot( 1,2,1 );
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( LH_v, LH_f ); view( -90, 0 );
            subplot( 1,2,2 );
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( RH_v, RH_f ); view( 90, 0 );
            
        end
        
        % Place data onto flood brain
        % Data will be at lower spatial resolution. Lookup tables are used to flood-fill brain
        % SplitIdx is used to delineate which sub-sample data vertices are in LH vs RH. This value is last vertex point in LH
        function dispFloodData( hL, hR, lookupL, lookupR, splitIdx, data )
            
            floodedL = pitt.cerebro.DispBrain.executeFloodLookup( lookupL, data(1:splitIdx) );
            floodedR = pitt.cerebro.DispBrain.executeFloodLookup( lookupR, data(splitIdx+1:end) );
            
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, floodedL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, floodedR )
            
        end
        
        function dispFullFloodData( hL, hR, splitIdx, data )
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, data(1:splitIdx) )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, data(splitIdx+1:end) )
            
        end
        
        function flooded = executeFloodLookup( lookupTable, data )
            
            flooded = zeros( 1, length(lookupTable ) );
            
            for i = 1:length(data)
                
                lst = find( lookupTable == i );
                flooded( lst ) = data( i );
                
            end
            
        end
        
        function downsampled = executeFloodDownsample( lookupTable, floodedData )
            
            downsampled = zeros( 1, length(unique(lookupTable )) );
            
            for i=1:length(downsampled)
                
                lst = find( lookupTable == i );
                downsampled = mean( floodedData( lst ) );
                
            end
            
        end
        
        function downsampled = executeDownsample( lh_lookup, rh_lookup, fullFloodData )
            
            downSamSize = length(unique(lh_lookup))+length(unique(rh_lookup ));
            downsampled = zeros( 1, downSamSize );
            
            full_lookup = [ lh_lookup; rh_lookup+max(unique(lh_lookup))+1];
            
            for i = 1:length(unique( fullFloodData ))
                
                lst = find( fullFloodData == i );
                downsampled(unique(full_lookup(lst))) = mean( fullFloodData( lst ) );
                
            end
            
        end
        
        
        % --------------------------------
        % Methods for use with Opt* structures of pitt.cerebro.Reg
        % --------------------------------
        
        function [hL, hR] = displaySurface( surface_struct )
            
            figure;
            subplot( 1,2,1 );
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.L.v, surface_struct.surface.L.f ); view( -90, 0 );
            subplot( 1,2,2 );
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.R.v, surface_struct.surface.R.f ); view( 90, 0 );
            
            surface_struct.handle.L = hL;
            surface_struct.handle.R = hR;
        end
        
        % Plot surface and store the handles within the opt object
        function opt = outputSurface( opt, surface_name )
            
            surface_struct = opt.(surface_name);
            
            h = figure;
            set(h, 'color', 'w');
            m = colormap;
            m(1,:) = .8; colormap(m)
            
            subplot( 1,2,1 );
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.L.v, surface_struct.surface.L.f ); view( -90, 0 );
            subplot( 1,2,2 );
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.R.v, surface_struct.surface.R.f ); view( 90, 0 );
            
            opt.(surface_name).handle.L = hL;
            opt.(surface_name).handle.R = hR;
            opt.(surface_name).handle.h = h;
        end
        
        % Plot surface and store the handles within the opt object
        function opt = outputSurfaceImgs( opt, surface_name )
            
            surface_struct = opt.(surface_name);
            
            h = figure;
            set(h, 'color', 'w');
            m = colormap;
            m(1,:) = .8; colormap(m)
            
            subplot( 2,2,1 );
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.L.v, surface_struct.surface.L.f ); 
            view( -90, 0 );
            camlight;
            subplot( 2,2,2 );
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.R.v, surface_struct.surface.R.f ); 
            view( 90, 0 );
            camlight;
            subplot( 2,2,3 );
            hLf = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.L.v, surface_struct.surface.L.f ); 
            view( -170, 0);
            camlight;
            subplot( 2,2,4 );
            hRf = pitt.exp.simu.GraphAnalysis.displayBrain( surface_struct.surface.R.v, surface_struct.surface.R.f ); 
            view( 170 , 0);
            camlight;
            
            m = colormap; m(1,:) = .8; colormap(m);
            
            opt.(surface_name).handle.L = hL;
            opt.(surface_name).handle.R = hR;
            opt.(surface_name).handle.Lf = hLf;
            opt.(surface_name).handle.Rf = hRf;
            opt.(surface_name).handle.h = h;
        end
        
        
        function opt = outputSingleSurface( opt, surface_name )
            
            surface_struct = opt.(surface_name);
            
            h = figure;
            set(h, 'color', 'w');
            m = colormap;
            m(1,:) = .8; colormap(m)
            
            
            hL = pitt.exp.simu.GraphAnalysis.displayBrain( pitt.cerebro.DispBrain.LHtransform( surface_struct.surface.L.v ), surface_struct.surface.L.f ); 
            hold on;
            hR = pitt.exp.simu.GraphAnalysis.displayBrain( pitt.cerebro.DispBrain.RHtransform( surface_struct.surface.R.v ), surface_struct.surface.R.f ); 
            
            opt.(surface_name).handle.hLS = hL;
            opt.(surface_name).handle.hRS = hR;
            
            view(0,0);
        end
        
        function v = LHtransform( v )
            
            transX = -.05;
            transY = .2;
            angl = deg2rad( -90 );
            
            v = v + repmat( [ transX, 0, 0], size(v,1), 1 );
            v = v + repmat( [ 0, transY, 0], size(v,1), 1 );
            
            for i = 1:size(v,1)
                %v(i,:) = v(i,:)*[ cos(angl) 0 -sin(angl); 0 1 0; sin(angl) 0 cos(angl)];
                v(i,:) = v(i,:)*[ cos(angl) -sin(angl) 0; sin(angl) cos(angl) 0; 0 0 1];
            end
        end
        
        function v = RHtransform( v )
            
            transX = .05;
            transY = .2;
            angl = deg2rad( 90 );
            
            v = v + repmat( [ transX, 0, 0], size(v,1), 1 );
            v = v + repmat( [ 0, transY, 0], size(v,1), 1 );
            
            for i = 1:size(v,1)
                %v(i,:) = v(i,:)*[ cos(angl) 0 -sin(angl); 0 1 0; sin(angl) 0 cos(angl)];
                v(i,:) = v(i,:)*[ cos(angl) -sin(angl) 0; sin(angl) cos(angl) 0; 0 0 1];
            end
        end
        
        function displayRawDataOntoSurface( opt, surf_name, data )
            
            hL      = opt.(surf_name).handle.L;
            hR      = opt.(surf_name).handle.R;
            
            maxIdx  = length( opt.(surf_name).surface.L.v ) + length( opt.(surf_name).surface.R.v );
            splitIdx= length( opt.(surf_name).surface.L.v );
            
            if( length(data) ~= maxIdx )
                maxIdx
                error('Data is wrong size, expected: 1 x %i\n', maxIdx);
            end
            
            vdL = data(1:splitIdx);
            vdR = data(splitIdx+1:end);
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            
        end
        
        function displayRawDataOntoSurfaceImgs( opt, surf_name, data )
            
            hL      = opt.(surf_name).handle.L;
            hR      = opt.(surf_name).handle.R;
            hLf     = opt.(surf_name).handle.Lf;
            hRf     = opt.(surf_name).handle.Rf;
            
            maxIdx  = length( opt.(surf_name).surface.L.v ) + length( opt.(surf_name).surface.R.v );
            splitIdx= length( opt.(surf_name).surface.L.v );
            
            if( length(data) ~= maxIdx )
                maxIdx
                error('Data is wrong size, expected: 1 x %i\n', maxIdx);
            end
            
            vdL = data(1:splitIdx);
            vdR = data(splitIdx+1:end);
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hLf, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hRf, vdR )
            
        end
        
        function displayRawDataOntoSingleSurface( opt, surf_name, data )
            
            hL      = opt.(surf_name).handle.hLS;
            hR      = opt.(surf_name).handle.hRS;
            
            maxIdx  = length( opt.(surf_name).surface.L.v ) + length( opt.(surf_name).surface.R.v );
            splitIdx= length( opt.(surf_name).surface.L.v );
            
            if( length(data) ~= maxIdx )
                maxIdx
                error('Data is wrong size, expected: 1 x %i\n', maxIdx);
            end
            
            vdL = data(1:splitIdx);
            vdR = data(splitIdx+1:end);
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            
        end
        
        % Display lower resolution data onto higher resolution brain space
        function displaySurfaceDataOntoHiRes( opt, loRes_surfName, hiRes_surfName, data )
            
            hiRes_struct    = opt.(hiRes_surfName);
            loRes_struct    = opt.(loRes_surfName);
            
            splitIdx        = length( loRes_struct.surface.L.v );
            maxIdx          = length( loRes_struct.surface.L.v ) + length( loRes_struct.surface.R.v );
            
            if( length(data) ~= maxIdx )
                maxIdx
                error('Data is wrong size, expected: 1 x %i\n', maxIdx);
            end
            
            hL = opt.(hiRes_surfName).handle.L;
            hR = opt.(hiRes_surfName).handle.R;
            
            vdL = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( loRes_struct.hiRes.L, data(1:splitIdx) );
            vdR = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( loRes_struct.hiRes.R, data(splitIdx+1:end) );
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            
        end
        
        function displaySurfaceDataOntoHiResRegistered( opt, optCommon, loRes_surfName, hiRes_surfName, data )
            
            hiRes_struct    = opt.(hiRes_surfName);
            loRes_struct    = opt.(loRes_surfName);
            
            splitIdx        = length( loRes_struct.surface.L.v );
            maxIdx          = length( loRes_struct.surface.L.v ) + length( loRes_struct.surface.R.v );
            
            if( length(data) ~= maxIdx )
                maxIdx
                error('Data is wrong size, expected: 1 x %i\n', maxIdx);
            end
            
            hL = optCommon.(hiRes_surfName).handle.L;
            hR = optCommon.(hiRes_surfName).handle.R;
            
            vdL = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( loRes_struct.reg.colin.L, data(1:splitIdx) );
            vdR = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( loRes_struct.reg.colin.R, data(splitIdx+1:end) );
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            
        end
        
        % Upsample lower resolution data into higher resolution space using a lookup table.
        % Tables can be generated in opt-obj using cerebro.Reg.smoothLoResToHiResInSphere()
        function vd = overlayLoResOntoHiRes( lookUp, loResData )
            
            vd = zeros( 1,length(lookUp) );
            
            % Quickly upsample into new space
            for i = 1:length( lookUp )
                vd(i) = loResData( lookUp(i) );
            end
            return;
            
            % This is a slower big O notation algorithm
            for i = 1:length( loResData )
                
                lst = find( lookUp == i );
                vd(lst) = loResData( i );
                
            end
            
        end
        
        function vd = overlayHiRestOntoLoRes( lookUp, splitIdx, HiResData )
            
            vd = zeros( 1,splitIdx );
            u  = unique(lookUp);
            
            for i = 1:max(length(u))
                
                lst = find( lookUp == u(i) );
                if( isempty(lst) ); continue; end;
                vd(u(i)) = mean( HiResData(lst) );
                
            end
            
        end
        
        function vd = upsampleToHiResRegisteredDualHemi( reg_struct, data, splitIdx )
            
            if( nargin == 2 )
                splitIdx = max( unique( reg_struct.reg.colin.L ) );
            end
            
            vdL = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( reg_struct.reg.colin.L, data(1:splitIdx) );
            vdR = pitt.cerebro.DispBrain.overlayLoResOntoHiRes( reg_struct.reg.colin.R, data(splitIdx+1:end) );
            
            vd = [vdL, vdR];
        end
        
        function vd = downsampleHiResRegisteredToLow( opt, surf_name, data )
            
            splitIdx = length( opt.(surf_name).reg.colin.L );
            splitIdxLo = max( unique( opt.(surf_name).reg.colin.L ) );
            
            vdL = pitt.cerebro.DispBrain.overlayHiRestOntoLoRes( opt.(surf_name).reg.colin.L, splitIdxLo, data(1:splitIdx) );
            vdR = pitt.cerebro.DispBrain.overlayHiRestOntoLoRes( opt.(surf_name).reg.colin.R, splitIdxLo, data(splitIdx+1:end) );
            %vdR = zeros(1, length(opt.(surf_name).surface.R.v) );
            vd = [vdL, vdR];
        end
        
        % LowRest should contain both 'surface' and 'orig' fields -- see cerebro.Reg.appendObservedData()
        function displayRawSurfaceOntoHiRes( opt, loRes_surfName, hiRes_surfName )
            
            hiRes_struct    = opt.(hiRes_surfName);
            loRes_struct    = opt.(loRes_surfName);
            
            hL = opt.(hiRes_surfName).handle.L;
            hR = opt.(hiRes_surfName).handle.R;
            
            vdL = zeros( 1, length( hiRes_struct.surface.L.v ) );
            vdR = zeros( 1, length( hiRes_struct.surface.R.v ) );
            
            vdL( loRes_struct.orig.L.v ) = 1;
            vdR( loRes_struct.orig.R.v ) = 1;
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hL, vdL )
            pitt.exp.simu.GraphAnalysis.overlayPLVData( hR, vdR )
            
            
        end
        
        function vd = downsampleHiResRegisteredToLowByMode( opt, surf_name, data )
            
            splitIdx    = length( opt.(surf_name).reg.colin.L );
            splitIdxLo  = max( unique( opt.(surf_name).reg.colin.L ) );
            
            vdL = pitt.cerebro.DispBrain.overlayHiRestOntoLoResByMode( opt.(surf_name).reg.colin.L, splitIdxLo, data(1:splitIdx) );
            vdR = pitt.cerebro.DispBrain.overlayHiRestOntoLoResByMode( opt.(surf_name).reg.colin.R, splitIdxLo, data(splitIdx+1:end) );
            %vdR = zeros(1, length(opt.(surf_name).surface.R.v) );
            vd = [vdL, vdR];
            
        end
        
        function vd = overlayHiRestOntoLoResByMode( lookUp, splitIdx, HiResData )
            
            vd = zeros( 1,splitIdx );
            u  = unique(lookUp);
            
            for i = 1:max(length(u))
                
                lst = find( lookUp == u(i) );
                if( isempty(lst) ); continue; end;
                vd(u(i)) = mode( HiResData(lst) );
                
            end
            
        end
        
        % Return a 1D cell array with indices for each cluster_idx found of pop sample of low resolution space given hiRes cluster idxs and specified idx of interest
        % use: pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode for downsampling to lower dim space
        function plv_cluster_map_data = loResPopSampleFromHiResCluster( clusters, plv_data_2d )
            
            num_clusters = unique( clusters );
            plv_cluster_map_data = cell( length(num_clusters) );
            
            for i = 1:length(num_clusters)
                
                lst_A        = find( clusters == i );
                %row_mean    = mean( plv_data_2d(lst, : ), 1 );
                
                for j = 1:length(num_clusters)
                    
                    lst_B    = find( clusters == j );
                    
                    plv_cluster_map_data{ i,j } = reshape( plv_data_2d( lst_A, lst_B ), 1, length(lst_A)*length(lst_B) );
                    
                end
                
            end
            
        end
        
        % Subplot all the combinations of the plv_cluster_map data from loResPopSampleFromHiResCluster
        % VERY SLOW FOR HIGH # OF CLUSTERS
        function plot_cell_cluster_histograms( plv_cluster_map_data )
            
            [r,c] = size( plv_cluster_map_data );
            
            count = 1;
            for i = 1:r
                for j =1:c
                    if( j > i ); count = count + 1; continue; end;
                    subplot( r,c, count );
                    hist( plv_cluster_map_data{i,j} );
                    count = count + 1;
                end
            end
            
        end
        
        % Go from data at cluster locations to optColin values for display on optColin surface
        function vd = upsampleClusterIdxs( cluster_idxs, cluster_values )
            
            u = unique( cluster_idxs );
            
            vd = zeros( 1,length(cluster_idxs) );
            
            for i = 1:length(u)
                
                lst = find( cluster_idxs == u(i) );
                vd(lst) = cluster_values(i);
                
            end
            
            
        end
        
        
        % --------------------------------
        % Overlay Spheres On Brain
        % --------------------------------
        
        function overlayGraph( adjmatrix )
            
            adjmatrix;
            
            vd = zeros(1,80);
            vd(8) = 1;
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( betaH.idx, vd );
            lst = find(i>0);
            v = optColin.inflated.surface.L.v(lst,:);
            mv = mean(v);
            
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( betaH.idx, vd );
            
        end
        
        function mv = getCentroid( clusterIdxs, idx, optColin )
            
            vd = zeros(1,80);
            vd(idx) = 1;
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( clusterIdxs, vd );
            lst = find(i>0);
            v = optColin.inflated.surface.L.v(lst,:);
            mv = mean(v);
            
        end
        
        function mv = getCentroidR( clusterIdxs, idx, optColin )
            
            vd = zeros(1,80);
            vd(idx) = 1;
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( clusterIdxs, vd );
            lst = find(i>0) - length(optColin.inflated.surface.L.v);
            v = optColin.inflated.surface.R.v(lst,:);
            mv = mean(v);
            
        end
        
        function normal = getNormal( clusterIdxs, idx, optColin )
            
            vd = zeros(1,80);
            vd(idx) = 1;
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( clusterIdxs, vd );
            lst = find(i>0);
            v = optColin.inflated.surface.L.v(lst,:);
            
            vn = get(optColin.inflated.handle.L, 'VertexNormals' );
            normal = mean(vn(lst,:));
        end
        
        function normal = getNormalR( clusterIdxs, idx, optColin )
            
            vd = zeros(1,80);
            vd(idx) = 1;
            i = pitt.cerebro.DispBrain.upsampleClusterIdxs( clusterIdxs, vd );
            lst = find(i>0) - length(optColin.inflated.surface.L.v);
            v = optColin.inflated.surface.R.v(lst,:);
            
            vn = get(optColin.inflated.handle.R, 'VertexNormals' );
            normal = mean(vn(lst,:));
        end
        
        function drawSphere( centroidCoords )
            
            [x,y,z] = sphere(100);
            
            factor = 100;
            
            x = x./factor + centroidCoords(1);
            y = y./factor + centroidCoords(2);
            z = z./factor + centroidCoords(3);
            
            %x = (x + centroidCoords(1))./factor;
            %y = (y + centroidCoords(2))./factor;
            %z = (z + centroidCoords(3))./factor;
            
            
            hold on;
            s = surf(x,y,z);
            
            set(s,'edgecolor','r');
            
        end
        
        function s = drawSphereL( centroidCoords, relativeSize )
            
            if( nargin == 1 )
                relativeSize = 1;
            end
            
            [x,y,z] = sphere(100);
            
            factor = 110-(relativeSize.*40);
            
            x = x./factor + -.12;
            y = y./factor + centroidCoords(2);
            z = z./factor + centroidCoords(3);
            
            %x = (x + centroidCoords(1))./factor;
            %y = (y + centroidCoords(2))./factor;
            %z = (z + centroidCoords(3))./factor;
            
            
            hold on;
            s = surf(x,y,z);
            
            set(s,'edgecolor','r');
            
            %text(x(1),y(1),z(1),'here');
        end
        
        function s = drawSphereR( centroidCoords, relativeSize )
            
            if( nargin == 1 )
                relativeSize = 1;
            end
            
            [x,y,z] = sphere(100);
            
            factor = 110-(relativeSize.*40);
            
            x = x./factor + .12;
            y = y./factor + centroidCoords(2);
            z = z./factor + centroidCoords(3);
            
            %x = (x + centroidCoords(1))./factor;
            %y = (y + centroidCoords(2))./factor;
            %z = (z + centroidCoords(3))./factor;
            
            
            hold on;
            s = surf(x,y,z);
            
            set(s,'edgecolor', 'r');
            %set(s,'FaceColor','interp','EdgeColor','none')
            
            %text(x(1),y(1),z(1),'here');
        end
        
        
        function s = drawSphereIL( centroidCoords, relativeSize )
            
            if( nargin == 1 )
                relativeSize = 1;
            end
            
            centroidCoords = pitt.cerebro.DispBrain.LHtransform( centroidCoords );
            
            [x,y,z] = sphere(100);
            
            factor = 110-(relativeSize.*40);
            
            x = x./factor + centroidCoords(1);
            y = y./factor + centroidCoords(2) - .2;
            z = z./factor + centroidCoords(3);
            
            
            
            %x = (x + centroidCoords(1))./factor;
            %y = (y + centroidCoords(2))./factor;
            %z = (z + centroidCoords(3))./factor;
            
            
            hold on;
            s = surf(x,y,z);
            
            set(s,'edgecolor', 'r');
            %set(s,'FaceColor','interp','EdgeColor','none')
            
            %text(x(1),y(1),z(1),'here');
        end
        
        function s = drawSphereIR( centroidCoords, relativeSize )
            
            if( nargin == 1 )
                relativeSize = 1;
            end
            
            centroidCoords = pitt.cerebro.DispBrain.RHtransform( centroidCoords );
            
            [x,y,z] = sphere(100);
            
            factor = 110-(relativeSize.*40);
            
            x = x./factor + centroidCoords(1);
            y = y./factor + centroidCoords(2) - .2;
            z = z./factor + centroidCoords(3);
            
            
            
            %x = (x + centroidCoords(1))./factor;
            %y = (y + centroidCoords(2))./factor;
            %z = (z + centroidCoords(3))./factor;
            
            
            hold on;
            s = surf(x,y,z);
            
            set(s,'edgecolor', 'r');
            %set(s,'FaceColor','interp','EdgeColor','none')
            
            %text(x(1),y(1),z(1),'here');
        end
        
        
        function h = drawLineL( coords1, coords2 )
            
            hold on;
            h = line( linspace(-.12,-.12,10000), linspace(coords1(2),coords2(2),10000), linspace(coords1(3),coords2(3),10000) );
            
        end
        
        function h = drawLineR( coords1, coords2 )
            
            hold on;
            h = line( linspace(.12,.12,10000), linspace(coords1(2),coords2(2),10000), linspace(coords1(3),coords2(3),10000) );
            
        end
        
        % RH / LH Point
        function h = drawLineI( coords1, coords2 )
            
            coords1 = pitt.cerebro.DispBrain.RHtransform( coords1 );
            coords2 = pitt.cerebro.DispBrain.LHtransform( coords2 );
            
            coords2(2) = coords2(2) - .2;
            coords1(2) = coords1(2) - .2;
            
            hold on;
            h = line( linspace(coords1(1),coords2(1),10000), linspace(coords1(2),coords2(2),10000), linspace(coords1(3),coords2(3),10000) );
            
        end
        
        
        function adjmatrix = processAdjMatrixFull( adjmatrix, idxs, optColin )
            
            [am, Lidx]  = pitt.cerebro.DispBrain.processAdjMatrixL( adjmatrix(1:40,1:40), idxs, optColin );
            [am, Ridx]  = pitt.cerebro.DispBrain.processAdjMatrixR( adjmatrix(41:end,41:end), idxs, optColin );
            
            for i=Lidx
                adjmatrix( i,: ) = 0;
                adjmatrix( :,i ) = 0;
            end
            
            
            for i=Ridx
                adjmatrix( i,: ) = 0;
                adjmatrix( :,i ) = 0;
            end
            
        end
        
        function [adjmatrix, ridx] = processAdjMatrixL( adjmatrix, idxs, optColin )
            
            ridx = [];
            
            for i = 1:size(adjmatrix,1)
                
                nv = pitt.cerebro.DispBrain.getNormal( idxs, i, optColin );
                
                nonv = nv./norm(nv);
                if( nonv(1) <= -.3 )
                    %fprintf('removing %i\n', i );
                    adjmatrix( i,: ) = 0;
                    adjmatrix( :,i ) = 0;
                    ridx = [ridx, i];
                end
                
            end

            
        end
        
        function [adjmatrix, ridx] = processAdjMatrixR( adjmatrix, idxs, optColin )
            
            ridx = [];
            
            for i = 40+1:40+size(adjmatrix,1)
                
                nv = pitt.cerebro.DispBrain.getNormalR( idxs, i, optColin );
                
                nonv = nv./norm(nv);
                if( nonv(1) >= .3 )
                    %fprintf('removing %i\n', i );
                    adjmatrix( i-40,: ) = 0;
                    adjmatrix( :,i-40 ) = 0;
                    ridx = [ridx, i];
                end
                
            end

            
        end
        
        function adjmatrix = processAdjMatrixI( adjmatrix, idxs, optColin )
            
            L = ones( 40 );
            Lt = pitt.cerebro.DispBrain.processAdjMatrixL( L, idxs, optColin );
            R = ones( 40 );
            Rt = pitt.cerebro.DispBrain.processAdjMatrixL( R, idxs, optColin );
            
            lst = find( sum(Lt) == 0 );
            adjmatrix( :,lst ) = 0;
            
            lst = find( sum(Rt) == 0 );
            adjmatrix( lst,: ) = 0;
            
            
            
        end
        
        
        
        function drawAdjMatrixAsGraphL( adjmatrix, idxs, optColin )
            
            verts = max(adjmatrix);
            lstverts = find(verts>0);
            
            %set(get(optColin.inflated.handle.L,'parent'), 'currentaxes', optColin.inflated.handle.L);
            %set(0, 'currenfigure', get(optColin.inflated.handle.L,'parent'));
            for i=lstverts
                
                mv = pitt.cerebro.DispBrain.getCentroid( idxs, i, optColin );
                s = pitt.cerebro.DispBrain.drawSphereL( mv, sum(adjmatrix(i,:))./max(sum(adjmatrix)) );
                %set(s,'edgecolor','r');
                
            end
            
            %shading interp;
            
            for i=1:size(adjmatrix,1)
                for j=1:size(adjmatrix,2)
                    
                    if( j > i ); continue; end;
                    
                    if( adjmatrix(i,j) > 0 )
                        mvI = pitt.cerebro.DispBrain.getCentroid( idxs, i, optColin );
                        mvJ = pitt.cerebro.DispBrain.getCentroid( idxs, j, optColin );
                        L = pitt.cerebro.DispBrain.drawLineL( mvI, mvJ );
                        set(L,'linewidth', .5);
                        set(L,'Color','k');
                    end
                    
                end
            end
            
        end
        
        function drawAdjMatrixAsGraphR( adjmatrix, idxs, optColin )
            
            verts = max(adjmatrix);
            lstverts = find(verts>0);
            
            %set(get(optColin.inflated.handle.L,'parent'), 'currentaxes', optColin.inflated.handle.L);
            %set(0, 'currenfigure', get(optColin.inflated.handle.L,'parent'));
            for i=lstverts
                i = i+40;
                mv = pitt.cerebro.DispBrain.getCentroidR( idxs, i, optColin );
                s = pitt.cerebro.DispBrain.drawSphereR( mv, sum(adjmatrix(i-40,:))./max(sum(adjmatrix)) );
                %set(s,'edgecolor','r');
                
            end
            
            %shading interp;
            
            for i=1:size(adjmatrix,1)
                
                for j=1:size(adjmatrix,2)
                    
                    if( j > i ); continue; end;
                    
                    if( adjmatrix(i,j) > 0 )
                        mvI = pitt.cerebro.DispBrain.getCentroidR( idxs, i+40, optColin );
                        mvJ = pitt.cerebro.DispBrain.getCentroidR( idxs, j+40, optColin );
                        L = pitt.cerebro.DispBrain.drawLineR( mvI, mvJ );
                        set(L,'linewidth', .5);
                        set(L,'Color','k');
                    end
                    
                end
            end
            
        end
        
        function drawAdjMatrixAsGraphI( adjmatrix, idxs, optColin )
            
            verts = max(adjmatrix);
            lhtverts = find(sum(adjmatrix)>0);
            rhtverts = find(sum(adjmatrix,2)>0)';
            
            %set(get(optColin.inflated.handle.L,'parent'), 'currentaxes', optColin.inflated.handle.L);
            %set(0, 'currenfigure', get(optColin.inflated.handle.L,'parent'));
            for i=lhtverts
                mv = pitt.cerebro.DispBrain.getCentroid( idxs, i, optColin );
                s = pitt.cerebro.DispBrain.drawSphereIL( mv, sum(adjmatrix(:,i))./max(sum(adjmatrix)) );

                %set(s,'edgecolor','r');
            end
            
            for i=rhtverts
                mv = pitt.cerebro.DispBrain.getCentroidR( idxs, i+40, optColin );
                s = pitt.cerebro.DispBrain.drawSphereIR( mv, sum(adjmatrix(i,:))./max(sum(adjmatrix)) );

                %set(s,'edgecolor','r');
            end
            
            %shading interp;
            
            for i=1:size(adjmatrix,1)
                % Row - RH
                
                for j=1:size(adjmatrix,2)
                    % Col - LH
                    
                    rhidx = i + 40;
                    lhidx = j;
                    
                    %if( j > i ); continue; end;
                    
                    if( adjmatrix(i,j) > 0 )
                        mvI = pitt.cerebro.DispBrain.getCentroidR( idxs, rhidx, optColin );
                        mvJ = pitt.cerebro.DispBrain.getCentroid( idxs, lhidx, optColin );
                        L = pitt.cerebro.DispBrain.drawLineI( mvI, mvJ );
                        set(L,'linewidth', .5);
                        set(L,'Color','k');
                    end
                    
                end
            end
            
        end
        
        % Poor attempt
        function adjmatrix = removeInteriorClusters( adjmatrix, idxs, optColin )
            
            for i=1:size(adjmatrix)
                mv = pitt.cerebro.DispBrain.getCentroid( idxs, i, optColin );
                if( mv(1) < -0.038 )
                    adjmatrix(i,:) = 0;
                    adjmatrix(:,i) = 0;
                end
            end
            
        end
        
        
    end
    
end

