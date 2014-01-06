classdef PLVAnalysis < handle
    
    methods (Static)
        %{ 
        For initializaiton use:
        
        [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
        
        %}
        
        
        %{ if it looks like it's not right, you might be viewing too much data. Look closer }%
        function full_plv_map = full_map( lower_tri_plv_map )
            
            for i = 1:size( lower_tri_plv_map, 1 )
                for j = 1:size( lower_tri_plv_map, 2 )
                    if( j > i ); lower_tri_plv_map(i,j) = 0; end;
                end
            end
            
            full_plv_map = lower_tri_plv_map + tril( lower_tri_plv_map, 1 )';
            lst = find( full_plv_map == 2 );
            full_plv_map( lst ) = 1;
        end
        
        
        function freq = read_freq_mat_half_file( mat_path )
            
            freq = load( mat_path );
            freq = freq.freq;
            freq = pitt.diss.PLVAnalysis.full_map( freq );
            
        end
        
        function freq = read_freq_mat_half_file_real( mat_path )
            
            freq = load( mat_path );
            freq = freq.output;
            freq = pitt.diss.PLVAnalysis.full_map( freq );
            
        end
        
        %{ 

        plot vertex indexes on the LH 
        
        figure; pitt.diss.PLVAnalysis.identify_vertices_lh( lh_coords, [3670, 560] );
        
        %}
        function identify_vertices_lh( coords_lh, vertices_lh )
            subplot( 2,3,[ 1 2 3] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3) );
            hold on;
            scatter3( coords_lh(vertices_lh,1), coords_lh(vertices_lh,2), coords_lh(vertices_lh,3), 300.*ones(size(vertices_lh,1),1), 'filled' );
            hold off;
            view( -90, 0 );
            
            subplot( 2,3,[4 5 6] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3) );
            hold on;
            scatter3( coords_lh(vertices_lh,1), coords_lh(vertices_lh,2), coords_lh(vertices_lh,3), 300.*ones(size(vertices_lh,1),1), 'filled' );
            hold off;
            view( -90, 90 );
            
        end
        
        
        function overlay_lh_values( coords_lh, values, seed_idx )
           
            if( size( values,1 ) ~= size( coords_lh, 1 ) )
                error( sprintf( 'Values must be %i x 1', size( coords_lh,1 ) ) );
            end
            
            if( nargin == 2 )
                do_seed = 0;
            else
                do_seed = 1;
            end
            
            lst = find( values <= 0 );
            if( length( lst ) > 0 )
                fprintf( '\n*** Zero Valued Entries encountered, fixing... ***\n\n' );
                values( lst) = .0001;
            end
            
            subplot( 2,3,[ 1 2 3] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3), values, values );
            if( do_seed )
                hold on;
                scatter3( coords_lh(seed_idx,1), coords_lh(seed_idx,2), coords_lh(seed_idx,3), 300.*ones(size(seed_idx,2),1), 'filled' );
                hold off;
            end
            view( -90, 0 );
            
            subplot( 2,3,[4 5 6] );
            scatter3( coords_lh(:,1), coords_lh(:,2), coords_lh(:,3), values, values );
            if( do_seed )
                hold on;
                scatter3( coords_lh(seed_idx,1), coords_lh(seed_idx,2), coords_lh(seed_idx,3), 300.*ones(size(seed_idx,2),1), 'filled' );
                hold off;
            end
            view( -90, 90 );
            
            
            
        end
        
        function false_positive_func_distance( distances, values )
            
            lst = find( values >  0 );
            size(lst)
            
            % Max distance:
            max_dist = 169.4076;
            
            scatter( distances, values );
            
            if( max( values ) > 1 )
                axis( [0 max_dist 0 max( values )] );
            else
                axis( [0 max_dist 0 1] )
            end
            
            xlabel( 'Distance' );
            ylabel( 'PLV Values' );
            title( 'PLV as a function of distance' );
            
        end
        
        %{
        
        Anfreq15 = load( 'raw_A_only_noise_plv/freq_015' );
        Anfreq6 = load( 'raw_A_only_noise_plv/freq_006' );
        
        freq    = Anfreq15.freq;
        Anfreq15 = pitt.diss.PLVAnalysis.full_map( freq );
        freq    = Anfreq6.freq;
        Anfreq6 = pitt.diss.PLVAnalysis.full_map( freq );
        
        Bnfreq15 = load( 'raw_B_only_noise_plv/freq_015' );
        Bnfreq6 = load( 'raw_B_only_noise_plv/freq_006' );
        
        freq    = Bnfreq15.freq;
        Bnfreq15 = pitt.diss.PLVAnalysis.full_map( freq );
        freq    = Bnfreq6.freq;
        Bnfreq6 = pitt.diss.PLVAnalysis.full_map( freq );
        
        
        ABnfreq15 = load( 'raw_A_B_combo_noise_plv/freq_015' );
        ABnfreq6 = load( 'raw_A_B_combo_noise_plv/freq_006' );
        
        freq    = ABnfreq15.freq;
        ABnfreq15 = pitt.diss.PLVAnalysis.full_map( freq );
        freq    = ABnfreq6.freq;
        ABnfreq6 = pitt.diss.PLVAnalysis.full_map( freq );
        
        
        % Check status of running plv job
        legion.stream.Master.status( 2, 100, '/synapse/logs/schmidtb/plvsim/ABC/raw_An1_only_noise_plv/freq_015/' )
        
        % Run to process after PLV calculation from wihtin directory
        An1_freq015 = legion.stream.Util.read_2d_split_as_array( './freq_015', 100 ); 
        An1_freq006 = legion.stream.Util.read_2d_split_as_array( './freq_006', 100 ); 
        freq = An1_freq006; 
        save( 'freq_006.mat', 'freq' );
        freq = An1_freq015;
        save( 'freq_015.mat', 'freq' );
        An1freq15 = An1_freq015;
        An1freq15 = pitt.diss.PLVAnalysis.full_map( An1freq15 );
        An1freq6 = An1_freq006;
        An1freq6 = pitt.diss.PLVAnalysis.full_map( An1freq6 );
        
        
        Ac_freq015 = legion.stream.Util.read_2d_split_as_array( './freq_015', 100 ); 
        Ac_freq006 = legion.stream.Util.read_2d_split_as_array( './freq_006', 100 ); 
        freq = Ac_freq006; 
        save( 'freq_006.mat', 'freq' );
        freq = Ac_freq015;
        save( 'freq_015.mat', 'freq' );
        Acfreq15 = Ac_freq015;
        Acfreq15 = pitt.diss.PLVAnalysis.full_map( Acfreq15 );
        Acfreq6 = Ac_freq006;
        Acfreq6 = pitt.diss.PLVAnalysis.full_map( Acfreq6 );
        
        
        idx1 = 3670; idx2 = 560;
        
        freq = Acfreq15;
        idx  = idx2;
        opidx = idx2;
        
        %freq = ABnfreq15 - Anfreq15 - Bnfreq15;
        %freq = ABnfreq15 - Anfreq6;
        data = freq(:, idx );
        
        
        figure(1); pitt.diss.PLVAnalysis.overlay_lh_values( lh_coords, data.*10, [idx] )
        figure(2); pitt.diss.PLVAnalysis.false_positive_func_distance( distances( idx , : ), data )
        
        
        
        %figure(2); pitt.diss.PLVAnalysis.false_positive_func_distance( distances( opidx , : ), data )
        
        
        
        maxes = max( ABnfreq15 );
        
        scalAnfreq6 = Anfreq6 ./ repmat(max(Anfreq6 ), 4098,1);
        freq = ABnfreq15 - smooth(scalAnfreq6');
        
        
        smo = smooth(scalAnfreq6( idx, : ));
        scalsmo = smo ./ (max(smo));
        data = ABnfreq15(:,idx) - scalsmo;
        
        
        
        out = ABnfreq15(idx,:); 
        num_close = 100;
        [S,I] = sort( distances( idx,: ) );
        for i = 1:num_close; 
            out = out - ABnfreq6(I(i),:); 
        end;
        figure; plot( out );
        
        
        
        
        
        out = ABnfreq15(idx,:); 
        num_close = 100;
        [S,I] = sort( distances( idx,: ) );
        for i = 1:4098; 
            [Q,out] = deconv( out, ABnfreq6( I(i),: ) ); 
        
        end;
        figure; plot( out );
        
        
        figure; plot( Bnfreq15( idx,: ) );
        figure; plot( Anfreq15( idx,: ) );
        
        
        
        
        
        
        
        
        
        
        
        freq    = load( 'freq_006.mat' );
        freq    = freq.freq;
        freq006 = pitt.diss.PLVAnalysis.full_map( freq );
        
        freq    = load( 'freq_015.mat' );
        freq    = freq.freq;
        freq015 = pitt.diss.PLVAnalysis.full_map( freq );
        
        figure(5); imagesc( freq006 ); title( 'Freq 006 ');
        figure(5); imagesc( freq015 ); title( 'Freq 015 ');
        
        idx1 = 3670; idx2 = 560;
        
        freq = freq015;
        idx  = idx2;
        
        figure(1); pitt.diss.PLVAnalysis.overlay_lh_values( lh_coords, (freq(:, idx ).*10)+.001, idx )
        figure(2); pitt.diss.PLVAnalysis.false_positive_func_distance( distances( idx , : ), freq( idx , : ) )
        
        %}
        
        function total = density_vs_distance( d, values )
            
            [H, R] = hist( d, 100 );
            R = [0 R];
            
            count = 0;
            total = [];
            
            for i = 2:length( R )-1
                
                tot_dist = H(i);
                % Find index which coresond to these distances
                plv_vals = find( d >= R(i-1) & d < R(i) );
                count = count + length( plv_vals );
                
                lg = length( plv_vals );
                if( lg <= 0 )
                    ttmpt = 0;
                else
                    ttmpt = sum( values( plv_vals ) ) ./ length( plv_vals );
                end
                %ttmpt = sum( values( plv_vals ).*10 );
                total = [total; ttmpt];
            end
            
            %plot( R(2:end-1), total );
            
            count;
        end
        
        function agg_total = all_density_plots( distances, phase_map, seed_idx )
            
            agg_total = [];
            
            for i = 1:size( distances,1 )
                
                d       = distances( i, : ); 
                data    = phase_map( seed_idx, : );
                
                total = pitt.diss.PLVAnalysis.density_vs_distance( d, data );
                agg_total = [agg_total; total'];
            end
            
        end
        
        
        function choice = pick_bad_cortical_spots( coords_lh )
            
            choice = zeros( size(coords_lh, 1), 1 );
            
            for i = 1:length( choice )
                
                figure( 17 ); 
                pitt.diss.PLVAnalysis.identify_vertices_lh( coords_lh, i );
                choice( i ) = input('Good? (1/0): ' );
            end
            
            
        end
        
        function decimate_surface()
            
            
            
            
        end
        
        
        function smoothed_data = average_smoothing_filter( plv_data_1d, distances_2d, width )
            
            smoothed_data = zeros( length(plv_data_1d),1 );
            
            % Loop over each point in data vector
            for i = 1:length( plv_data_1d )
                
                % Sort by the current point and get all it's neighbors
                [S,I] = sort( distances_2d(i,:) );
                
                new_data = mean( plv_data_1d( I(1:width) ) );
                smoothed_data( i ) = new_data;
                
            end
            
            
            
        end
        
        function smoothed_data = average_std_filter( plv_data_1d, distances_2d, width )
            
            smoothed_data = zeros( length(plv_data_1d),1 );
            
            % Loop over each point in data vector
            for i = 1:length( plv_data_1d )
                
                % Sort by the current point and get all it's neighbors
                [S,I] = sort( distances_2d(i,:) );
                
                new_data = std( plv_data_1d( I(1:width) ) );
                smoothed_data( i ) = new_data;
                
            end
            
            
            
        end
        
        
        function [Q,R] = bounded_spatial_deconv_kernel( plv_data_1d, err_plv_data_1d, distances_2d, width, idx )
            
            % All points close to the point of interest
            [S,I] = sort( distances_2d(idx,:) );
            
            % Bounded plv values
            bound_data  = plv_data_1d( I(1:width) );
            bound_err   = err_plv_data_1d( I(1:width) );
            
            [Q,R] = deconv( bound_data, bound_err );
            
            
        end
        
        function [Q,R] = bounded_spatial_deconv_kernel_mse( plv_data_1d, err_plv_data_1d, distances_2d, width, idx )
            
            % All points close to the point of interest
            [S,I] = sort( distances_2d(idx,:) );
            
            % Bounded plv values
            bound_data  = plv_data_1d( I(1:width) );
            bound_err   = err_plv_data_1d( I(1:width) );
            
            [Q,R] = deconv( bound_data, bound_err );
            
            Q = sum((R).^2);
            
        end
        
        function deconv_data = bounded_spatial_deconv_plv( plv_data_2d, plv_empty_2d, distances, width )
            
            deconv_data = zeros(size(plv_data_2d,1),1);
            
            for i = 1:size(plv_data_2d,1)
                
                [Q,R] = pitt.diss.PLVAnalysis.bounded_spatial_deconv_kernel( plv_data_2d(i,:), plv_empty_2d(i,:), distances, width, i );
                deconv_data(i) = Q;
                
            end
            
        end
        
        function deconv_data = bounded_spatial_deconv_plv_mse( plv_data_2d, plv_empty_2d, distances, width )
            
            deconv_data = zeros(size(plv_data_2d,1),1);
            
            for i = 1:size(plv_data_2d,1)
                
                [Q,R] = pitt.diss.PLVAnalysis.bounded_spatial_deconv_kernel_mse( plv_data_2d(i,:), plv_empty_2d(i,:), distances, width, i );
                deconv_data(i) = Q;
                
            end
            
        end
        
        % Matlab sucks. Long lines get moved in comments sections. Boo.
        function donothing
            
            A = 3673;
            B = 534;
            C = 2029;
            d1 = pitt.diss.PLVAnalysis.average_smoothing_filter( pitt.exp.simu.GraphAnalysis.boundPLV( freq15( A,: ), 0.2,1.1), distances, 30 );
            d2 = pitt.diss.PLVAnalysis.average_smoothing_filter( pitt.exp.simu.GraphAnalysis.boundPLV( freq6( A,: ), 0.2,1.1), distances, 30 );
            d3 = pitt.diss.PLVAnalysis.average_smoothing_filter( pitt.exp.simu.GraphAnalysis.boundPLV( freq6( C,: ), 0.2,1.1), distances, 30 );
            d4 = pitt.diss.PLVAnalysis.average_smoothing_filter( pitt.exp.simu.GraphAnalysis.boundPLV( freq6( B,: ), 0.2,1.1), distances, 30 );
            
            
            kernel1 = legion.Kernel();
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', 0.2, 1.1 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', .15, 10.282 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel1.initial( freq15opt(B,:) );
            data = kernel1.execute();
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( h, data );
            
            %pitt.exp.simu.GraphAnalysis.displayVertexLocation( h, 2015 );
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( h, pitt.diss.PLVAnalysis.isolate_k_clusters( data, lh_coords, 20 ) );
            
            kernel1.initial( freq15empty(A,:) );
            
            
            kernel2 = legion.Kernel();
            kernel2.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', 0.2, 1.1 );
            kernel2.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel2.initial( freq15(A,:) );
            
            
            
            
            kernel1 = legion.Kernel();
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', 0.4, 1.1 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 40 );
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', .15, 10.282 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel1.initial( freq15(B,:) );
            data = kernel1.execute();
            
            
            
            d6 = pitt.diss.PLVAnalysis.bounded_spatial_deconv_plv_mse( kernel1.execute(), kernel2.execute(), distances, 30 );
            
            pitt.exp.simu.GraphAnalysis.overlayPLVData( h, d6 );
            
        end
        
        
        function cluster_vertex_idxs = isolate_k_clusters( plv_data_1d, lh_coords, k_clusters )
            
            if( size(plv_data_1d,1) < size(plv_data_1d,2) )
                plv_data_1d =plv_data_1d';
            end
            
            kdata   = [plv_data_1d, lh_coords];
            lst     = find( plv_data_1d > 0 );
            idx     = kmeans( kdata(lst,:), k_clusters );
            vd      = zeros( 4098,1 );
            vd(lst) = idx;
            
            cluster_vertex_idxs = vd;
        end
        
        function ave_cluster_value = average_k_means_cluster( plv_data_1d, lh_coords, k_clusters )
            
            cluster_idxs = pitt.diss.PLVAnalysis.isolate_k_clusters( plv_data_1d, lh_coords, k_clusters );
            
            ave_cluster_value = zeros( length(plv_data_1d), 1 );
            
            for i = 1:k_clusters
                
                lst = find( cluster_idxs == i );
                ave_cluster_value(lst) = mean( plv_data_1d(lst) );
                
            end
            
        end
        
        
        function t_test_values = t_test_k_means_cluster( plv_data_1d, plv_data_base_1d, lh_coords, k_clusters )
            
            cluster_idxs = pitt.diss.PLVAnalysis.isolate_k_clusters( plv_data_1d, lh_coords, k_clusters );
            
            t_test_values = zeros( length(plv_data_1d), 1 );
            
            fwer_alpha = .05/k_clusters;
            
            for i = 1:k_clusters
                
                lst = find( cluster_idxs == i );
                
                data_cluster = plv_data_1d(lst);
                base_cluster = plv_data_base_1d(lst);
                
                [H,P] = ttest2( data_cluster, base_cluster, fwer_alpha );
                
                if( H==1 )
                    t_test_values(lst) = 1-P;
                else
                    % Zero for non-sig
                    t_test_values(lst) = H;
                end
                
            end
            
        end
        
        
        function t_test_values = t_test_k_means_cluster_cross_wise( plv_data_1d, plv_data_base_1d, lh_coords, k_clusters )
            
            % Calculate T-tests from seed to all
            t_tests = pitt.diss.PLVAnalysis.t_test_k_means_cluster( plv_data_1d, plv_data_base_1d, lh_coords, k_clusters );
            
            
            
        end
        
        
        function t_test_mat = t_test_mat( plv_data_2d, plv_data_base_2d, lh_coords, k_clusters )
            
            t_test_mat = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            w = waitbar(0,'Init');
            for i = 1:size(plv_data_2d,1)
                waitbar(i/size(plv_data_2d,1),w, num2str(i));
                t_test_mat(i,:) = pitt.diss.PLVAnalysis.t_test_k_means_cluster( plv_data_2d(i,:), plv_data_base_2d(i,:), lh_coords, k_clusters );
            end
            close(w);
            
        end
        
        
        % Calculate PLV from idx_seed to all vertices in data
        % lh_meg_data = mne_read_stc_file('fwdsubj-sim-cortex-lh.stc');
        function plv_values = plv_at_points( data_2d, idx_seed, Freqs, Fs )
            
            plv_values = zeros( size(data_2d,1), length(Freqs) );
            
            for i = 1:size(data_2d,1)
                
                output = pitt.exp.plv.PLV.resting_plv_pairewise_multiple_freq( data_2d(idx_seed,:), data_2d(i,:), Fs, Freqs );
                
                for j = Freqs
                    plv_values( i, j ) = output{j};
                end
                
            end
            
            
        end
        
        
        
        function degraded_plv_map = degradePLVMap( cluster_idxs, plv_data_2d )
            
            num_clusters = length(unique(cluster_idxs));
            degraded_plv_map = zeros( num_clusters );
            
            for i = 1:num_clusters
                
                idxs_i = find( cluster_idxs == i );
                
                for j = 1:num_clusters
                    
                    idxs_j = find( cluster_idxs == j );
                    
                    if( i==j)
                        degraded_plv_map( i,j ) = mean(mean( [plv_data_2d(idxs_i,idxs_i)] ));
                    end
                    
                    degraded_plv_map( i,j ) = mean(mean( [plv_data_2d(idxs_i, idxs_j)] ));
                end
                
            end
            
        end
        
        function unified_degraded_plv_map = unifyDegradePLVMap( cluster_idxs, plv_data_2d )
            
            num_clusters = length(unique(cluster_idxs));
            unified_degraded_plv_map = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:num_clusters
                
                idxs_i = find( cluster_idxs == i );
                
                for j = 1:num_clusters
                    
                    idxs_j = find( cluster_idxs == j );
                    
                    if( i==j)
                        unified_degraded_plv_map( idxs_i,idxs_i ) = mean(mean( [plv_data_2d(idxs_i,idxs_i)] ));
                    end
                    
                    unified_degraded_plv_map( idxs_i,idxs_j ) = mean(mean( [plv_data_2d(idxs_i, idxs_j)] ));
                end
                
            end
            
        end
        
    end
    
    
end