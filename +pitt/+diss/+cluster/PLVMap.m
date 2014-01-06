classdef PLVMap < handle
    
    
    
    methods (Static)
        
        
        % --------------------------------------------
        %   MASKING
        % --------------------------------------------
        
        % Threshold 2D PLV Map
        function bounded_plv_2d = boundPLVMap( plv_data_2d, lower, upper )
            
            lst = find( plv_data_2d > lower & plv_data_2d < upper );
            bounded_plv_2d = zeros(size(plv_data_2d,1), size(plv_data_2d,2));
            bounded_plv_2d(lst) = plv_data_2d(lst);
            
        end
        
        function masked_plv_map2d = maskPLVMap( plv_data_2d, mask )
            
            L = length( find( mask(1,:) > 0 ));
            
            masked_plv_map2d = plv_data_2d * mask ./L;
        end
        
        function hard_mask = hardMask( data_2d, mask_2d )
            
            lst = find( mask_2d > 0 );
            hard_mask = data_2d;
            hard_mask(lst) = 0;
            
        end
        
        function plv_1d = boundPLV( plv_data_1d, lower, upper )
            
            if( size(plv_data_1d,1) > size(plv_data_1d,2) )
                plv_data_1d = plv_data_1d';
            end
            
            lst         = find( plv_data_1d < lower & plv_data_1d >= upper );
            plv_1d      = plv_data_1d;
            plv_1d(lst) = 0;
            
        end
        
        % --------------------------------------------
        %   CLUSTERING
        % --------------------------------------------
        
        % Identify clusters in data and coordinate space
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
        
        % Mean unify given cluster areas
        function unified_degraded_plv_map = unifyMeanDegradePLVMap( cluster_idxs, plv_data_2d )
            
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
        
        % Max unify given cluster areas
        function unified_degraded_plv_map = unifyMaxDegradePLVMap( cluster_idxs, plv_data_2d )
            
            num_clusters = length(unique(cluster_idxs));
            unified_degraded_plv_map = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:num_clusters
                
                idxs_i = find( cluster_idxs == i );
                
                for j = 1:num_clusters
                    
                    idxs_j = find( cluster_idxs == j );
                    
                    if( i==j)
                        unified_degraded_plv_map( idxs_i,idxs_i ) = max(max( [plv_data_2d(idxs_i,idxs_i)] ));
                    end
                    
                    unified_degraded_plv_map( idxs_i,idxs_j ) = max(max( [plv_data_2d(idxs_i, idxs_j)] ));
                end
                
            end
            
        end
        
        function unif_t_test_degrade = unifyTTestDegrade( cluster_idxs, plv_data_2d, plv_empty_2d, alpha )
            
            if (nargin == 3 )
                alpha = .05;
            end
            
            num_clusters = length(unique(cluster_idxs));
            unif_t_test_degrade = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:num_clusters
                
                idxs_i = find( cluster_idxs == i );
                
                for j = 1:num_clusters
                    
                    idxs_j = find( cluster_idxs == j );
                    
                    if( i==j)
                        
                        data = [plv_data_2d(idxs_i,idxs_i)];
                        rest = [plv_empty_2d(idxs_i,idxs_i)];
                        
                        data = reshape( data, size(data,1)*size(data,2), 1 );
                        rest = reshape( rest, size(rest,1)*size(rest,2), 1 );
                        
                        [H,R] = ttest2( data, rest, alpha );
                        unif_t_test_degrade( idxs_i,idxs_i ) = H;
                    end
                    
                    data = [plv_data_2d(idxs_i,idxs_j)];
                    rest = [plv_empty_2d(idxs_i,idxs_j)];
                    
                    data = reshape( data, size(data,1)*size(data,2), 1 );
                    rest = reshape( rest, size(rest,1)*size(rest,2), 1 );
                        
                    [H,R] = ttest2( data, rest, alpha );
                    
                    unif_t_test_degrade( idxs_i,idxs_j ) = H;
                end
                
            end
            
            
        end
        
        function degrade_plv_clusters = degradePLVClustering( cluster_idxs, plv_data_2d )
            
            num_clusters = length(unique(cluster_idxs));
            degrade_plv_clusters = zeros( num_clusters );
            
            for i = 1:num_clusters
                
                idxs_i = find( cluster_idxs == i );
                
                for j = 1:num_clusters
                    
                    idxs_j = find( cluster_idxs == j );
                    
                    if( i==j)
                        degrade_plv_clusters( i,i ) = mean(mean( [plv_data_2d(idxs_i,idxs_i)] ));
                    end
                    
                    degrade_plv_clusters( i,j ) = mean(mean( [plv_data_2d(idxs_i, idxs_j)] ));
                end
                
            end
            
        end
        
        function degradePLV = enforceBiDirectClustersDegraded( cluster_idxs, plv_data_2d )
            
            degradePLV = pitt.diss.cluster.PLVMap.degradePLVClustering( cluster_idxs, plv_data_2d );
            
            for i = 1:size(degradePLV,1)
                for j = 1:size(degradePLV,1)
                    
                    if( j > i ); continue; end;
                    
                    if( degradePLV(i,j) ~= degradePLV(j,i) )
                        degradePLV(i,j) = 0;
                        degradePLV(j,i) = 0;
                    end
                    
                end
            end
            
        end
        
        function diag_enforced = enforceBiDirectCluster( cluster_idxs, plv_data_2d )
            
            degradePLV = pitt.diss.cluster.PLVMap.enforceBiDirectClustersDegraded( cluster_idxs, plv_data_2d );
            
            num_clusters = length(unique(cluster_idxs));
            diag_enforced = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:size(degradePLV,1)
                for j = 1:size(degradePLV,2)
                    lsti = find( cluster_idxs == i );
                    lstj = find( cluster_idxs == j );
                    diag_enforced( lsti, lstj ) = degradePLV( i,j );
                end
            end
            
            
        end
        
        
        
        function tri_calc = triagEnforce( cluster_idxs, plv_data_2d )
            
            num_clusters = length(unique(cluster_idxs));
            tri_calc = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:num_clusters
                lsti = find( cluster_idxs == i );
                for j = 1:num_clusters
                    lstj = find( cluster_idxs == j );
                    
                    ratio = mean(mean( plv_data_2d( lsti,lstj ) ) ) ./ mean(mean( [reshape( plv_data_2d(lsti,lsti), length(lsti).^2,1); reshape( plv_data_2d(lstj,lstj), length(lstj).^2,1)] ));
                    
                    tri_calc( lsti,lstj ) = ratio;
                    
                end
            end
            
            
        end
        
        
        
        % I don't think these work presently
        function cluster = recursiveSearchBackWards( plv_data_2d, distances, idx, max_cluster_radius, cutoff_distance )
            
            
            data = plv_data_2d(idx,:);
            dist = distances(idx,:);
            
            [S,I] = sort( dist );
            
            cluster = zeros( 1, length(data) );
            clcounter = 1;
            
            for i = length(data):-1:1
                
                if( S(i) < cutoff_distance )
                    continue;
                end
                
                if( cluster(I(i)) > 0 )
                    continue;
                end
                
                if( data(I(i)) > 0 )
                    cluster(I(i)) = clcounter;
                    
                    % Neighborhood search
                    [Sn,In] = sort( distances(I(i),:) );
                    
                    for j = 1:max_cluster_radius
                       
                        if( data(In(j)) > 0 )
                            cluster(In(j)) = clcounter;
                        end
                    end
                    
                end
                
                clcounter = clcounter + 1;
                
            end
            
            
        end
        
        function cluster_2d = recursiveSearchBackWards2d( plv_data_2d, distances, max_cluster_radius, cutoff_distance )
            
            cluster_2d = zeros( size(plv_data_2d,1), size(plv_data_2d,2) );
            
            for i = 1:size(plv_data_2d,1)
                cluster_2d(i,:) = pitt.diss.cluster.PLVMap.recursiveSearchBackWards( plv_data_2d, distances, i, max_cluster_radius, cutoff_distance );
            end
            
            
        end
        
    end
    
    
end