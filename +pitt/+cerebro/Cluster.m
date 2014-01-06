classdef Cluster < handle
    
    methods( Static )
        
        function [cluster_vertex_idxs, cluster_centroids] = isolate_k_clusters( plv_data_1d, coords, k_clusters )
            
            if( size(plv_data_1d,1) < size(plv_data_1d,2) )
                plv_data_1d =plv_data_1d';
            end
            
            kdata   = [plv_data_1d, coords];
            lst     = find( plv_data_1d > 0 );
            opts    = statset( 'MaxIter', 10000 );
            [idx,C] = kmeans( kdata(lst,:), k_clusters, 'Options', opts );
            vd      = zeros( length(plv_data_1d),1 );
            vd(lst) = idx;
            
            cluster_vertex_idxs = vd;
            cluster_centroids = C(:,2:4);
        end
        
        function cluster_avgs = avgInCluster( plv_data_1d, clusters )
            
            % Enforce column vector
            if( size(plv_data_1d,1) < size(plv_data_1d,2) )
                plv_data_1d =plv_data_1d';
            end
            
            num_clusters = unique( clusters );
            cluster_avgs = zeros(1,length(plv_data_1d));
            
            for i=1:length(num_clusters)
                
                cluster_idxs = find( clusters == i );
                cluster_avgs( cluster_idxs ) = mean( plv_data_1d(cluster_idxs) );
                
            end
            
        end
        
        function bounded_plv_2d = boundPLVMap( plv_data_2d, lower, upper )
            
            if(nargin == 2)
                upper = 100000000;
            end
            
            lst = find( plv_data_2d > lower & plv_data_2d < upper );
            bounded_plv_2d = zeros(size(plv_data_2d,1), size(plv_data_2d,2));
            bounded_plv_2d(lst) = plv_data_2d(lst);
            
        end
        
        function bound_vector = boundVector( vector, lower, upper )
            
            lst = find( vector > lower & vector < upper );
            bound_vector = zeros(1, length(vector));
            bound_vector(lst) = vector(lst);
            
        end
        
        function mean_cluster_to_all = meanClusterToAll( plv_data_2d, clusters, cluster_idx )
            
            lst = find( clusters == cluster_idx );
            mean_cluster_to_all = mean( plv_data_2d(lst,:),1 );
            
        end
        
        function mean_clusters_from_cluster = meanClusterToClusters( plv_data_2d, clusters, cluster_idx )
            
            lst = find( clusters == cluster_idx );
            mean_cluster_to_all = mean( plv_data_2d(lst,:),1 );
            mean_clusters_from_cluster = pitt.cerebro.Cluster.avgInCluster( mean_cluster_to_all, clusters );
            
        end
        
        function displayCluster( hL, hR, clusters, idxs )
            lst = find( clusters == idxs );
            vd = zeros( 1,8196 );
            vd(lst) = 1;
            pitt.cerebro.DispBrain.dispVertices( hL, hR, lst );
        end
        
        function plv_cluster_map = reducePLVByMeanClusters( plv_data_2d, clusters )
            
            num_clusters = unique( clusters );
            plv_cluster_map = zeros( length(num_clusters) );
            
            for i = 1:length(num_clusters)
                
                lst = find( clusters == i );
                row_mean = mean( plv_data_2d(lst, : ), 1 );
                
                for j = 1:length(num_clusters)
                    
                    lstj = find( clusters == j );
                    plv_cluster_map( i,j ) = mean( row_mean( lstj ) );
                    
                end
                
            end
            
        end
        
        function plv_map_ttest = reducePLVByTTestCluster( plv_rest_data_2d, plv_er_data_2d, clusters, pvalue )
            
            num_clusters = unique( clusters );
            plv_map_ttest = zeros( length(num_clusters) );
            
            for i = 1:length(num_clusters)
                
                lsti = find( clusters == i );
                
                for j = 1:length(num_clusters)
                    
                    lstj = find( clusters == j );
                    
                    restData = plv_rest_data_2d( lsti, lstj );
                    emprData = plv_er_data_2d( lsti, lstj );
                    
                    restData = pitt.cerebro.Cluster.zScored( pitt.cerebro.Cluster.makeOneD(sqrt(restData)) );
                    emprData = pitt.cerebro.Cluster.zScored( pitt.cerebro.Cluster.makeOneD(sqrt(emprData)) );
                    
                    restData = restData - mean(emprData);
                    emprData = emprData - mean(emprData);
                    
                    [H,P] = ttest( restData, emprData, pvalue, 'both' );
                    
                    
                    [Ne Xe] = hist( emprData, 100 );
                    [Nr Xr] = hist( restData, 100 );
                    
                    figure(6); 
                    clf;
                    plot( Xe, Ne' );
                    hold on;
                    plot( Xr, Nr' );
                    hold off;
                    
                    
                    figure(8); hist( restData );
                    figure(9); hist( emprData );
                    
                    disp(H);
                    %disp(P);
                    pause;
                    
                    if( ~isnan(H) & H )
                        plv_map_ttest( i,j ) = H;
                    end
                    
                end
                
            end
            
        end
        
        function zscored = zScored( plv_data_1d )
            
            lst = find( plv_data_1d == 1 );
            plv_data_1d(lst) = 0.99;
            zscored = (plv_data_1d - .5 ).*2;
            
        end
        
        % Assume 8196 points
        % Expand num_clusters and their indexes into 1d plv data
        function plv_data_1d = expandClusterValues( plv_cluster_data_1d, clusters, num_values )
            
            if( nargin == 2 )
                num_values = 8196;
            end
            
            num_clusters = unique( clusters );
            plv_data_1d = zeros( 1,num_values );
            
            for i = 1:length(num_clusters)
                
                lst = find( i == clusters );
                plv_data_1d( lst ) = plv_cluster_data_1d( i );
                
            end
            
        end
        
        function networks = floodFill( plv_map )
            
            num_nodes = size(plv_map,1);
            networks = zeros(num_nodes);
            accessed = zeros(1,num_nodes);
            
            for i = 1:num_nodes
                
                accessed(i) = 1;
                
                for j = 1:num_nodes
                    
                    if( accessed(j) == 1 ); continue; end;
                    
                    if( plv_map(i,j) > 0 )
                        networks(i,j) = 1;
                        accessed(j) = 1;
                    end
                    
                    
                end
                
            end
            
            
        end
        
        function networks = floodFill2( plv_map )
            
            num_nodes = size(plv_map,1);
            networks = zeros(num_nodes);
            accessed = zeros(num_nodes);
            usednodes = zeros(1,num_nodes);
            
            for i = 1:num_nodes
                
                usednodes(i) = 1;
                
                lst = find( plv_map( i,: ) > 0 );
                
                if( empty(lst) ); continue; end;
                
                for j = 1:length(lst)
                    if( usednodes(j) ); continue; end;
                    
                end
                
            end
            
        end
        
        function involved_nodes = recursiveDepthSearch( plv_map, idx )
            
            new_involved_nodes = find( plv_map(idx,:) > 0 );
            
            
            for i=1:length(involved_nodes)
                depth_nodes = pitt.cerebro.Cluster.recursiveFindGraphs( plv_map, idx, access_history );
            end
            
            
        end
        
        function oneD = makeOneD( plv_map_2d )
            
            oneD = reshape( plv_map_2d, 1, size(plv_map_2d,1)*size(plv_map_2d,2));
            
        end
        
        
        
        
        % --------------------------------
        % Methods for use with Opt* structures of pitt.cerebro.Reg
        % --------------------------------
        
        function idx = kClusterOnSurface( opt, surf_name, num_clusters, extra_data )
            
            vL = opt.(surf_name).surface.L.v;
            vR = opt.(surf_name).surface.R.v;
            
            opts = statset('MaxIter', 10000 );
            
            % LEFT HEMISPHERE
            cluster_data = [vL];
            
            if( nargin == 4 )
                if( size( extra_data, 1 ) ~= size(vL,1)+size(vR,1) )
                    error('Extra data is wrong size' );
                end
                
                lst          = find( extra_data(1:length(vL),:) > 0 );
                cluster_data = [extra_data(1:length(vL),:), [vL]];
                cluster_data = cluster_data(lst,:);
                
                idxLT = kmeans( cluster_data, num_clusters, 'Distance', 'cosine', 'Options', opts );
                
                idxL = zeros(1,length(vL) );
                idxL(lst) = idxLT;
                
            else
                idxL = kmeans( cluster_data, num_clusters, 'Options', opts );
            end
            
            % RIGHT HEMISPHERE
            cluster_data = [vR];
            
            if( nargin == 4 )
                if( size( extra_data, 1 ) ~= size(vL,1)+size(vR,1) )
                    error('Extra data is wrong size' );
                end
                
                lst          = find( extra_data(length(vL)+1:end,:) > 0 );
                cluster_data = [extra_data(length(vL)+1:end,:), [vR]];
                cluster_data = cluster_data(lst,:);
                
                idxRT = kmeans( cluster_data, num_clusters, 'Distance', 'cosine', 'Options', opts );
                %link = linkage( cluster_data );
                %idxRT = cluster( link, 'cutoff', 1.2 );
                
                idxR = zeros(1,length(vR) );
                idxR(lst) = idxRT;
                
            else
                idxR = kmeans( cluster_data, num_clusters, 'Options', opts );
            end
            
            
            
            
            idx = [idxL'; idxR'+num_clusters];
        end
        
        function [cluster_means, cluster_mean_by_idx] = reduceClustersByMean( cluster_idxs, data )
            
            if( length(cluster_idxs) ~= length(data) )
                error('Size mismatch');
            end
            
            u = unique( cluster_idxs );
            cluster_means = zeros(1,length(cluster_idxs) );
            cluster_mean_by_idx = zeros(1,length(u));
            
            for i =1:length(u)
                
                lst = find( cluster_idxs == u(i) );
                cluster_means( lst ) = mean( data(lst) );
                cluster_mean_by_idx(i) = mean( data(lst) );
            end
            
        end
        
    end
    
end

