classdef AnalyzePLVMap < handle
    
    methods( Static )
        
        %{
        Init Enviornment:
        
        [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
        [A, B, C] =  pitt.diss.cluster.AnalyzePLVMap.loadVertexCoords();
        h = pitt.exp.simu.GraphAnalysis.displayBrain( v, f ); view( -90, 0 );
        [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
        
        
        Load AB data:
        freq6rest = pitt.diss.PLVAnalysis.read_freq_mat_half_file_real( '/synapse/logs/schmidtb/plvdata/opt001/empty_room_plv/freq_006' );
        freq6rest = freq6rest( 1:4098, 1:4098 );
        freq15rest = pitt.diss.PLVAnalysis.read_freq_mat_half_file_real( '/synapse/logs/schmidtb/plvdata/opt001/empty_room_plv/freq_015' );
        freq15rest = freq15rest( 1:4098, 1:4098 );
        
        freq6ab = pitt.diss.cluster.AnalyzePLVMap.load_simu_freq_data( '/synapse/logs/schmidtb/plvsim/point_sim/AB_n1_sim01_plv/freq_006' );
        freq15ab = pitt.diss.cluster.AnalyzePLVMap.load_simu_freq_data( '/synapse/logs/schmidtb/plvsim/point_sim/AB_n1_sim01_plv/freq_015' );
        
        %}
        
        
        
        % Given plv data and coordinates, identify k_cluster number of clusters using kmeans
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
        
        % Calculate T-Stat for each cluster compared to empty room
        function cluster_t_values = t_test_cluster_regions( plv_data_1d, lh_coords, cluster_vertex_idxs, k_clusters, plv_data_empty_2d )
            
            % Calculate for each cluster
            for i = 1:k_clusters
                
                lst = find( cluster_vertex_idxs == i );
                
                data_pop = plv_data_1d(lst);
                empty_pop = [];
                
                for j = 1:length(lst)
                    empty_pop = [empty_pop; plv_data_empty_2d(j,lst)'];
                end
                
                [H,P] = ttest2( data_pop, empty_pop );
                
                if( H==1 )
                    cluster_t_values(lst) = i.*(1-P);
                else
                    % Zero for non-sig
                    cluster_t_values(lst) = H;
                end
                
                %cluster_t_values(lst) = length(lst);
            end
            
        end
        
        
        function plv_data_2d = process_plv_data( plv_data_2d, distances )
            
            kernel1 = legion.Kernel();
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', 0.2, 1.1 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel1.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', .15, 10.282 );
            kernel1.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            tic;
            for i = 1:size(plv_data_2d,1)
                if( mod(i,100) == 0); fprintf('%i\n',i); end;
                kernel1.initial( plv_data_2d(i,:) );
                plv_data_2d(i,:) = kernel1.execute();
            end
            toc;
        end
        
        function sort_mask = sort_mask( distances, radius )
            
            sort_mask = zeros(size(distances,1), size(distances,2));
            %w = waitbar(0,'');
            for i = 1:size(distances,1)
                %waitbar(i/size(distances,1),w, num2str(i));
                [S,I] = sort(distances(i,:));
                lst = I(1:radius);
                sort_mask(i,lst) = 1;
                
            end
            %close(w);
        end
        
        
        function z = fisherz(r);

            % Z=fisherz(R)
            %
            % this function performs a fisher Z-transform of vector R and output the
            % transformed vector Z. This function is used to modify for example
            % correlations in [-1;+1] into a near gaussian population. It could be used
            % to transform any uniform population in the interval [-1;+1] into a
            % gaussian population. The inverse operation is done by function inverse_fisherz.m
            %
            % Vincent MORON
            % feb. 2005

            z=0.5*log((1+r)./(1-r));

        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        % freq15 = pitt.diss.cluster.AnalyzePLVMap.load_simu_freq_data( 'freq_015' );
        function freq = load_simu_freq_data( simu_data_path )
            
            freq = pitt.diss.PLVAnalysis.read_freq_mat_half_file( simu_data_path );
            
        end
        
        % [A, B, C] =  pitt.diss.cluster.AnalyzePLVMap.loadVertexCoords();
        function [A, B, C] = loadVertexCoords()
           A = 3673;
           B = 534;
           C = 2029;
        end
        
        % kernel = pitt.diss.cluster.AnalyzePLVMap.thresh_smooth_thresh_smooth( distances );
        function kernel = thresh_smooth_thresh_smooth( distances )
            
            kernel = legion.Kernel();
            kernel.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', 0.2, 1.1 );
            kernel.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            kernel.add( @pitt.exp.simu.GraphAnalysis.boundPLV, 'X', .15, 10.282 );
            kernel.add( @pitt.diss.PLVAnalysis.average_smoothing_filter, 'X', distances, 30 );
            
            
        end
        
        function ind = morton(n)
            % MORTON(N) return the morton permutation order for array of size 2^N
            % e.g.
            % n=2; % matrix size is 2^n
            % ind=morton(n);
            % d=fix(rand(2^n,2^n)*10);
            % disp(d)
            % disp(d(ind))
            
            linind=(1:4^n)-1; %start index count for array at zero
            ind2str=dec2bin(linind); %convert indices to base-2
            rb=ind2str(:,1:2:end); %take alternating bits for row and column
            cb=ind2str(:,2:2:end);
            r=bin2dec(rb)+1; %convert the row from bit to decimal
            c=bin2dec(cb)+1; %convert column
            ind=[2^n*(c-1)+r]'; %make a linear index into array for easy addressing
        end
        
    end
    
end