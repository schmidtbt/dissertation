classdef Stats < handle

	methods (Static)
        
        
        % --------------------------------
        % PLV Normalization
        % --------------------------------
        
        % Accepted convention for normalizing rayleigh distributed data
        function norm_data = rayleighToNorm( data_1d )
            norm_data = sqrt(data_1d);
        end
        
        % Fisher's Transformation: http://en.wikipedia.org/wiki/Fisher_transformation
        % This will output a normalized mean zero system
        % Probably redundant with normalization of rayleigh (norm of norm)
        function mean_zero = boundZeroOneToMeanZero( data_1d )
            % A hack is included because atanh(1) = infinity
            % In addition, we mean normalize before transforming
            mean_zero = atanh( 2.*(data_1d-.5) .*.9999 );
        end
        
        % Convert matrix into 1D data
        function oneD = makeOneD( plv_map_2d )
            oneD = reshape( plv_map_2d, 1, size(plv_map_2d,1)*size(plv_map_2d,2));
        end
        
        
        
        % --------------------------------
        % Permutation Testing
        % --------------------------------
        
        % Calculate the permutation test statistics between each cell array of rest data dn empty room data.
        % See Oostenveld 2007
        function [tstats, output, tpops] = permute_calc_stats( cell_arr_cluster_data_rest, cell_arr_cluster_data_er, num_sim )
            
            if( nargin == 2 )
                num_sim = 1;
            end
            
            
            num_rest        = length( cell_arr_cluster_data_rest );
            num_er          = length( cell_arr_cluster_data_er );
            num_clusters    = size( cell_arr_cluster_data_rest{1},1 );
            
            
            tpops   = zeros( num_clusters, num_clusters, num_sim );
                
            tstats  = zeros( num_clusters );
            output  = zeros( num_clusters );
            
            for i = 1:num_clusters
                for j = 1:num_clusters
                    
                    if( j > i ); continue; end; %skip other diagonal
                    %{
                    dataRest = cell(0);
                    for k = 1:num_rest
                        dataRest{k} = cell_arr_cluster_data_rest{k}{i,j};
                    end
                    
                    dataER = cell(0);
                    for k = 1:num_er
                        dataER{k} = cell_arr_cluster_data_er{k}{i,j};
                    end
                    %}
                    
                    allData = cell(0);
                    for k = 1:num_rest
                        allData{k} = cell_arr_cluster_data_rest{k}{i,j};
                    end
                    
                    for k = 1:num_er
                        allData{k+num_rest} = cell_arr_cluster_data_er{k}{i,j};
                    end
                    
                    %[H,P,Ci,stat] = ttest2( dataRest{1}, dataER{1} );
                    
                    % Calculate population statistics for this pairing
                    [tpop, tvalue] = pitt.cerebro.Stats.singlePermute( allData, num_rest, num_sim );
                    
                    % Calculate critical value for this population of statistics
                    crit_t_value = pitt.cerebro.Stats.tstatCriticalValue( tpop, .05 );
                    
                    tpops( i,j, : ) = tpop';
                    
                    % Display significant observed groupings
                    if( tvalue > crit_t_value )
                        output(i,j) = 1;
                    end
                    
                    % Report back all tstats of observed groupings
                    tstats(i,j) = [tvalue ];
                    
                end
            end
            
        end
        
        % Run a single cluster -> cluster permutation test
        function [tpop, tvalue] = singlePermute( allData, num_rest, num_sim )
            
            tot_length  = length(allData);
            tpop        = [];
            
            % Loop over number of simulations per cluster pair
            for i = 1:num_sim
                
                % Calculate random permutation order
                perm_order = randperm( tot_length );
                
                % Calculate statistics for this new grouping of data
                stat = pitt.cerebro.Stats.singleTest( allData, num_rest, perm_order );
                
                % Store population statistics of permutation
                tstat = stat.tstat;
                tpop  = [ tpop, tstat ];
            end
            
            % Calculate the observed statistics (not permuted)
            tvalue = pitt.cerebro.Stats.singleTest( allData, num_rest, [1:num_rest,num_rest+1:length(allData)] );
            tvalue = tvalue.tstat;
        end
        
        % Run a single t-test against two populations given a permutation order
        function stat = singleTest( allData, num_rest, perm_order )
                
            % We first resort data to be single giant matrices
            dRest = [];
            for i = perm_order( 1:num_rest )
                dRest = [dRest; allData{i}'];
            end
            
            dEmpty = [];
            for i = perm_order( num_rest+1:end )
                dEmpty = [dEmpty; allData{i}'];
            end

            % Now run stats against two populations
            [H,P,Ci,stat] = ttest2( dRest, dEmpty );
                
        end
        
        
        
        
        % --------------------------------
        % Hemisphere Testing - Inter/Intra
        % --------------------------------
        
        function sigT = hemisphereTstats( tstats, critvalue, numLHclusters )
            
            numRHclusters = size(tstats,1) - numLHclusters;
            
            if( length(critvalue) == 0 )
                cL = critvalue;
                cR = critvalue;
                cI = critvalue;
            else
                cL = critvalue(1);
                cR = critvalue(2);
                cI = critvalue(3); 
            end
            
            intraLcrit = pitt.cerebro.Stats.tstatCriticalValue( tstats(1:numLHclusters,1:numLHclusters), cL );
            intraRcrit = pitt.cerebro.Stats.tstatCriticalValue( tstats(numLHclusters:end,numLHclusters:end), cR );
            intercrit  = pitt.cerebro.Stats.tstatCriticalValue( tstats(numLHclusters:end,1:numLHclusters), cI );
            
            sigT = zeros(size(tstats,1));
            
            lstL = find( tstats > intraLcrit );
            lstR = find( tstats > intraRcrit );
            lstI = find( tstats > intercrit );
            
            L = zeros(size(tstats,1));
            R = L;
            I = L;
            
            L(lstL) = 1;
            R(lstR) = 1;
            I(lstI) = 1;
            
            sigT = I;
            sigT( 1:numLHclusters,1:numLHclusters ) = L( 1:numLHclusters,1:numLHclusters );
            sigT( numLHclusters:end,numLHclusters:end ) = R( numLHclusters:end,numLHclusters:end );
            
            
        end
        
        
        function cvertex = clusterHemi( tstats, type, numLHclusters )
            
            cvertex = zeros(numLHclusters, size(tstats,1) );
            rclusterStar = size(tstats,1) - numLHclusters + 1;
            
            switch( type )
                case 'L'
                    cvertex(1:numLHclusters,1:numLHclusters)   = tstats( 1:numLHclusters,1:numLHclusters );
                case 'R'
                    cvertex(1:numLHclusters,rclusterStar:end) = tstats( rclusterStar:end,rclusterStar:end );
                case 'I'
                    cvertex(1:numLHclusters,1:numLHclusters)   = tstats( rclusterStar:end,1:numLHclusters );
                    cvertex(1:numLHclusters,rclusterStar:end) = tstats( 1:numLHclusters,rclusterStar:end );
                case 'LR'
                    cvL = pitt.cerebro.Stats.clusterHemi( tstats, 'L', numLHclusters );
                    cvR = pitt.cerebro.Stats.clusterHemi( tstats, 'R', numLHclusters );
                    cvertex = cvL + cvR;
            end
            
            
        end
        
        
        
        % --------------------------------
        % Utility
        % --------------------------------
        
        % From distribution of tstats, find crit-value cutoff (e.g. 5% tail = .05) expressed as tstat value for the
        % cutoff
        function cutoff = tstatCriticalValue( tstats, crit_value )
            
            [N,X]   = hist( pitt.cerebro.Stats.makeOneD(tstats), 100 );
            r       = cumsum( N );
            idx     = dsearchn( r', (1-crit_value)*max(r) );
            cutoff  = X(idx);
        end
        
        
	end

end
