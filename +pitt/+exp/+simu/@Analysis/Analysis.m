classdef Analysis < handle
   
    methods( Static )
        
        
        %{
        
        Calculate for each vertex, the number of non-zero connections
        coming into that vertex. Additionally calculate the number of
        non-zero connections going out of this node
        
        Inflow is number causing that vertex
        Outflow is number that this vertex causes
        
        Columns cause rows
        
        Use thresholding on input to restrict gc values to those of
        interested (eg significant only according to prb)
        
        %}
        function [inflow, outflow] = flowAnalysisGC( gc )
            
            [rows, cols] = size( gc );
            
            % Initialize two empty output variables
            inflow = zeros( rows,1 );
            outflow = zeros( rows,1 );
            
            
            
            for i = 1:cols
               outflow( i ) = sum( find(gc( :,i ) > 0 ) ); 
            end
            
            for i = 1:rows
               inflow( i ) = sum( find( gc(i,:) > 0 ) ); 
            end
            
            % We traverse columns first
            %{
            for i = 1:cols
                
                for j = 1:rows
                    
                    % Skip self terms
                    if( i == j ) continue; end;
                    
                    if( gc(j,i) > 0 )
                        inflow( j )     = inflow( j ) + 1;
                        outflow( i )    = outflow( i ) + 1;
                    end
                    
                end
            end
            %}
            
        end
        
        %{
        
        Calculate for each vertex, the number of non-zero connections
        coming into that vertex. Additionally calculate the number of
        non-zero connections going out of this node
        
        Inflow is sum of causing that vertex
        Outflow is sum of that this vertex causes
        
        Columns cause rows
        
        Use thresholding on input to restrict gc values to those of
        interested (eg significant only according to prb)
        
        %}
        function [inflow, outflow] = flowMagnitudeAnalysisGC( gc )
            
            [rows, cols] = size( gc );
            
            % Initialize two empty output variables
            inflow = zeros( rows,1 );
            outflow = zeros( rows,1 );
            
            for i = 1:cols
               outflow( i ) = sum( gc( :,i ) ); 
            end
            
            for i = 1:rows
               inflow( i ) = sum( gc(i,:) ); 
            end
            
            % We traverse columns first
            %{
            for i = 1:cols
                
                
                for j = 1:rows
                    
                    % Skip self terms
                    if( i == j ) continue; end;
                    
                    if( gc(j,i) > 0 )
                        inflow( j )     = inflow( j ) + gc(j,i);
                        outflow( i )    = outflow( i ) + gc(j,i);
                    end
                    
                end
            end
            %}
            
        end
        
        
        function [inflow, outflow, count] = flowMagnitudeByLabel( gc )
            
            [rows, cols] = size( gc );
            
            % Initialize two empty output variables
            inflow = zeros( 3,3 );
            outflow = zeros( 3,3 );
            count = zeros(3,3);
            
            
            % We traverse columns first
            for i = 1:cols
                
                label_col = pitt.exp.simu.Analysis.chooseLabel( 130, 192, 384, i );
                
                for j = 1:rows
                    
                    label_row = pitt.exp.simu.Analysis.chooseLabel( 130, 192, 384, j );
                    
                    % Skip self terms
                    if( i == j ) continue; end;
                    
                    count(label_row, label_col ) = count(label_row, label_col ) + 1;
                    
                    if( gc(j,i) > 0 )
                        inflow( label_row, label_col )     = inflow( label_row, label_col ) + gc(j,i);
                        outflow( label_row,label_col )    = outflow( label_row,label_col ) + gc(j,i);
                        
                    end
                    
                end
            end
            
            % Normalize by number of parts
            inflow = inflow ./count;
            outflow = outflow ./count;
            
        end
        
        function label = chooseLabel( endL1, endL2, endDat, idx )
            
            if( idx <= endL1 )
                label = 1;
            elseif( idx<= endL2 )
                label = 2;
            elseif( idx <= endDat )
                label = 3;
            end
        end
        
        
        %{
        Convert the flow numbers to per label normalized by number of
        vertices in each channel type (L1, L2, rand)
        %}
        function fa = flowLabelConversion( flow )
            fa = [sum( flow(1:130) ), sum( flow(131:192) ), sum( flow(193:end) )]./[130, 62, 384-192];
        end
        
        %{
        Remove any terms relating to same L1->L1 or L2->L2
        %}
        function gcst = removeSelfTerms( gcst )
           
            gcst( 1:130, 1:130 ) = 0;
            gcst( 131:192, 131:192 ) = 0;
            
            
        end
        
        function [gc, prb] = removeNonSig( prb, gc )
           
            lst     = find( prb > .05 );
            gc(lst) = 0;
            prb(lst) = 1;
            
        end
        
        % PAGERANK  Google's PageRank
        % pagerank(U,G,p) uses the URLs and adjacency matrix produced by SURFER,
        % together with a damping factory p, (default is .85), to compute and plot
        % a bar graph of page rank, and print the dominant URLs in page rank order.
        % x = pagerank(U,G,p) returns the page ranks instead of printing.
        % See also SURFER, SPY.
        function x = pagerank(U,G,p)
            
            if nargin < 3, p = .85; end

            % Eliminate any self-referential links

            G = G - diag(diag(G));

            % c = out-degree, r = in-degree

            [n,n] = size(G);
            c = sum(G,1);
            r = sum(G,2);

            % Scale column sums to be 1 (or 0 where there are no out links).

            k = find(c~=0);
            D = sparse(k,k,1./c(k),n,n);

            % Solve (I - p*G*D)*x = e

            e = ones(n,1);
            I = speye(n,n);
            x = (I - p*G*D)\e;

            % Normalize so that sum(x) == 1.

            x = x/sum(x);

            % Bar graph of page rank.

            shg
            bar(x)
            title('Page Rank')

            % Print URLs in page rank order.

            if nargout < 1
               [ignore,q] = sort(-x);
               disp('     page-rank  in  out  url')
               k = 1;
               while (k <= n) & (x(q(k)) >= .005)
                  j = q(k);
                  disp(sprintf(' %3.0f %8.4f %4.0f %4.0f  %s', j,x(j),full(r(j)),full(c(j)),U{j}))
                  k = k+1;
               end
            end
        end
        
        
        function calcAllAllGCvsDistance( gc_values, dist_values )
            points = [pitt.exp.simu.Util.makeOneD(gc_values); pitt.exp.simu.Util.makeOneD(dist_values)]';
            lst = find( points( :,1 ) ~= 0 );
            
            [X1,N1] = hist( points(lst,2), 100 );
            [X2,N2] = hist( pitt.exp.simu.Util.makeOneD(dist_values), N1 );
            
            figure; bar( N1, X1 );
            title('Histogram of GCvalues occruing at distances'); xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
            figure; bar( N2, X2 );
            title('Histogram of distance values' );xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
            figure; bar( N2, X1./X2 );
            title('Histogram of GC values normalized by count of distance values' );xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
        end
        
        %{
        Calculate a 2D plot with one vertex being distance and the other
        being a histogram of GC values at that distance.
        %}
        function calcVerticesAllGCvsDistance( idxs, gc_values, dist_values )
            
            gc_values = [ gc_values( idxs, : ); gc_values( :, idxs )' ];
            dist_values = [ dist_values( idxs, : ); dist_values( :, idxs )' ];
            
            points = [pitt.exp.simu.Util.makeOneD(gc_values); pitt.exp.simu.Util.makeOneD(dist_values)]';
            lst = find( points( :,1 ) ~= 0 );
            
            [X1,N1] = hist( points(lst,2), 100 );
            [X2,N2] = hist( pitt.exp.simu.Util.makeOneD(dist_values), N1 );
            
            figure; bar( N1, X1 );
            title('Histogram of GCvalues occruing at distances'); xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
            figure; bar( N2, X2 );
            title('Histogram of distance values' );xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
            figure; bar( N2, X1./X2 );
            title('Histogram of GC values normalized by count of distance values' );xlabel('Bins of Distance Values'); ylabel('Number of GC values occuring at specified bin' );
            
        end
        
        %{
        Calculate a 2D plot with one vertex being distance and the other
        being a histogram of GC values at that distance. Each point in the
        second dimension is normalized by the total number of GC values
        appearing at that distance.
        %}
        function output = calcAllAllGCvsDistance2D( idxs, gc_values, dist_values )
            
            gc_values   = [ gc_values( idxs, : ); gc_values( :, idxs )' ];
            dist_values = [ dist_values( idxs, : ); dist_values( :, idxs )' ];
            
            lst         = find( gc_values ~= 0 );
            [X,N]       = hist( pitt.exp.simu.Util.makeOneD(gc_values(lst)), 100 );
            
            points      = [pitt.exp.simu.Util.makeOneD(gc_values); pitt.exp.simu.Util.makeOneD(dist_values)]';
            lst         = find( points( :,1 ) ~= 0 );
            points      = points(lst,:);
            
            [X1,N1]     = hist( points(:,2), 100 );
            
            output      = zeros( length(N), length(N1) );
            
            N1          = [0,N1];
            
            for i = 1:length(N1)-1
                
                lst     = find( points( :,2 ) > N1(i) & points( :,2 ) < N1(i+1) );
                length(lst);
                [X2,N2] = hist( points(lst,1 ), 100 );
                size(X2);
                output( :,i ) = (X2);
                
            end
            
            imagesc( output );
            
        end
        
        %{
        Calculate a 2D plot with one vertex being distance and the other
        being a histogram of GC values at that distance. Each point in the
        second dimension is normalized by the total number of GC values
        appearing at that distance.
        %}
        function output = calcAllAllGCvsDistance2DNormalized( idxs, gc_values, dist_values )
            
            gc_values   = [ gc_values( idxs, : ); gc_values( :, idxs )' ];
            dist_values = [ dist_values( idxs, : ); dist_values( :, idxs )' ];
            
            lst         = find( gc_values ~= 0 );
            [X,N]       = hist( pitt.exp.simu.Util.makeOneD(gc_values(lst)), 100 );
            
            points      = [pitt.exp.simu.Util.makeOneD(gc_values); pitt.exp.simu.Util.makeOneD(dist_values)]';
            lst         = find( points( :,1 ) ~= 0 );
            points      = points(lst,:);
            
            [X1,N1]     = hist( points(:,2), 100 );
            
            output      = zeros( length(N), length(N1) );
            
            N1          = [0,N1];
            
            for i = 1:length(N1)-1
                
                lst     = find( points( :,2 ) > N1(i) & points( :,2 ) < N1(i+1) );
                [X2,N2] = hist( points(lst,1 ), 100 );
                size(X2);
                output( :,i ) = (X2./X1(i));
                
            end
            
            imagesc( output );
            
        end
        
        function output = calcAllAllGCvsDistanceCrossIdxs2D( idxs1, idxs2,  gc_values, dist_values, hist_bins )
            
            idxs        = pitt.exp.simu.Util.getIdxValuesAtIntersection( gc_values, idxs1, idxs2 );
            gc_values   = gc_values( idxs );
            dist_values = dist_values( idxs );
            
            lst         = find( gc_values ~= 0 );
            [X,N]       = hist( pitt.exp.simu.Util.makeOneD(gc_values(lst)), hist_bins );
            
            points      = [pitt.exp.simu.Util.makeOneD(gc_values), pitt.exp.simu.Util.makeOneD(dist_values)];
            lst         = find( points( :,1 ) ~= 0 );
            points      = points(lst,:);
            
            [X1,N1]     = hist( points(:,2), hist_bins );
            
            output      = zeros( length(N), length(N1) );
            
            N1          = [0,N1];
            
            for i = 1:length(N1)-1
                
                lst     = find( points( :,2 ) > N1(i) & points( :,2 ) < N1(i+1) );
                [X2,N2] = hist( points(lst,1 ), 100 );
                size(X2);
                output( :,i ) = (X2);
                
            end
            
            imagesc( output );
            
        end
        
        function awesomeDisplay( twoDhist )
            
            
            subplot( 4,3,[ 2 3 5 6] );
            imagesc( twoDhist );
            title( 'Image of GC values distribution as function of distances' );
            
            subplot( 4,3,[1 4] );
            plot( flipud( sum( twoDhist, 2 ) ), [1:length(twoDhist)] );
            title( 'Marginal sum of GC values across all distances' );
            
            subplot( 4,3, [8, 9] );
            plot( sum( twoDhist,1 ) );
            title( 'Marginal sum of GC values at distances' );
            
            subplot( 4,3, [11 12] );
            counts = zeros( 1, size(twoDhist, 2 ) );
            for i = 1: size(twoDhist, 2)
                counts(i) = length( find( twoDhist(:,i) > 0 ) );
            end
            plot( counts );
            title( 'Count of non-zero GC values -- WRONG' ); % Actually count of non-zero histogram values
            
        end
        
        %{
        
        calculate the pwelch at each vertex across all timepoints in data
        
        %}
        function output = chan_Pwelch( data )
           
            [row, col ] = size( data );
            output = [];
            for i = 1:row;
                d = pwelch( data(i,:) );
                output = [output; d'];
            end
        end
        
        
        function displayGCPRBGeneric( gc, prb )
           
            figure;
            subplot( 2,2,1  );
            imagesc( gc );
            title('Raw map of GC values' );
            
            subplot( 2,2,2 );
            imagesc( prb );
            title('Raw map of PRB values' );
            
            subplot( 2,2,3 );
            pitt.exp.simu.Util.hist2D( gc );
            title( 'Histogram of GC values' );
            
            
            subplot( 2,2,4 );
            pitt.exp.simu.Util.hist2D( prb );
            title( 'Histogram of prb values' );
            
        end
        
    end
    
end