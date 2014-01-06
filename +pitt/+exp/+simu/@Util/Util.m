classdef Util < handle
    
    methods( Static )
        
        %{
        Given idxs1 as vertices in orig, find all intersecting indexes in idxs2
        and REPLACE with values
        
        Values can be either single scalar, or size of intersection
        (length of idxs1 * idxs2)
        %}
        function orig = replace2Dcoords( orig, idxs1, idxs2, values )
            
            idx = pitt.exp.SimTesting.getIdxValuesAtIntersection( orig, idxs1, idxs2 );
            orig(idx) = values;
        end
        
        %{
        Given idxs1 as vertices in orig, find all intersecting indexes in idxs2
        and return their VALUES
        
        Values can be either single scalar, or size of intersection
        (length of idxs1 * idxs2)
        %}
        function values = getValuesAtIntersection( orig, idxs1, idxs2 )
            
            idx = pitt.exp.simu.Util.getIdxValuesAtIntersection( orig, idxs1, idxs2 );
            values = orig(idx);
        end
        
        %{
        Given idxs1 as vertices in orig, find all intersecting indexes in idxs2
        and return their INDEXES
        
        Values can be either single scalar, or size of intersection
        (length of idxs1 * idxs2)
        %}
        function idx = getIdxValuesAtIntersection( orig, idxs1, idxs2 )
            coords = [];
            for i = 1:length(idxs1)
               for j = 1:length(idxs2)
                  coords = [coords; [idxs1(i),  idxs2(j)] ];
               end
            end
            
            idx = sub2ind( size( orig ), coords(:,1), coords(:,2) );
        end
        
        %{
        Plot a histogram with n bins of each array value in array
        Accepts both 1D and 2D array input
        Noramlizes the X axis to be [0,1] interval
        %}
        function normXHist2D( array, n )
            array = pitt.exp.simu.Util.makeOneD( array );
            if( nargin == 1 )
                n = 100;
            end
           hist( reshape( array, 1, size(array,1)*size(array,2) ), n );
           
           pitt.exp.simu.Util.normalizeXAxis();
           
        end
        
        %{
        Plot a histogram with n bins of each array value in array
        Accepts both 1D and 2D array input
        %}
        function hist2D( array, n )
            array = pitt.exp.simu.Util.makeOneD( array );
            if( nargin == 1 )
                n = 100;
            end
           hist( reshape( array, 1, size(array,1)*size(array,2) ), n );
           
        end
        
        function normalizeXAxis()
           ax = axis();
           axis([0 1 ax(3) ax(4)]); 
        end
        
        % Given two 1-D arrays, find the values of needle occuring in
        % haystack
        % Returns the co-occuring values
        function occur = findOverlapVertices( needle, haystack )
            
            occur = [];
            
            for i = 1:length(needle)
               val = find( haystack == needle(i) );
               if( length( val ) > 0 )
                   occur = [occur; needle(i)];
               end
            end
           
        end
        
        % Given two 1-D arrays, find the values of needle occuring in
        % haystack
        % Returns the co-occuring values
        function occur = findOverlapVerticesIndexes( needle, haystack )
            
            occur = [];
            
            for i = 1:length(needle)
               val = find( haystack == needle(i) );
               if( length( val ) > 0 )
                   occur = [occur; val];
               end
            end
           
        end
        
        
        
        
        % Load the sim data object, assumed in ~/simu directory
        function sim = loadSimData( sim_directory )
            sim = pitt.Simu( 'subj', sprintf('~/simu/%s', sim_directory), 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' ); 
            sim.post_simulate();
        end
        
        % Get label 1 vertices
        function vertices = label1Vertices( sim )
            vertices = sim.dipole_loc_label1;
        end
        
        % Get label 2 vertices
        function vertices = label2Vertices( sim )
            vertices = sim.dipole_loc_label2;
        end
        
        % Given two arrays, plot as 'r' and 'g' histograms on same axis
        % Accepts 2D or 1D inputs
        function overlayHistograms( data1, data2 )
            
            data1 = pitt.exp.simu.Util.makeOneD( data1 );
            data2 = pitt.exp.simu.Util.makeOneD( data2 );
            
            [N1,X1] = hist( data1, 100 );
            [N2,X2] = hist( data2, 100 );
            
            bar( X1, N1, 'r' );
            hold on;
            bar( X2, N2, 'g' );
            hold off;
            
            
        end
        
        function overlayHistogramsRemoveZeros( data1, data2 )
            
            data1 = pitt.exp.simu.Util.makeOneD( data1 );
            data2 = pitt.exp.simu.Util.makeOneD( data2 );
            
            lst1  = find(data1 > 0);
            lst2  = find(data2 > 0);
            
            
            [N1,X1] = hist( data1(lst1), 100 );
            [N2,X2] = hist( data2(lst2), 100 );
            
            bar( X1, N1, 'r' );
            hold on;
            bar( X2, N2, 'g' );
            hold off;
            
            
        end
        
        function overlayHistogramsRemoveZerosNormY( data1, data2 )
            
            data1 = pitt.exp.simu.Util.makeOneD( data1 );
            data2 = pitt.exp.simu.Util.makeOneD( data2 );
            
            lst1  = find(data1 > 0);
            lst2  = find(data2 > 0);
            
            
            [N1,X1] = hist( data1(lst1), 100 );
            [N2,X2] = hist( data2(lst2), 100 );
            
            N1 = N1 ./sum(sum(N1));
            N2 = N2 ./sum(sum(N2));
            
            bar( X1, N1, 'r' );
            hold on;
            bar( X2, N2, 'g' );
            hold off;
            
            
        end
        
        % Check if OneD or TwoD array is provided, and return a 1D array [1
        % x N]
        function oneD = makeOneD( possTwoDArray )
            
            if( isvector( possTwoDArray ) )
                oneD = possTwoDArray;
            else
                oneD = reshape( possTwoDArray, 1, size( possTwoDArray, 1 ) * size( possTwoDArray, 2 ) );
            end
            
        end
        
        
        %{
        
        given an input gc, remove all values which occur across the
        diagonal (they cause each other)
        
        %}
        function gct = removeCrossDuplicates( gc )
            
            gct = gc;
            [row, col] = size( gct );
            
            for i = 1:row
                for j = 1:col
                    if( i > j ); continue; end;
                    
                    if( gct(i,j) > 0 && gc(j,i) > 0 )
                        gct( i,j ) = 0;
                        gct( j,i ) = 0;
                    end
                    
                end
            end
            
        end
        
        %{
        
        Given two equally sizes matrices, return a 2x(sizeof input) array
        of points of matrix1 value by matrix 2 value
        
        %}
        function points = scatterPointsAsIndices( matrix1, matrix2 )
           
            [row, col] = size( matrix1 );
            points = [];
            counter = 1;
            for i = 1:row
                for j = 1:col
                    
                    if( matrix1(i,j) == 0 ); continue; end; 
                    
                    points( counter, : ) = [ matrix1(i,j), matrix2(i,j) ];
                    counter = counter + 1;
                end
            end
            
            
        end
        
    end
    
end