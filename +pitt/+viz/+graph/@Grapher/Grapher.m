classdef Grapher < handle
    %Example:
    %pitt.viz.graph.Grapher.genDotFileWithPosLabels( 'filename.dot', randn(3,3), [0], floor([rand(3,2)*5]), Labels );
    
    methods (Static)
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        %}
        function genDotFile( Outfilename, X, threshold)
            
            if( nargin == 2 )
                threshold = -10e10;
            end
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'digraph g{\n' );
            
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                     if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [num2str(i), '->',num2str(j),';\n'] );
                    end
                end
            end
            
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        function genUndirFileFromList( Outfilename, X, Labels, Colors )
            
            
            if( nargin == 2 )
                Labels = [];
                Colors = [];
            end
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'graph g{\n' );
            
            uniq = unique( [X(:,1),X(:,2)] );
            
            for i = 1:length(uniq)
                fprintf( fid, '%s %s;\n', pitt.viz.graph.Grapher.convertToLabel( uniq(i), Labels ),  pitt.viz.graph.Grapher.convertColor( uniq(i), Colors ) );
            end
            
            for i = 1:size(X,1)
                fprintf( fid, [ pitt.viz.graph.Grapher.convertToLabel( X(i,1), Labels ), '--',pitt.viz.graph.Grapher.convertToLabel( X(i,2), Labels ),';\n'] );
            end
            
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        function genUndirFileFromListLineThickness( Outfilename, X, Labels, Colors )
            
            
            if( nargin == 2 )
                Labels = [];
                Colors = [];
            end
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'graph g{\n' );
            
            uniq = unique( [X(:,1),X(:,2)] );
            
            for i = 1:length(uniq)
                fprintf( fid, '%s %s;\n', pitt.viz.graph.Grapher.convertToLabel( uniq(i), Labels ),  pitt.viz.graph.Grapher.convertColor( uniq(i), Colors ) );
            end
            
            for i = 1:size(X,1)
                
                fprintf( fid, [ pitt.viz.graph.Grapher.convertToLabel( X(i,1), Labels ), '--',pitt.viz.graph.Grapher.convertToLabel( X(i,2), Labels ),';\n'] );
                
            end
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        
        
        
        function genStuff( Outfilename, X )
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'graph g{\n' );
            
            fprintf( fid, 'splines=true;\nsep="+25,25";\noverlap=scalexy;\nnodesep=0.9;\ngraph [ dpi = 300 ];\n')
            
            lst = find( sum(X) > 0 );
            
            for i = 1:length(lst)
                
                if( mod( lst(i), 2 ) == 0 )
                    HEMI = 'LH';
                else
                    HEMI = 'RH';
                end
                
                if strcmp( HEMI, 'LH' )
                    color = 'blue';
                else
                    color = 'red';
                end
                
                vi = pitt.viz.graph.Grapher.convertVal(lst(i));
                
                fprintf( fid, '"%s" [fillstyle=fill, color=%s];\n', vi, color );
                
            end
            
            X = tril(X);
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                    
                    if( X(i,j) > 0 )
                        
                        vi = pitt.viz.graph.Grapher.convertVal(i);
                        vj = pitt.viz.graph.Grapher.convertVal(j);
                        
                        fprintf( fid, '"%s" -- "%s" [penwidth=%i];\n', vi, vj, X(i,j).*2 );
                    end
                    
                end 
            end
            
            fprintf( fid, '}\n' );
            fclose( fid );
            
        end
        
        function v = convertVal( val )
            val = num2str(val);
            v = val;
            v = v(1:end-1);
            h = val(end);
            if( h == '0' )
                hem = 'LH';
            else
                hem = 'RH';
            end
            v = sprintf( 'BA %s | %s', v, hem );
        end
        
        
        
        
        
        
        
        
        
        
        function label = convertToLabel( idx, Labels )
            
            label = idx;
            
            if( ~isempty(Labels) )
                label = sprintf( '"%s"', Labels{idx});
            end
            
        end
        
        function color = convertColor( idx, Colors )
            
            color = '';
            
            if( ~isempty(Colors) )
                
                color = sprintf( '[ style=filled, color=%s ]', Colors{idx} );
            end
            
        end
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        %}
        function genNeatoFile( Outfilename, X, threshold)
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'graph g{\n' );
            
            
            for i = 1:size(X,1)
                if( mod(i,100) == 0); fprintf('%i\n', i); end;
                for j = 1:size(X,2)
                    
                     if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [num2str(i), '--',num2str(j)] );
                        %fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n');
                    end
                end
            end
            
            fprintf('\n');
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        %}
        function genNeatoTwoHemiFreqFile( Outfilename, X, threshold)
            
            
            if( nargin == 2 )
                threshold = -10e10;
            end
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'graph g{\n' );
            
            lhidx = 1:floor(size(X,1)/2);
            rhidx = floor(size(X,1)/2)+1:size(X,1);
            
            
            fprintf( fid, 'subgraph clusterLH {\n' );
            fprintf( fid, 'node [style=filled];\n');
            for i=lhidx
                fprintf( fid, ['\t',num2str(i)] );
                fprintf( fid, ';\n');
            end
            
            for i = lhidx
                for j = lhidx
                    
                    if( j >= i ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, ['\t',num2str(i), '--',num2str(j)] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n');
                    end
                    
                end
            end
            fprintf( fid, '}\n');
            
            fprintf( fid, 'subgraph clusterRH {\n' );
            
            
            for i=rhidx
                fprintf( fid, ['\t',num2str(i)] );
                fprintf( fid, ';\n');
            end
            
            for i = rhidx
                for j = rhidx
                    
                    if( j >= i ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, ['\t',num2str(i), '--',num2str(j)] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n');
                    end
                    
                end
            end
            fprintf( fid, '}\n');
            
            for i = rhidx
                for j = lhidx
                    
                    %if( j > i ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [num2str(i), '--',num2str(j)] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n');
                    end
                    
                end
            end
            
            fprintf('\n');
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        function lst = freqFileToList( X, threshold)
            
            
            if( nargin == 1 )
                threshold = -10e10;
            end
            
            lst = [];
            
            for i = 1:size(X,1)
                
                for j = 1:size(X,2)
                    
                    if( j >= i ) continue; end;
                    
                    if( X(i,j) > threshold )
                        lst = [lst; [i,j]];
                    end
                    
                end
                
            end
            
            
        end
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        Pos is a numChannels x 2 matrix of [X,Y] pairs for each channel in
        X
        %}
        function genDotFileWithPos(Outfilename, X, threshold, Pos)
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'digraph g{\n' );
            
            for i = 1:size(Pos,1)
                fprintf( fid, [num2str(i), '[ pos= "',num2str(Pos(i,1)),',',num2str(Pos(i,2)),'!" ]\n']);
            end
            
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                     if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [num2str(i), '->',num2str(j),';\n'] );
                    end
                end
            end
            
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        Pos is a numChannels x 2 matrix of [X,Y] pairs for each channel in
        X
        Labels is a cellarray of strings with label names
        %}
        function genDotFileWithPosLabels(Outfilename, X, threshold, Pos, Labels)
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'digraph g{\n' );
            
            for i = 1:size(Pos,1)
                fprintf( fid, [Labels{i}, '[ pos= "',num2str(Pos(i,1)),',',num2str(Pos(i,2)),'!" ]\n']);
            end
            
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                    
                    if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [Labels{i}, '->',Labels{j}] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n' );
                    end
                end
            end
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        Pos is a numChannels x 2 matrix of [X,Y] pairs for each channel in
        X
        Labels is a cellarray of strings with label names
        %}
        function genDotFileWithLabels(Outfilename, X, threshold, Labels)
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'digraph g{\n' );
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                    
                    if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [Labels{i}, '->',Labels{j}] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n' );
                    end
                end
            end
            
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        
        %{
        Outfilename is a string for the file to output to (usually with
        .dot extension)
        X is a numChannels x numChannels matrix with values corresponding
        to weights for that pair
        threshold is the minimum value at a i,j point in X to draw a line
        Pos is a numChannels x 3 matrix of [X,Y,Z] values for each channel in
        X
        Labels is a cellarray of strings with label names
        %}
        function gen3DDotFileWithPosLabels(Outfilename, X, threshold, Pos,Labels)
            
            fid = fopen( Outfilename, 'w' );
            fprintf( fid, 'digraph g{\n' );
            
            for i = 1:size(Pos,1)
                fprintf( fid, [Labels{i}, '[ pos= "',num2str(Pos(i,1)), ',' ,num2str(Pos(i,2)), ',', num2str(Pos(i,3)),'!" ]\n']);
            end
            
            for i = 1:size(X,1)
                for j = 1:size(X,2)
                     if( i == j ); continue; end;
                    if( X(i,j) > threshold )
                        fprintf( fid, [Labels{i}, '->',Labels{j}] );
                        fprintf( fid, ['[ color=', pitt.viz.graph.Grapher.colorSwatch(X(i,j)), ']']);
                        fprintf( fid, ';\n' );
                    end
                end
            end
            
            
            fprintf( fid, '}\n' );
            fclose( fid );
        end
        
        
        function colorString = colorSwatch( value )
            out = 'black';
            if( value < .1 )
                out = 'blue';
            elseif( value <.2 )
                out = 'green';
            elseif( value <.8 )
                out = 'orange';
            else
                out = 'red';
            end
            colorString = out;
        end
        
    end
    
    
end