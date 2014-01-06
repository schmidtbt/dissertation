classdef GC < handle
    
    
    
    methods (Static)
        
        %{
        @param X numChans x numTimepoints
        @return image 2D array of significant GCs with their magnitude
        displayed
        %}
        function image = allBivarGC( X, MODELORDER )
            numChans = size(X,1);
            
            % Square by number of channels
            image = zeros( size(X,1) );
            
            progressbar(0,0,0);
            
            %Loop over channels
            for i = 1:size(X, 1 )
                progressbar( i/numChans, [],[] );
                for j = 1:size(X,1)
                    % continue if crossing diagonal or equal
                    if( i >= j )
                        continue;
                    end
                    progressbar( [],j/numChans,[] );
                    ret = cca_granger_regress( X([i,j], :), MODELORDER,1 );
                    
                    if( ret.prb( 2,1 ) > .1 )
                        image( j,i ) = ret.gc( 2,1 );
                    end
                    
                    if( ret.prb( 1,2 ) > .1 )
                        image( i,j ) = ret.gc( 1,2 );
                    end
                    
                    progressbar( [],[], i / size(X,1) );
                    
                end
                
            end
            
        end
        
        function jobArray = runJobBivarGC( X, MODELORDER )
            
            numChans = size(X,1);
            jobArray = [];
            % Square by number of channels
            image = zeros( size(X,1) );
            progressbar( 0,0 );
            for i = 1:size(X, 1 )
                progressbar( i/numChans,[]);
                for j = 1:size(X,1)
                    progressbar( [], j/numChans );
                    
                    % continue if crossing diagonal or equal
                    if( i >= j )
                        continue;
                    end
                    
                    job = pitt.jobs.GCJob( [num2str(i),' ',num2str(j)], i,j, X(i,:), X(j,:), MODELORDER );
                    jobArray = [jobArray, job];
                    
                end
                
            end
            
            
            jr = pitt.JobRunner;
            jr.parallel( true );
            jobArray = jr.runJob( jobArray );
            
        end
        
    end
    
end