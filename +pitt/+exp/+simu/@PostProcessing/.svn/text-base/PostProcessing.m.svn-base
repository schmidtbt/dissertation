classdef PostProcessing < handle
   
    methods( Static )
        
        % Provided a distributed master job which used the
        % pitt.exp.simu.Processing.gregresskernel, reintegrate to obtaina
        % single output file which is a struct containing a prb and gc
        % component.
        function [prb, gc, completed] = reintegrateSimuKernel( master, DISPLAY )
            
            
            % Default Don't Display
            if( nargin == 1 )
                DISPLAY = 0;
            end
            
            
            if( DISPLAY )
                fprintf('Progress will be monitored\n' );
                figure( 143 ); 
            else
                fprintf( 'Progress will not be monitored\n' );
            end
            
            % Allocate two storage variables
            numvertices     = size( master.data, 1 );
            prb             = -1.*ones( numvertices,numvertices );
            gc              = -1.*ones( numvertices,numvertices );
            completed       = zeros( numvertices,numvertices );
            cycle           = 1;
            files           = dir([master.jobPath, 'output*']);
            fprintf( 'Found %i output files\n', length(files ) );
            tic;
            for i = 1:length(files)
                fprintf( 'Loading: %s\n', files(i).name );
                loaded = load( files(i).name );
                fprintf( '\t%i of %i\n', i, length(files) );
                
                [rows, cols] = size( loaded.output );
                fprintf( '\tBlock size: %i x %i\n', rows, cols );
                
                % Loop over block output data
                for r = 1:rows
                    for c = 1:cols
                        
                        ridx = loaded.rowidxs(r);
                        cidx = loaded.colidxs(c);
                        
                        TwoOneGC = loaded.output{r,c}.gc(2,1);
                        OneTwoGC = loaded.output{r,c}.gc(1,2);
                        
                        TwoOnePrb = loaded.output{r,c}.prb(2,1);
                        OneTwoPrb = loaded.output{r,c}.prb(1,2);
                        
                        prb( ridx, cidx ) = TwoOnePrb;
                        prb( cidx, ridx ) = OneTwoPrb;
                        
                        gc( ridx, cidx ) = TwoOneGC;
                        gc( cidx, ridx ) = OneTwoGC;
                        
                        % Update the completion matrix to check and ensure
                        % progress
                        completed( ridx, cidx ) = cycle;
                        completed( cidx, ridx ) = cycle;
                        
                    end
                    
                end
                
                
                if( DISPLAY )
                    if( mod( i, 100 ) == 0 || i == 1 )
                        % Update teh cycle value, to display differences
                        % between updates
                        cycle = cycle+1;
                        fprintf( 'Displaying Progress\n' );
                        figure( 143 ); imagesc( completed );
                        elapsedtime = toc;
                        complete = i/length(files);
                        estcompletion = (elapsedtime/complete); % in seconds
                        estremaining = estcompletion - elapsedtime; %in seconds
                        title( sprintf( 'Progress: %i of %i. Elapsed time: %2.3f (seconds), Est Complete in: %2.3f (seconds)', i, length(files), elapsedtime, estremaining ) );
                        drawnow;
                    end
                end
                
            end
            
        end
        
        
    end
    
end