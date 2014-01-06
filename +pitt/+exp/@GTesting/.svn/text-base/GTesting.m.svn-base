classdef GTesting < handle
    
    properties
        
    end
    
    methods (Static)
        
        function output = testSingleCasual
            
            R = randn(1, 10000 );
            X = filter( 1, [1, .5, -.2], R );
            Y = .2*randn( 1,10000 ) + [randn(1), X(2:end)] + filter( 1, [1, .3, -.2], R );
            
            output = cca_granger_regress( [X; Y], 5 )
            
        end
        
        function dispSigHistograms( dsize )
            
            if( nargin ==  0 )
                dsize = 100;
            end
            
            OneTwo = 0;
            TwoOne = 0;
            OneTwoGC = [];
            TwoOneGC = [];
            
            for i = 1:dsize
                sig = pitt.exp.GTesting.simulateTwoSources;
                ret = cca_granger_regress( sig, 5 );
                
                
                if( ret.prb(1,2) <.05 )
                    TwoOne = TwoOne + 1;
                    TwoOneGC = [TwoOneGC; ret.gc(1,2) ];
                end
                
                if( ret.prb(2,1) <.05 )
                    OneTwo = OneTwo + 1;
                    OneTwoGC = [OneTwoGC; ret.gc(2,1) ];
                end
                
                
            end
            
            
            figure(1); hist( OneTwoGC, 100 );
            title( 'One->Two GC Distribution - Causal' );
            
            figure(2); hist( TwoOneGC, 100 );
            title( 'Two->One GC Distribution - Causal' );
            
            
        end
        
        function [out] = simulateTwoSources
            
            R = randn(1, 10000 );
            X = filter( 1, [1, .5, -.2], R );
            Y = .2*randn( 1,10000 ) + [randn(1,3), X(1:end-3)] + filter( 1, [1, .3, -.2], R );
            %Y = .2*randn( 1,10000 ) + [randn(1,3), X(1:end-3)];
            %Y = [randn(1,3), X(1:end-3)];
            out = [X;Y];
            % Here the signal X should GC the signal Y
        end
        
        % Rsize = Number of ranks to use
        function output = testRandData( Rsize, nlags )
            
            if( nargin == 0 )
                Rsize = 5;
                nlags = 5;
            end
            
            pitt.Depend.CCAadd
            
            Rdata = pitt.exp.GTesting.genRandTestData( Rsize );
            
            output = cell( 0 );
            
            for i = 1:Rsize
                for j = 1:Rsize
                    if( i == j ); continue; end;
                    fprintf( '%i - %i | \n', i, j );
                    output{i,j} = cca_granger_regress( [Rdata( i, : ); Rdata(j,:)], nlags );
                end
            end
            
        end
        
        
        function data = genRandTestData( Rsize )
            
            if( nargin == 0 )
                Rsize = 5;
            end
            
            data = randn( Rsize, 10000 );
            
        end
        
        % Rsize is number of pairs of casual signals
        function output = testCasual( Rsize, nlags )
            
            if( nargin == 0 )
                Rsize = 5;
                nlags = 5;
            end
            
            pitt.Depend.CCAadd
            
            Rdata = [];
            for i = 1:Rsize
                fprintf( 'Simulating %i\n', i );
                Rdata = [Rdata; pitt.exp.GTesting.simulateTwoSources];
            end
            
            output = cell( 0 );
            
            for i = 1:Rsize*2
                for j = 1:Rsize*2
                    if( i == j ); continue; end;
                    fprintf( '%i - %i | \n', i, j );
                    output{i,j} = cca_granger_regress( [Rdata( i, : ); Rdata(j,:)], nlags );
                end
            end
            
        end
        
        function output = calcOutputMatrix( Rdata, nlags )
            Rsize = size(Rdata,1);
            for i = 1:Rsize
                for j = 1:Rsize
                    if( i == j ); continue; end;
                    fprintf( '%i - %i | \n', i, j );
                    output{i,j} = cca_granger_regress( [Rdata( i, : ); Rdata(j,:)], nlags );
                end
            end
        end
        
        function kernel = gregressKernel
            kernel = legion.Kernel;
            kernel.add( @pitt.signal.Util.varnorm, 'X' );
            kernel.add( @cca_granger_regress, 'X', 10 );
            kernel.addEnviro( @pitt.Depend.CCAadd );
        end
        
        %Given an output system of size x size of an all-all calculation
        function prb = extractPrb( output )
            
            [x,y] = size(output);
            
            prb = zeros( x,y );
            
            for i = 1:x
                for j = 1:y
                    if( i >= j ); continue; end;
                    prb( i,j ) = output{j,i}.prb( 1,2 );
                    prb( j,i ) = output{j,i}.prb( 2,1 );
                end
            end
        end
        
        %Given an output system of size x size of an all-all calculation
        function gc = extractGc( output )
            
            [x,y] = size(output);
            
            gc = zeros( x,y );
            
            for i = 1:x
                for j = 1:y
                    if( i >= j ); continue; end;
                    gc( i,j ) = output{j,i}.gc( 1,2 );
                    gc( j,i ) = output{j,i}.gc( 2,1 );
                end
            end
        end
        
        function histGC( output )
            gc = pitt.exp.GTesting.extractGc( output );
            hist( reshape( gc, 1, size(gc,1)*size(gc,2) ), 100 );
        end
        
        function histPrb( output )
            prb = pitt.exp.GTesting.extractPrb( output );
            hist( reshape( prb, 1, size(prb,1)*size(prb,2) ), 50 );
        end
        
        function histGCSigOnly( output )
            
            [x,y] = size(output);
            
            dataOT = [];
            dataTO = [];
            
            for i = 1:x
                for j = 1:y
                    if( i>=j ); continue; end;
                    
                    if( output{j,i}.prb(1,2) < .05 )
                        dataTO = [dataTO; output{j,i}.gc(1,2)];
                    end
                    
                    if( output{j,i}.prb(2,1) < .05 )
                        dataOT = [dataOT; output{j,i}.gc(2,1)];
                    end
                end
                
            end
            
            figure;
            hist( dataOT, 100 );
            title('Two->One GC values at Sig level .05');
            
            figure;
            hist( dataTO, 100 );
            title('One->Two GC values at Sig level .05');
            
            lst = find( dataOT > .1 );
            figure;
            hist( dataOT(lst), 100 );
            title('Two->One GC values at Sig level .05 -- High GC');
            
            lst = find( dataTO > .1 );
            figure;
            hist( dataTO(lst), 100 );
            title('One->Two GC values at Sig level .05 -- High GC');
            
            figure; 
            hist( [dataTO; dataOT], 100 );
            title('Significant GC values between label regions' )
            
        end
        
        function histGCSigOnlyPrbGCInput( gc, prb )
            
            lst = find( prb <= .05 );
            figure;
            hist( gc(lst), 100 );
            
            figure;
            hist( dataOT, 100 );
            title('Two->One GC values at Sig level .05');
            
            figure;
            hist( dataTO, 100 );
            title('One->Two GC values at Sig level .05');
            
            lst = find( dataOT > .1 );
            figure;
            hist( dataOT(lst), 100 );
            title('Two->One GC values at Sig level .05 -- High GC');
            
            lst = find( dataTO > .1 );
            figure;
            hist( dataTO(lst), 100 );
            title('One->Two GC values at Sig level .05 -- High GC');
            
            figure; 
            hist( [dataTO; dataOT], 100 );
            title('Significant GC values between label regions' )
            
        end
        
        % given the return of a single cca_granger_regress output
        function gcreport( ret )
            
            fprintf( 'Source 1 -> Source 2\n\n' );
            if( ret.prb(1,2) < .05 )
                fprintf( '\t  ||-- SIG --||\n' );
            else
                fprintf( '\n' );
            end
            fprintf( '\t prb: %2.7f\n' , ret.prb(1,2));
            fprintf( '\t gc : %2.7f\n\n', ret.gc(1,2) );
            
            fprintf( '\nSource 2 -> Source 1\n\n' );
            if( ret.prb(2,1) < .05 )
                fprintf( '\t  ||-- SIG --||\n' );
            else
                fprintf( '\n' );
            end
            fprintf( '\t prb: %2.7f\n' , ret.prb(2,1));
            fprintf( '\t gc : %2.7f\n\n', ret.gc(2,1) );
            
        end
        
        function [leg, data] = testDistributedJobRunner( seed )
            
            if( nargin == 0 )
                seed = randi( 1000, 1 );
            end
            fprintf('Executing seed: %i', seed );
            % use a cca_granger kernel
            kernel = pitt.exp.GTesting.gregressKernel;
            data = pitt.exp.GTesting.genRandTestData( 5 );
            numcores = 2;
            leg = legion.DistrMaster( kernel, data, numcores, sprintf('test_distr_runner_integrate_%i',seed) );
            leg.run();
            
        end
        
        function result = orderDependence( X )
            result = X(1) - X(2);
        end
        
    end
    
end