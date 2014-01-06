classdef FParam < handle
    
    
    methods (Static)
        
        function [alpha, rsquare] = calculate( input )
            
            if( size( input,1 ) > size( input,2 ) )
                input = input';
            end
            %=figure; pwelch( input );
            data = pwelch( input );
            data = data';
            
            X               = mag2db( 1:length(data) );
            magData  = mag2db( data );
            fit = polyfit( X, magData, 1 );
            pred = polyval( fit, X );
            cc = corrcoef( magData, pred );
            
            %{
            figure;
            scatter(X, magData ); hold on;
            plot( X, pred,'k' ); hold off;
            %}
            alpha = -fit(1);
            rsquare = cc(1,2).^2;
            
        end
        
        
        function [alphas, rsquares] = validateWhiteNoise( NTrials )
            
            alphas = [];
            rsquares = [];
            
            for i = 1:NTrials
                r = randn(1,10000);
                [alpha, rsquare] = pitt.exp.FParam.calculate( r );
                alphas = [alphas; alpha];
                rsquares = [rsquares; rsquare];
            end
            
        end
        
        
        function [alphas, rsquares] = validateRandomWalk( NTrials )
            
            alphas = [];
            rsquares = [];
            
            for i = 1:NTrials
                r = randn(1,10000);
                walk = filter( 1, [1 -1], r );
                [alpha, rsquare] = pitt.exp.FParam.calculate( walk );
                alphas = [alphas; alpha];
                rsquares = [rsquares; rsquare];
            end
            
        end
        
        function output = runFParam( data, numScales )
            s = pitt.jobs.Scale( data, pitt.exp.FParam.genFParamPipe() );
            output = s.runByWindow( numScales, 0 );
        end
        
        % Detrend each part
        function pipe = genFParamPipe()
            
            pipe = pitt.jobs.Pipe; 
            pipe.add( @cca_detrend, 'X');
            pipe.add( @pitt.exp.FParam.calculate, 'X');
            
        end
        
        %Execute on local host to begin distributed job
        function distrRun( subject, jobName, numscales )
            grid = pitt.grid.Grid( jobName );
            grid.setExecutable( [ 'pitt.exp.FParam.remoteRun( ''',subject,''', ', num2str(numscales),' );' ] );
            grid.executePBS();
        end
        
        % Executed on remote host to begin processing there
        function remoteRun( subject, numscales )
            
            disp('loading');
            dobj = pitt.data.meg.RestData( subject );
            data = dobj.getData();
            disp('executing');
            output = pitt.exp.FParam.runFParam( data, numscales );
            disp('saving');
            save( 'output.mat', 'output', '-v7.3'  );
            
        end
        
        % Save data then re-load it later. Executed locally
        function distrRun2( data, jobName, numscales )
            grid = pitt.grid.Grid( jobName );
            save( ['/synapse/logs/schmidtb/',jobName,'/data.mat'], 'data' );
            grid.setExecutable( [ 'pitt.exp.FParam.remoteRun2( ', num2str(numscales),' );' ] );
            grid.executePBS();
        end
        
        % Executed on remote host to begin processing there (via load
        % command of previously saved data
        function remoteRun2( numscales )
            
            disp('loading');
            load('data.mat');
            disp('executing');
            output = pitt.exp.FParam.runFParam( data, numscales );
            disp('saving');
            save( 'output.mat', 'output', '-v7.3'  );
            
        end
        
    end
    
end