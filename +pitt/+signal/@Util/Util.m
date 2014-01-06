classdef Util < handle
    
    
    methods (Static)
        
        function output = varnorm( input )
            output = input ./ (ones( length(input),1 ) * var(input'))';
        end
        
        function output = addpinc( X )
            output = X + pincnoise(length(X))';
        end
        
        function output = addPinkNoise( X )
            Nx = 2^16;  % number of samples to synthesize
            B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
            A = [1 -2.494956002   2.017265875  -0.522189400];
            x = filter(B,A,X);    % Apply 1/F roll-off to PSD
            output = x;
        end
        
        function output = notchFilter( Data, Fs, RemoveFreq )
            output = Data';
            for i = 1: length(RemoveFreq )
                Wo = RemoveFreq(i)/(Fs/2);  BW = Wo/40;
                [b,a] = iirnotch(Wo,BW);
                output = filter( b,a, output );
                %figure(501); clf; pwelch( output );
            end
            output = output';
            
        end
        
        function output = addWgnNoise( X, P )
            output = X + wgn( size(X,1), size(X,2), P );
        end
        
        % Good order ~ 8, recursions ~ 3
        function model = getPinkNoiseModel( recursions, order )
            p = pitt.Pipe;
            for i = 1:recursions
                p.add( @pitt.signal.Util.addPinkNoise, 'X' );
            end
            p.initial( randn(1,10000) );
            fnoise = p.execute;
            model = ar( fnoise, order );
        end
        
        
        function output = decimateSeries( series, order )
            output = [];
            % Row loop
            for i = 1:size(series, 1 )
                output= [output; decimate( series(i, :), order )];
            end
            
        end
        
        
        function avePwelch( data, fs )
            allp = [];
            for i = 1:size(data, 1 )
                p = pwelch( data(i,:), [],[],[],fs );
                allp = [ allp; p' ];
            end
            size(allp);
            figure(734); clf; plot( 20*log10(mean( allp, 1 )) ); title('Ave PSD');
            figure(735); clf; plot( 20*log10(std( allp, 1 ) ) ); title(' Std PSD');
            
        end
        
    end
    
end