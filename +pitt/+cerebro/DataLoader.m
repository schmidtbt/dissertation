classdef DataLoader < handle
   
    %{
    
    freq6empty = pitt.cerebro.DataLoader.load6Empty();
    freq15empty = pitt.cerebro.DataLoader.load15Empty();
    freq15ab = pitt.cerebro.DataLoader.load15AB();
    freq6ab = pitt.cerebro.DataLoader.load6AB();
    
    distances = pitt.cerebro.DataLoader.loadDistances();
    
    [A,B,C] = pitt.cerebro.DataLoader.loadVertexIdxs();
    
    %}
    
    
    methods (Static)
        
        function freq = loadOptEmpty( pathToEmpty )
            
            freq = load( pathToEmpty );
            freq = freq.output;
            freq = pitt.diss.PLVAnalysis.full_map( freq );
            
        end
        
        % Use for loading generic all-all data
        % pathToRest is *.mat file of data
        function freq = loadOptRest( pathToRest )
            
            freq = load( pathToRest );
            
            % Backwards compatible output formats
            if( isfield( freq, 'freq' ) )
                freq = freq.freq;
            else
                freq = freq.output;
            end
            
            % Flip over diag to make symmetric
            freq = pitt.diss.PLVAnalysis.full_map( freq );
            
        end
        
        
        
        
        % --------------------------------
        %  Specific Data Loaders
        % --------------------------------
        function freq6empty = load6Empty()
            
            freq6empty = load( '~/data/empty/freq_006.mat' );
            freq6empty = freq6empty.output( 1:4098, 1:4098 );
            freq6empty = pitt.diss.PLVAnalysis.full_map( freq6empty );
            
        end
        
        function freq15empty = load15Empty()
            
            freq15empty = load( '~/data/empty/freq_015.mat' );
            freq15empty = freq15empty.output( 1:4098, 1:4098 );
            freq15empty = pitt.diss.PLVAnalysis.full_map( freq15empty );
            
        end
        
        function freq15ab = load15AB()
            
            freq15ab = load( '~/data/absim/freq_015.mat' );
            freq15ab = freq15ab.freq;
            freq15ab = pitt.diss.PLVAnalysis.full_map( freq15ab );
            
        end
        
        function freq6ab = load6AB()
            
            freq6ab = load( '~/data/absim/freq_006.mat' );
            freq6ab = freq6ab.freq;
            freq6ab = pitt.diss.PLVAnalysis.full_map( freq6ab );
            
        end
        
        function distances = loadDistances()
           
            distances = load( '~/data/distances.mat' );
            distances = distances.distances;
            
        end
        
        function [A,B,C] = loadVertexIdxs()
            A = 3673;
            B = 534;
            C = 2029;
        end
        
    end
    
end