classdef PinkGeneration < handle
    
    methods (Static)
        
        function signal = simulate( numARProcesses, numsamples )
            
            if( nargin == 2 )
                numsamples = 10000;
            end
                        
            coeff = rand(1, numARProcesses);
            signal = pitt.exp.PinkGeneration.simWithCoeff( coeff, numsamples );
        end
        
        function signal = simWithCoeff( coeff, numsamples )
            
            if( nargin == 2 )
                numsamples = 10000;
            end
            
            M = length(coeff);
            Xt = [];
            for i = 1:M
                r = randn( 1,numsamples );
                filt = filter( 1, [1 -coeff(i)], r );
                Xt = [Xt; filt];
            end
            signal = mean(Xt,1);
            
        end
        
    end
    
end