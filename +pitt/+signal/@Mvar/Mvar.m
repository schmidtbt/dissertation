classdef Mvar
    
    
    methods( Static )
        
        function X = jointProcess( AR1, AR2, NumSamples )
            
            if( length( AR1 ) ~= length( AR2 ) )
                %error('AR1 and AR2 must have same length');
            end
            
            X = rand( 1, NumSamples );
            for i = length( AR1 )+1:NumSamples
                
                tX = X(i);
                for j = 1:length(AR1)
                    tX = tX + X(i-j)*AR1(j);
                end
                for j = 1:length(AR2)
                    tX = tX +X(i-j)*AR2(j);
                end
                X(i) = tX;
            end
            
            
            
        end
        
    end
    
end