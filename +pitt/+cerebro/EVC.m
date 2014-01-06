classdef EVC < handle
    
    methods( Static )
        
        function [ V, D] = eigenVectors( freq )
            [V,D] = eig( freq );
        end
        
        % Calculate Eigen Vector centrality of all-by-all network
        % EVC values correspond to "Hubness"
        function evc = eigenVectorCentrality( freq )
            
            [V,D]   = pitt.cerebro.EVC.eigenVectors( freq );
            eigvals = D(1:size(D,1)+1:end);
            [M,I]   = max( eigvals );
            
            evc     = V(:,I);
        end
        
    end
    
end

