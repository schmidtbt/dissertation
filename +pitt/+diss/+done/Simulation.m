classdef Simulation < handle
    
    methods (Static)
        
        
        %{
        False positive as a function of distance from given index
        
        Radial fall off of the plvMap(A,:) from the idx seed point
        
        %}
        function fp_from_idx = falsePositiveAroundB( cluster_idxs, distances, plv_data_1d, idx )
            
            
        end
        
        %{
        
        Calculate the confusion matrix (TP/FP etc) for the given data (plvMap(idx1,:))
        eminanting from idx1 (usually A) and idx2 (usually B)
        
        %}
        function confusion_matrix = falsePositiveRateForDistance( cluster_idxs, plv_data_1d, idx1, idx2 )
            
        end
        
        
        function tri_calc = triangulateAndNormalizeData( cluster_idxs, plv_data_2d, plv_empty_2d )
            
            
        end
        
    end
    
    
end