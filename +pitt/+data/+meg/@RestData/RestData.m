classdef RestData < pitt.data.meg.MEGData
    
    properties (Constant )
        dpath =  '/synapse/home/schmidtb/data/meg/';
    end
    
    methods
        
        function obj = RestData( datapath )
            path = [pitt.data.meg.RestData.dpath, datapath ];
            obj = obj@pitt.data.meg.MEGData( path );
        end
        
    end
    
    methods( Static )
        function showFiles
            dir( [pitt.data.RestData.dpath, '*mat'] );
        end
    end
    
    
end