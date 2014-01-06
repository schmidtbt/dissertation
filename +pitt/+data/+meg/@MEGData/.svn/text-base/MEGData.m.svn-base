classdef MEGData < pitt.data.Data
    
    properties
        datalocation
    end
    
    methods
        
        function obj = MEGData( datapath )
            obj.datalocation = datapath;
            D = spm_eeg_load( datapath )
            obj.data = D;
        end
        
        function dat = getMagData(obj)
            
            idx = []; 
            for i = 1:length(obj.data.chantype); 
                if(  strcmp(obj.data.chantype{i},'MEGMAG') ); 
                    idx = [idx; i]; 
                end; 
            end;
            
            dat = obj.data(idx, :,:);
        end
        
        function dat = getGradData(obj)
            
            idx = []; 
            for i = 1:length(obj.data.chantype); 
                if(  strcmp(obj.data.chantype{i},'MEGPLANAR') ); 
                    idx = [idx; i]; 
                end;
            end;
            
            dat = obj.data(idx, :,:);
        end
        
        function dat = getAllData(obj)
            dat = [ obj.getMagData(); obj.getGradData() ];
        end
        
        function X = getData( obj )
            X = obj.getAllData();
        end
        
    end
    
end