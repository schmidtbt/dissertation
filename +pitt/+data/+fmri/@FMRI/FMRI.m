classdef FMRI < pitt.data.Data
    
    properties
        sliceX;
        sliceY;
        sliceZ;
        timepts;
    end
    
    methods
        
        %Loads *.nii files of fMRI data
        function obj = FMRI( path )
            
            pitt.Depend.FreeSurferAdd;
            dobj = load_nifti( path );
            obj.sliceX = size( dobj.vol, 1 );
            obj.sliceY = size( dobj.vol, 2 );
            obj.sliceZ = size( dobj.vol, 3 );
            obj.timepts = size( dobj.vol, 4 );
            obj.data = reshape( dobj.vol, obj.sliceX*obj.sliceY*obj.sliceZ, obj.timepts );
            
        end
        
    end
    
    methods (Static)
        
        function fourD = loadFmriData( vertexByTime, sX, sY, sZ )
            fourD = reshape( vertexByTime, sX, sY, sZ );
        end
        
    end
    
end