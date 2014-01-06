classdef Depend < handle 
    
    methods ( Static )
        
        function SPMadd
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/spm8/');
            path(g,path );
        end
        
        function SPMremove
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/spm8/');
            rmpath(g);
        end
        
        function TFadd
            g = genpath('~/gc/lib/tftb-0.2/');
            path(g,path );
        end
        
        function TFremove
            g = genpath('~/gc/lib/tftb-0.2/');
            rmpath(g,path );
        end
        
        function MNEadd
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/mne/matlab/toolbox/');
            path(g,path);
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/mne/share/matlab/');
            path(g,path);
        end
        
        function MNEremove
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/mne/matlab/toolbox/');
            rmpath(g );
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/mne/share/matlab/');
            rmpath(g );
        end
        
        function FreeSurferAdd
            g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/free_surfer_matlab');
            path(g,path);
        end
        
        function FreeSurferRemove
             g = genpath('/synapse/labresources/AnalysisPipeline/thirdparty/free_surfer_matlab');
            rmpath(g );
        end
        
        function CCAadd
             g = genpath('~/gc/gc_seth_toolbox/');
            path(g,path);
        end
        
        function CCAremove
            g = genpath('~/gc/gc_seth_toolbox/');
            rmpath(g );
        end
        
        function SurfaceCodeAdd
            g = genpath('/synapse/labresources/AnalysisPipeline/fNIRS_pipeline/display/showhead/');
            path(g,path);
        end
    end
    
    
end