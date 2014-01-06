classdef Stream < handle
    
    methods (Static)
        
        % Iterate over opt* directories, and find freq maps. Then print them after imagesc()
        function [] = disp_freq_maps_dir_search()
            
            d = dir('opt*');
            
            for i=1:length(d)
                
                d2 = dir([d(i).name,'/rest*']);
                for j=1:length(d2)
                    
                    fprintf('Working On: %s\n', [d(i).name,'/',d2(j).name] );
                    
                    d3 = dir([d(i).name,'/',d2(j).name,'/*.mat']);
                    
                    for k=1:length(d3)
                        fprintf( '%s\n', [d(i).name,'/',d2(j).name,'/',d3(k).name] );
                        freq = pitt.cerebro.DataLoader.loadOptRest( [d(i).name,'/',d2(j).name,'/',d3(k).name] );
                        fprintf('Saving Image: %s\n', [d(i).name,'_',d2(j).name,'_',d3(k).name] );
                        figure(1); imagesc( freq );
                        title( [d(i).name,'_',d2(j).name,'_',d3(k).name],'interpreter','none' );
                        print(1,[d(i).name,'_',d2(j).name,'_',d3(k).name(1:8)], '-djpeg' );
                    end
                    
                end
                
            end
            
        end
        
        
    end
    
end

