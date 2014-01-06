classdef Kernels < handle
   
    methods (Static)
        
        %--------------------------------------------------%
        %       PLV Calculation
        %--------------------------------------------------%
        
        % For raw input
        function kernel = plvkernel( Fs, Freqs)
            
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.PLV.resting_plv_pairwise_stacked, 'X', Fs, Freqs );
            
        end
        
        % For complex input
        function kernel = plvcomplexkernel( Fs,Freqs )
           
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.PLV.resting_plv_pairwise_complex_stacked, 'X', Fs, Freqs );
            
        end
        
        
        %--------------------------------------------------%
        %       Preprocessing
        %--------------------------------------------------%
        
        function kernel = preprocess_to_complex_single_vertex( Fs, Freqs )
           
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.PLV.preprocess_single_vertex, 'X', Fs, Freqs );
            
        end
        
        function kernel = preprocess_to_complex_block( Fs, Freqs )
           
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.PLV.preprocess_data, 'X', Fs, Freqs );
            
        end
        
        
        %--------------------------------------------------%
        %       Save By Freq Kernels
        %--------------------------------------------------%
        
        function kernel = save_by_frequency_kernel()
            
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.Kernels.splitOutputByFrequency, 'X' );
            
        end
        
        function result = splitOutputByFrequency( input )
            
            dataPath        = input.dataPath;
            jrID            = input.jrID;
            r_block_idx     = input.r_block_idx;
            c_block_idx     = input.c_block_idx;
            orig_output     = input.output;
            rowidxs         = input.rowidxs;
            colidxs         = input.colidxs;
            
            
            numFreq         = length( orig_output{1,1} );
            
            for k = 1:numFreq
                
                intermed    = -1.*ones( size(orig_output,1), size(orig_output,2) );
                
                for i = 1:size(orig_output,1)
                    for j = 1:size(orig_output,2)
                        if( orig_output{i,j}{k} )
                            intermed( i,j ) = orig_output{i,j}{k};
                        end
                            
                    end
                end
                
                output = intermed;
                
                freq_folder = [dataPath,'freq_',num2str(k)];
                [SUCCESS,MESSAGE,MESSAGEID] = mkdir( freq_folder );
                
                filename = [freq_folder,'/output',num2str( r_block_idx ), '_', num2str( c_block_idx ),'.mat'];
                fprintf( 'Saving: %s\n', filename );
                freq = k;
                save( filename, 'output', 'r_block_idx', 'c_block_idx', 'jrID', 'rowidxs', 'colidxs', 'freq' , '-v7.3' );
            end
            
            result = input;
            
        end
        
        function kernel = save_by_frequency_kernel_v2()
            
            kernel = legion.Kernel;
            kernel.add( @pitt.exp.plv.Kernels.splitOutputByFrequencyUsingDefault, 'X' );
            
        end
        
        function result = splitOutputByFrequencyUsingDefault( input )
            
            dataPath        = input.dataPath;
            jrID            = input.jrID;
            r_block_idx     = input.r_block_idx;
            c_block_idx     = input.c_block_idx;
            orig_output     = input.split;
            rowidxs         = input.rowidxs;
            colidxs         = input.colidxs;
            
            numFreq         = length( orig_output{1,1} );
            
            for k = 1:numFreq
                
                intermed    = -1.*ones( size(orig_output,1), size(orig_output,2) );
                
                for i = 1:size(orig_output,1)
                    for j = 1:size(orig_output,2)
                        if( ~isempty(orig_output{i,j}{k}) )
                            intermed( i,j ) = orig_output{i,j}{k};
                        else
                            doSkip = 1;
                        end 
                    end
                end
                
                % Skip empty frequencies
                if( doSkip )
                    doSkip = 0;
                    continue;
                end
                
                
                freq_folder                 = [dataPath,sprintf('/freq_%03i',k)];
                [SUCCESS,MESSAGE,MESSAGEID] = mkdir( legion.stream.Master.normPath(freq_folder) );
                
                freq_save_struct            = input;
                freq_save_struct.split      = intermed;
                freq_save_struct.dataPath   = freq_folder;
                
                save_kernel = legion.stream.Kernel.default_save_kernel();
                save_kernel.initial( freq_save_struct );
                save_kernel.execute();
                
            end
            
            result = input;
            
        end
        
        %--------------------------------------------------%
        %       Block Readers With Complex Pre-Processing
        %--------------------------------------------------%
        
        function kernel = read_and_complex_preprocess_kernel( Fs, Freqs )
            
            kernel = legion.Kernel();
            kernel.add( @pitt.exp.plv.Kernels.readAndComplexPreProcess, 'X', Fs, Freqs );
            
        end
        
        function result = readAndComplexPreProcess( filename, Fs, Freqs )
            
            r_kernel = legion.stream.Kernel.default_read_kernel();
            r_kernel.initial( filename );
            result = r_kernel.execute();
            
            p_kernel = pitt.exp.plv.Kernels.preprocess_to_complex_block( Fs, Freqs );
            p_kernel.initial( result.split );
            result.split = p_kernel.execute();
            
        end
        
        
    end
    
end