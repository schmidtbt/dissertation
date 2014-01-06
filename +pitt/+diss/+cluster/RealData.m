classdef RealData < handle
    
    methods (Static)
        
        function write_real_opt001_rest_2_lh_only()
            
            source = '~/data/Opt001_resting2_projon_ds4_filt_0_50-lh.stc';
            
            pitt.Depend.MNEadd;
            lh_meg_data = mne_read_stc_file( source );
            
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            sim_name = 'opt001_rest2';
            
            proj_path   = 'plvopt';
            subj_path   = 'opt001';
            data_type   = [sim_name,'_raw'];
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        function launch_plv_process_real_opt001_rest2_lh_only()
            
            freqs = [6 10 15 20];
            
            sim_name = 'opt001_rest2';
            
            master = legion.stream.Master(2);
            
            grid_job        = ['job',sim_name,'_plv_calc'];
            for i = freqs;
               grid_job     = [grid_job,sprintf('_%i',(i)) ]; 
            end
            
            proj_path       = 'plvopt';
            subj_path       = 'opt001';
            in_data_type    = [sim_name,'_raw'];
            out_data_type   = [sim_name,'_plv'];
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 6*8;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 150;
            Freqs           = freqs;
            
            % Complex inputs from Read kernel
            preproc_kern= pitt.exp.plv.Kernels.plvcomplexkernel(Fs, Freqs);
            master.setThreadKernel( preproc_kern );
            
            save_kern = pitt.exp.plv.Kernels.save_by_frequency_kernel_v2();
            master.setSaveKernel( save_kern );
            
            % Block level pre-processing
            read_kern = pitt.exp.plv.Kernels.read_and_complex_preprocess_kernel( Fs, Freqs);
            master.setReadKernel( read_kern );
            
            master.setGridWorkingDirectory( grid_job );
            master.setNumThreads( NumProcs );
            
            master.submit_job();
            
            
        end
        
    end
    
end