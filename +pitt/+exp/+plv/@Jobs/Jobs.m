classdef Jobs < handle
    
    properties
    end
    
    methods (Static)
        
        
        
        %--------------------------------------------------%
        %       Non-Streaming
        %--------------------------------------------------%
        
        function master = plv_lagsin_multiple_patched_4_avg1000_all_all()
            
            % Timepoints
            tpts        = [1:9000];
            
            Fs          = 150;
            Freqs       = 5:40;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1000';
            
            % Number of cores to run on
            NCores      = 60;
            
            % Legion job name
            jobName     = 'plv_lagsin_multiple_patched_4_avg1000_all_all';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Construct the job data
            job_data    = sim.lhdata.data;
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.plv.Kernels.plvkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_empty_room_opt001_test( lhdat )
            
            % Timepoints
            tpts        = [1:15000];
            
            Fs          = 250;
            Freqs       = 5:40;
                        
            % Number of cores to run on
            NCores      = 2;
            
            % Legion job name
            jobName     = 'plv_opt001_test_input_complex';
            
            % Construct the job data
            job_data    = lhdat;
            
            kernel      = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            master.setPreProcKernel( pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs) );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_empty_room_opt001_lh_only()
            
            filename    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-lh.stc';
            
            % Timepoints
            tpts        = [1:15000];
            
            Fs          = 250;
            Freqs       = 5:40;
                        
            % Number of cores to run on
            NCores      = 60;
            
            % Read STC file
            fprintf( 'Reading STC file\n' );
            dat         = mne_read_stc_file( filename );
            
            % Legion job name
            jobName     = 'plv_opt001_lh_empty_room';
            
            % Construct the job data
            job_data    = dat.data(:,tpts);
            
            kernel      = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            master.setPreProcKernel( pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs) );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_opt001_lh_only()
            
            
            
            filename    = '~/data/Opt001_resting2_projon_ds4_filt_0_50-lh.stc';
            
            % Timepoints
            tpts        = [1:15000];
            
            Fs          = 250;
            Freqs       = 5:40;
            
            % Number of cores to run on
            NCores      = 60;
            
            % Read STC file
            fprintf( 'Reading STC file\n' );
            dat         = mne_read_stc_file( filename );
            
            % Legion job name
            jobName     = 'plv_opt001_lh';
            
            % Construct the job data
            job_data    = dat.data(:,tpts);
            
            kernel      = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            master.setPreProcKernel( pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs) );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_opt001_all()
            
            
            
            filenameL    = '~/data/Opt001_resting2_projon_ds4_filt_0_50-lh.stc';
            filenameR    = '~/data/Opt001_resting2_projon_ds4_filt_0_50-rh.stc';
            
            % Timepoints
            tpts        = [1:15000];
            
            Fs          = 250;
            Freqs       = 5:40;
            
            % Number of cores to run on
            NCores      = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            % Legion job name
            jobName     = 'plv_opt001_all';
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            kernel      = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            master.setPreProcKernel( pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs) );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_opt001_empty_room_all()
            
            
            
            filenameL    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-lh.stc';
            filenameR    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-rh.stc';
            
            % Timepoints
            tpts        = [1:15000];
            
            Fs          = 250;
            Freqs       = 5:40;
            
            % Number of cores to run on
            NCores      = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            % Legion job name
            jobName     = 'plv_opt001_empty_room_all';
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            kernel      = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            master.setPreProcKernel( pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs) );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = plv_opt001_empty_room_all_spmd()
            
            
            
            filenameL    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-lh.stc';
            filenameR    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-rh.stc';
            
            % Timepoints
            tpts        = [1:75000];
            
            Fs          = 250;
            Freqs       = 5:40;
            
            % Number of cores to run on
            NProcCores  = 30;
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            % Legion job name
            jobName     = 'plv_opt001_empty_room_whole_head_spmd';
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            proc_kernel = pitt.exp.plv.Kernels.plvcomplexkernel( Fs, Freqs );
            save_kernel = pitt.exp.plv.Kernels.save_by_frequency_kernel();
            preproc_kern= pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs);
            
            master      = legion.MasterPairWiseSPMD( proc_kernel, job_data, NProcCores, NDatSplits, jobName );
            
            master.setPreProcKernel( preproc_kern );
            master.setSaveKernel( save_kernel );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        
        
        
        %--------------------------------------------------%
        %       Test Streamining
        %--------------------------------------------------%
        
        function master = testing_streaming_plv_preproc()
            
            master = legion.stream.Master(1);
            
            master.setInputPath(  'opt001_emptyroom_lh_raw' ); 
            master.setOutputPath( 'opt001_emptyroom_lh_raw_proc' );
            
            Fs          = 250;
            Freqs       = 5:40;
            
            preproc_kern= pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs);
            master.setThreadKernel( preproc_kern );
            
            master.setGridWorkingDirectory( 'grid/proc_opt001_emptyroom_lh' );
            master.setNumThreads( 5 );
            
            master.submit_job();
            
            
        end
        
        function master = testing_streaming_plv_calc()
            
            Fs          = 250;
            Freqs       = 5:40;
            
            master = legion.stream.Master(2); 
            
            proj_path       = 'plvdata';
            subj_path       = 'testing';
            in_data_type    = 'test_read_preproc';
            out_data_type   = 'test_read_preproc_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 9;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            % Complex inputs
            preproc_kern= pitt.exp.plv.Kernels.plvcomplexkernel(Fs, Freqs);
            master.setThreadKernel( preproc_kern );
            
            save_kern = pitt.exp.plv.Kernels.save_by_frequency_kernel_v2();
            master.setSaveKernel( save_kern );
            
            % Block level pre-processing
            read_kern = pitt.exp.plv.Kernels.read_and_complex_preprocess_kernel( Fs, Freqs);
            master.setReadKernel( read_kern );
            
            master.setGridWorkingDirectory( 'job_test_preproc_plv' );
            master.setNumThreads( NumProcs );
            
            master.submit_job();
            
        end
        
        function master = testing_streaming_plv_no_preproc_calc()
            
            Fs          = 250;
            Freqs       = 5:40;
            
            master = legion.stream.Master(2); 
            
            proj_path       = 'plvdata';
            subj_path       = 'testing';
            in_data_type    = 'test_read_preproc';
            out_data_type   = 'test_plv_nopreproc';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 9;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            % Raw inputs
            preproc_kern= pitt.exp.plv.Kernels.plvkernel(Fs, Freqs);
            master.setThreadKernel( preproc_kern );
            
            save_kern = pitt.exp.plv.Kernels.save_by_frequency_kernel_v2();
            master.setSaveKernel( save_kern );
            
            % No preprocessing during reading
            %read_kern = pitt.exp.plv.Kernels.read_and_complex_preprocess_kernel( Fs, Freqs);
            %master.setReadKernel( read_kern );
            
            master.setGridWorkingDirectory( 'job_test_no_preproc_plv' );
            master.setNumThreads( NumProcs );
            
            master.submit_job();
            
        end
        
        
        
        
        %--------------------------------------------------%
        %       File Writers (Single-threaded)
        %--------------------------------------------------%
        
        % --- Opt001 real data
        function write_raw_file_opt001_empty_room()
            
            filenameL    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-lh.stc';
            filenameR    = '~/data/Opt001_emptyroom_projon_ds4_filt_0_50-rh.stc';
            
            % Timepoints
            tpts        = [1:75000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'opt001';
            data_type   = 'raw_all_empty_room';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        function write_raw_file_opt001_rest2()
            
            filenameL    = '~/data/Opt001_resting2_projon_ds4_filt_0_50-lh.stc';
            filenameR    = '~/data/Opt001_resting2_projon_ds4_filt_0_50-rh.stc';
            
            % Timepoints
            tpts        = [1:75000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'opt001';
            data_type   = 'raw_all_rest2';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        % --- Simu of plv data -- but used BAD labels
        function write_raw_file_avniel_clone_simu()
            
            filenameL    = '~/simu/avniel_clone/subj-sim.fif-lh.stc';
            filenameR    = '~/simu/avniel_clone/subj-sim.fif-rh.stc';
            
            % Timepoints (5 mins of data)
            tpts        = [1:45000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'simu';
            data_type   = 'dual_s1_locked_raw';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        function write_raw_file_avniel_clone_empty_room_simu()
            
            filenameL    = '~/simu/avniel_clone_empty_room/subj-sim.fif-lh.stc';
            filenameR    = '~/simu/avniel_clone_empty_room/subj-sim.fif-rh.stc';
            
            % Timepoints (5 mins of data)
            tpts        = [1:45000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'simu';
            data_type   = 'empty_room_raw';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        
        % --- Files from simu using opt001 sufrace and labels
        function write_raw_file_simu_opt001_empty_room_simu()
            
            filenameL    = '~/simu/opt001_empty_sim/subj-sim.fif-lh.stc';
            filenameR    = '~/simu/opt001_empty_sim/subj-sim.fif-rh.stc';
            
            % Timepoints (5 mins of data)
            tpts        = [1:45000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'simu';
            data_type   = 'opt001_empty_room_raw';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        % --- Files from simu using opt001 sufrace and labels
        function write_raw_file_simu_opt001_dual_s1_simu()
            
            filenameL    = '~/simu/avniel_clone_empty_room/subj-sim.fif-lh.stc';
            filenameR    = '~/simu/avniel_clone_empty_room/subj-sim.fif-rh.stc';
            
            % Timepoints (5 mins of data)
            tpts        = [1:45000];
            
            NDatSplits  = 240;
            
            % Read STC file
            fprintf( 'Reading STC file LH\n' );
            datL         = mne_read_stc_file( filenameL );
            fprintf( 'Reading STC file RH\n' );
            datR         = mne_read_stc_file( filenameR );
            
            proj_path   = 'plvdata';
            subj_path   = 'simu';
            data_type   = 'opt001_dual_s1_raw';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = [datL.data(:,tpts); datR.data(:,tpts)];
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        
        
        
        %--------------------------------------------------%
        %       Subject Processing
        %--------------------------------------------------%
        function master = pre_process_to_complex_opt001_rest2()
            
            master = legion.stream.Master(1);
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_rest2';
            out_data_type   = 'raw_all_rest2_pre_processed';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = 5:40;
            
            preproc_kern    = pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs);
            master.setThreadKernel( preproc_kern );
            
            master.setGridWorkingDirectory( 'opt001_rest2_pre_processing' );
            master.setNumThreads( NumProcs );
            
            master.submit_job();
            
        end
        
        function master = pre_process_to_complex_opt001_empty_room()
            
            master = legion.stream.Master(1);
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_empty_room';
            out_data_type   = 'raw_all_empty_room_pre_processed_by_freq';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = 5:40;
            
            preproc_kern    = pitt.exp.plv.Kernels.preprocess_to_complex(Fs, Freqs);
            save_kern       = pitt.exp.plv.Kernels.save_by_frequency_kernel_v2();
            
            master.setThreadKernel( preproc_kern );
            master.setSaveKernel( save_kern );
            
            
            master.setGridWorkingDirectory( 'opt001_emptyroom_pre_processing' );
            master.setNumThreads( NumProcs );
            
            master.submit_job();
            
        end
        
        function master = opt001_empty_room_plv()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_opt001_empty_room_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_empty_room';
            out_data_type   = 'empty_room_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = [6, 15];
            
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
        
        function master = opt001_empty_room_plv_10_25hz()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_opt001_empty_room_plv_calc_10_25';
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_empty_room';
            out_data_type   = 'empty_room_plv_10_25';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = [10 25];
            
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
        
        function master = opt001_rest2_plv()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_opt001_rest2_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_rest2';
            out_data_type   = 'rest2_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = [6, 15];
            
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
        
        function master = opt001_rest2_plv_10_25hz()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_opt001_rest2_plv_calc_10_25';
            
            proj_path       = 'plvdata';
            subj_path       = 'opt001';
            in_data_type    = 'raw_all_rest2';
            out_data_type   = 'rest2_plv_10_25';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 250;
            Freqs           = [10 25];
            
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
        
        
        
        
        %--------------------------------------------------%
        %       Simu Processing
        %--------------------------------------------------%
        function master = simu_dual_s1_plv_calc()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_simu_dual_s1_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'simu';
            in_data_type    = 'dual_s1_locked_raw';
            out_data_type   = 'dual_s1_locked_raw_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 150;
            Freqs           = [6, 15];
            
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
        
        function master = simu_empty_room_plv_calc()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_simu_empty_room_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'simu';
            in_data_type    = 'empty_room_raw';
            out_data_type   = 'empty_room_raw_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 150;
            Freqs           = [6, 15];
            
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
        
        function master = simu_opt001_empty_room_plv_calc()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_simu_opt001_empty_room_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'simu';
            in_data_type    = 'opt001_empty_room_raw';
            out_data_type   = 'opt001_empty_room_raw_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 150;
            Freqs           = [6, 15];
            
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
        
        function master = simu_opt001_dual_s1_plv_calc()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_simu_opt001_dual_s1_plv_calc_6_15';
            
            proj_path       = 'plvdata';
            subj_path       = 'simu';
            in_data_type    = 'opt001_dual_s1_raw';
            out_data_type   = 'opt001_dual_s1_raw_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 80;
            
            master.setInputPath(  input_path ); 
            master.setOutputPath( output_path );
            
            Fs              = 150;
            Freqs           = [6, 15];
            
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

