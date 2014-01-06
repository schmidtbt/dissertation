classdef SimJob < handle 
    
    
    methods( Static )
        
        %{ Ideal, no noise %}
        function testing_job()
            
            cd( '/synapse/home/schmidtb/simu/testing_diss_sim' );
            
            [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
            
            simdef = pitt.diss.SimDef;
            
            idx1 = 3670; % A
            idx2 = 560;  % B
            
            labelA = 'A-lh.label';
            labelB = 'B-lh.label';
            labelAB = 'AB-lh.label';
            
            simdef.numaverages  = 1000;
            simdef.end_time     = 60000;
            
            
            % Create Labels
            pitt.diss.Sim.write_new_label( simdef, labelA, [lh_vertices(idx1,:), lh_coords(idx1,:)] );
            pitt.diss.Sim.write_new_label( simdef, labelB, [lh_vertices(idx2,:), lh_coords(idx2,:)] );
            pitt.diss.Sim.write_new_label( simdef, labelAB, [lh_vertices([idx1 idx2],:), lh_coords([idx1 idx2],:)] );
            
            % Fwd Simulate A
            simdef.labelfile = labelA;
            pitt.diss.Sim.forward_simulate( simdef, 'A' );
            
            % Fwd Simulate B
            simdef.labelfile = labelB;
            pitt.diss.Sim.forward_simulate( simdef, 'B' );
            
            % Fwd Simulate AB
            simdef.labelfile = labelAB;
            pitt.diss.Sim.forward_simulate( simdef, 'AB' );
            
            
            % Process Add Combo File
            pitt.diss.Sim.combine_fwd_simulation( simdef, 'A', 'B', 'A_B_combo.fif' );
            pitt.diss.Sim.process_event_timings( simdef, 'A_B_combo.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'A_B_combo.fif' );
            
            % Process A only
            pitt.diss.Sim.process_event_timings( simdef, 'Asubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'Asubj-sim.fif' );
            
            % Process B only
            pitt.diss.Sim.process_event_timings( simdef, 'Bsubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'Bsubj-sim.fif' );
            
            % Process AB label combo
            pitt.diss.Sim.process_event_timings( simdef, 'ABsubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'ABsubj-sim.fif' );
            
        end
        
        %{
        
        Noise added
        
        A simulation was run using numave = 1 for label A
        saved as An1.
        
        %}
        function ab_with_noise()
            
            cd( '/synapse/home/schmidtb/simu/ab_with_noise' );
            
            [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
            
            simdef = pitt.diss.SimDef;
            
            idx1 = 3670; % A
            idx2 = 560;  % B
            
            labelA = 'A-lh.label';
            labelB = 'B-lh.label';
            labelAB = 'AB-lh.label';
            labelAcluster = 'Acluster-lh.label';
            simdef.numaverages  = 4;
            simdef.end_time     = 60000;
            
            
            % Create Labels
            pitt.diss.Sim.write_new_label( simdef, labelAcluster, [lh_vertices(I(1:15),:), lh_coords(I(1:15),:)] );
            pitt.diss.Sim.write_new_label( simdef, labelA, [lh_vertices(idx1,:), lh_coords(idx1,:)] );
            pitt.diss.Sim.write_new_label( simdef, labelB, [lh_vertices(idx2,:), lh_coords(idx2,:)] );
            pitt.diss.Sim.write_new_label( simdef, labelAB, [lh_vertices([idx1 idx2],:), lh_coords([idx1 idx2],:)] );
            
            % Fwd Simulate A
            simdef.labelfile = labelA;
            pitt.diss.Sim.forward_simulate( simdef, 'A' );
            
            % Fwd Simulate B
            simdef.labelfile = labelB;
            pitt.diss.Sim.forward_simulate( simdef, 'B' );
            
            % Cluster A simulate
            simdef.labelfile = labelAcluster;
            pitt.diss.Sim.forward_simulate( simdef, 'Acluster' );
            
            
            % Fwd Simulate AB
            simdef.labelfile = labelAB;
            pitt.diss.Sim.forward_simulate( simdef, 'AB' );
            
            % Process Add Combo File
            pitt.diss.Sim.combine_fwd_simulation( simdef, 'A', 'B', 'A_B_combo.fif' );
            pitt.diss.Sim.process_event_timings( simdef, 'A_B_combo.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'A_B_combo.fif' );
            
            % Process A only
            pitt.diss.Sim.process_event_timings( simdef, 'Asubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'Asubj-sim.fif' );
            
            % Process B only
            pitt.diss.Sim.process_event_timings( simdef, 'Bsubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'Bsubj-sim.fif' );
            
            %{not run yet}%
            % Process AB label combo
            pitt.diss.Sim.process_event_timings( simdef, 'ABsubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'ABsubj-sim.fif' );
            
            % Cluster simulation
            pitt.diss.Sim.process_event_timings( simdef, 'Aclustersubj-sim.fif' );
            pitt.diss.Sim.inverse_simulate( simdef, 'Aclustersubj-sim.fif' );
            
            
        end
        
        %--------------------------------------------------%
        %       Data Writers - No Noise
        %--------------------------------------------------%
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_A_only()
            
            
            lh_meg_data = mne_read_stc_file('Asubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_A_only';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_B_only()
            
            
            lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_B_only';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_A_B_combo()
            
            
            lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_A_B_combo';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        
        %--------------------------------------------------%
        %       PLV Calculation - No Noise
        %--------------------------------------------------%
        
        %{ executed 11/27/2012 %}
        function master = A_B_combo_plv()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_A_B_combo_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_A_B_combo';
            out_data_type   = 'raw_A_B_combo_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/27/2012 %}
        function master = A_only_plv()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_A_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_A_only';
            out_data_type   = 'raw_A_only_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/27/2012 %}
        function master = B_only_plv()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_B_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_B_only';
            out_data_type   = 'raw_B_only_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        
        
        
        %--------------------------------------------------%
        %       Data Writers - With Noise
        %--------------------------------------------------%
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_A_only_noise()
            
            
            lh_meg_data = mne_read_stc_file('Asubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_A_only_noise';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/29/2012 %}
        function testing_writing_plv_An1_only_noise()
            
            
            lh_meg_data = mne_read_stc_file('An1subj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_An1_only_noise';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_B_only_noise()
            
            
            lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_B_only_noise';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/27/2012 %}
        function testing_writing_plv_A_B_combo_noise()
            
            
            lh_meg_data = mne_read_stc_file('A_B_combo-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_A_B_combo_noise';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        %{ executed 11/29/2012 %}
        function testing_writing_plv_A_cluster_noise()
            
            
            lh_meg_data = mne_read_stc_file('Aclustersubj-sim-cortex-lh.stc');
            
            %tpts        = [1:size( lh_meg_data.data,2 )];
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'ABC';
            data_type   = 'raw_A_cluster_noise';
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
            %lh_meg_data = mne_read_stc_file('Bsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('ABsubj-sim-cortex-lh.stc');
            %lh_meg_data = mne_read_stc_file('A_B_combosubj-sim-cortex-lh.stc');
            
        end
        
        
        
        %--------------------------------------------------%
        %       PLV Calculation - With Noise
        %--------------------------------------------------%
        
        %{ executed 11/27/2012 %}
        function master = A_B_combo_plv_noise()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_A_B_combo_noise_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_A_B_combo_noise';
            out_data_type   = 'raw_A_B_combo_noise_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/27/2012 %}
        function master = A_only_plv_noise()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_A_noise_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_A_only_noise';
            out_data_type   = 'raw_A_only_noise_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/29/2012 %}
        function master = An1_only_plv_noise()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_An1_noise_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_An1_only_noise';
            out_data_type   = 'raw_An1_only_noise_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/29/2012 %}
        function master = A_cluster_only_plv_noise()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_A_cluster_noise_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_A_cluster_noise';
            out_data_type   = 'raw_A_cluster_noise_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{ executed 11/27/2012 %}
        function master = B_only_plv_noise()
            
            master = legion.stream.Master(2);
            
            grid_job        = 'job_B_noise_plv_calc_6_15';
            
            proj_path       = 'plvsim';
            subj_path       = 'ABC';
            in_data_type    = 'raw_B_only_noise';
            out_data_type   = 'raw_B_only_noise_plv';
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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