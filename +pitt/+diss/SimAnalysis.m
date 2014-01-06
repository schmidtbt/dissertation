classdef SimAnalysis < handle
    
    methods (Static)
        %{
        [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
        [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
        %}
        function [v, f, lh_coords, lh_vertices] = load_init_data()
            
            base_path = '~/simu/';
            
            fwd_vert_data = load( [base_path,'fwd_v_f.mat']);
            
            v = fwd_vert_data.v;
            f = fwd_vert_data.f;
            
            distances = load( [base_path,'distances.mat']);
            distances = distances.distances;
            
            lh_coords = load( [base_path,'lh_coords.mat']);
            lh_coords = lh_coords.lh_coords;
            
            lh_vertices = load( [base_path,'lh_vertex.mat']);
            lh_vertices = double(lh_vertices.lh_vertex);
            
        end
        
        function new_lh_point_simulation( v, f, lh_coords, lh_vertices, sim_name, label_idxs, num_averages, end_time )
            
            base_path = '~/simu/';
            
            cd( base_path );
            
            if( exist( sim_name, 'dir' ) == 7 )
                error( 'Sim Directory Already Exists' );
            end
            mkdir( sim_name );
            
            local_path = [base_path, sim_name];
            cd( local_path );
            unix( 'cp ../base_diss/* .' );
            
            fid = fopen( 'sim.log', 'w' );
            
            fprintf( fid,'\n\nSimulation:\t%s\n\n', sim_name );
            fprintf( fid,'CWD: \t\t%s\n', pwd );
            fprintf( fid,'NumAve: \t%i\n', num_averages );
            fprintf( fid,'End-time: \t%i\n', end_time );
            for i = 1:length( label_idxs )
               fprintf( fid, 'Idx:\t\t%i\n', label_idxs(i) ); 
            end
            fprintf( fid, '\n\n\n' );
            
            simdef              = pitt.diss.SimDef;
            simdef.numaverages  = num_averages;
            simdef.end_time     = end_time;
            simdef.labelfile    = 'label-lh.label';
            
            save( 'simdef.mat', 'simdef' );
            
            figure(5); h = pitt.exp.simu.GraphAnalysis.displayLHBrain( v, f );
            pitt.exp.simu.GraphAnalysis.displayVertexLocation( h, label_idxs );
            print( 5, [sim_name,'-labels'], '-djpeg' );
            close 5;
            
            fprintf( fid, 'Writing new label: %s\n', simdef.labelfile );
            pitt.diss.Sim.write_new_label( simdef, simdef.labelfile, [lh_vertices(label_idxs,:), lh_coords(label_idxs,:)] );
            fprintf( fid, '%s\n', pitt.diss.Sim.forward_simulate( simdef ) );
            fprintf( fid, '%s\n', pitt.diss.Sim.process_event_timings( simdef, 'fwdsubj-sim.fif' ) );
            fprintf( fid, '%s\n', pitt.diss.Sim.inverse_simulate( simdef, 'fwdsubj-sim.fif',1 ) );
            
            fclose( fid );
            
        end
        
        %{
        
        Simulations To Date:
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'A_only_n1_sim01', 3673, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'A_only_n1_sim02', 3673, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'A_only_n1_sim03', 3673, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'A_only_n1_sim04', 3673, 1, 60000 );
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'B_only_n1_sim01', 534, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'B_only_n1_sim02', 534, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'B_only_n1_sim03', 534, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'B_only_n1_sim04', 534, 1, 60000 );
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AB_n1_sim01', [3673 534], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AB_n1_sim02', [3673 534], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AB_n1_sim03', [3673 534], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AB_n1_sim04', [3673 534], 1, 60000 );
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'C_n1_sim01', 2029, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'C_n1_sim02', 2029, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'C_n1_sim03', 2029, 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'C_n1_sim04', 2029, 1, 60000 );
        
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AC_n1_sim01', [3673 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AC_n1_sim02', [3673 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AC_n1_sim03', [3673 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'AC_n1_sim04', [3673 2029], 1, 60000 );
        
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'BC_n1_sim01', [534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'BC_n1_sim02', [534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'BC_n1_sim03', [534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'BC_n1_sim04', [534 2029], 1, 60000 );
        
        
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'ABC_n1_sim01', [3673 534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'ABC_n1_sim02', [3673 534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'ABC_n1_sim03', [3673 534 2029], 1, 60000 );
        pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, 'ABC_n1_sim04', [3673 534 2029], 1, 60000 );
        
        %}
        
        %{
        Assumes used: new_lh_point_simulation()
        
        Writes blocks to disk for PLV analysis
        %}
        function write_lh_point_simulation(sim_name)
            
            
            base_path = '~/simu/';
            
            cd( [base_path, sim_name] );
            
            pitt.Depend.MNEadd;
            lh_meg_data = mne_read_stc_file('fwdsubj-sim-cortex-lh.stc');
            
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvsim';
            subj_path   = 'point_sim';
            data_type   = [sim_name,'_raw'];
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        %{
        
        Writers Executed:
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'A_only_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'A_only_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'A_only_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'A_only_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'B_only_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'B_only_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'B_only_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'B_only_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AB_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AB_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AB_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AB_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'C_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'C_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'C_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'C_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AC_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AC_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AC_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'AC_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'BC_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'BC_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'BC_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'BC_n1_sim04' );
        
        
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'ABC_n1_sim01' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'ABC_n1_sim02' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'ABC_n1_sim03' );
        pitt.diss.SimAnalysis.write_lh_point_simulation( 'ABC_n1_sim04' );
        
        
        %}
        
        
        
        function plv_analysis_lh_point_simulation( sim_name, freqs )
            
            
            master = legion.stream.Master(2);
            
            grid_job        = ['job',sim_name,'_plv_calc'];
            for i = freqs;
               grid_job     = [grid_job,sprintf('_%i',(i)) ]; 
            end
            
            proj_path       = 'plvsim';
            subj_path       = 'point_sim';
            in_data_type    = [sim_name,'_raw'];
            out_data_type   = [sim_name,'_plv'];
            
            % Legion job name
            input_path      = sprintf( '%s/%s/%s', proj_path, subj_path, in_data_type ) ;
            output_path     = sprintf( '%s/%s/%s', proj_path, subj_path, out_data_type ) ;
            
            NumProcs        = 40;
            
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
        
        %{
        
        PLV calculations
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'A_only_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'A_only_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'A_only_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'A_only_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'B_only_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'B_only_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'B_only_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'B_only_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'C_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'C_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'C_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'C_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AB_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AB_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AB_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AB_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AC_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AC_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AC_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'AC_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'BC_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'BC_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'BC_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'BC_n1_sim04', [6 15] );
        
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'ABC_n1_sim01', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'ABC_n1_sim02', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'ABC_n1_sim03', [6 15] );
        pitt.diss.SimAnalysis.plv_analysis_lh_point_simulation( 'ABC_n1_sim04', [6 15] );
        
        %}
        
        % Method used for plv_post_analysis
        % Assumes a 1D cell string array ( a{1} = 'mystring' )
        function success = plv_post_analysis_analyze( path )
            path = path{1};
            class( path )
            
            fprintf( 'Entering: %s\n', path );
            cd( path );
            
            poss = dir( 'freq_*' );
            
            for i =1:length(poss)
                
                if( poss(i).isdir )
                    fprintf( 'Working on: %s\n', poss(i).name );
                    freq = legion.stream.Util.read_2d_split_as_array( ['./',poss(i).name], 100 );
                    save( [poss(i).name,'.mat'], 'freq' );
                    
                end
                
            end
            success = 1;
        end
        
        % Given a cell array of paths to completed _plv analysis
        % reintegrate all frequencies and save as freq_XXX.mat in that
        % directory
        function plv_post_analysis( cell_plv_data, jobName )
            
            
            kernel = legion.Kernel();
            kernel.add( @pitt.diss.SimAnalysis.plv_post_analysis_analyze, 'X' );
            
            master = legion.Master( kernel, cell_plv_data, length( cell_plv_data ), jobName );
            master.run();
            
        end
        
        %{
        
                
        % Find all _plv files in the current directory and use them to
        % regroup any found freq folders
        base = '/synapse/logs/schmidtb/plvsim/point_sim/';
        plv_loc = {};
        f = dir( '*_plv' );
        for i = 1:length(f)
            plv_loc{i} = [base, f(i).name];
        end
        plv_loc = plv_loc';
        pitt.diss.SimAnalysis.plv_post_analysis( plv_loc, 'job_plvsim_point_sim_regroup' );
        
        %}
        
        % Use: [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
        function success = do_basic_post_analysis( path )
            
            [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
            
            idxs = [3673, 563, 2029];
            
            class( path )
            
            fprintf( 'Entering: %s\n', path );
            cd( path );
            
            poss = dir( 'freq_*.mat' );
            
            currentDirectory = pwd;
            [upperPath, deepestFolder, ~] = fileparts(currentDirectory);
            % The name of the simulation
            simname = deepestFolder;
            
            
            for i =1:length(poss)
                freqmap_hz = poss(i).name;
                
                lst = find( freqmap_hz == '.' );
                freqmap_hs_name = freqmap_hz(1:lst-1);
                
                fullsimname = [simname, '-', freqmap_hs_name];
                
                fprintf( 'Working on: %s\n', fullsimname );
                
                % Load data
                freq = load( freqmap_hz );
                freq = freq.freq;
                
                % complete both halves
                freq = pitt.diss.PLVAnalysis.full_map( freq );
                
                figure(1);
                imagesc( freq );
                title( ['Phase map: ', fullsimname], 'interpreter', 'none' );
                xlabel( 'vertex' );
                ylabel( 'vertex' );
                
                print( 1, [fullsimname, '_phase_map' ], '-djpeg' );
                
                for j = 1:length(idxs)
                    figure(2);
                    cla;
                    h = pitt.exp.simu.GraphAnalysis.displayBrain( v, f );
                    view( -90, 0 );
                    
                    pitt.exp.simu.GraphAnalysis.overlayPLVData( h, freq( idxs(j), : ) );
                    caxis( [0 1] )
                    title( [fullsimname, '_seed_idx_', num2str(idxs(j)) ], 'interpreter', 'none' );
                    print( 2, [fullsimname, '_seed_idx_', num2str(idxs(j)) ], '-djpeg' );
                end
                
                
            end
            success = 1;
            
        end
        
        
        %{
        
        d = dir( '*_plv' );
        cwd = pwd;
        for i = 1:length(d)
            path = [cwd,'/',d(i).name];
            pitt.diss.SimAnalysis.do_basic_post_analysis( path );
        end
        
        %}
        
        
    end
    
end