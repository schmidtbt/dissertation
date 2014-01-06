classdef SimulateABWithDist < handle
    
    methods(Static)
        
        
        %{
        [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
        [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
        %}
        
        
        %-------------------------------------------------------------------
        %       Simulation of A+B with a distance dependence
        %-------------------------------------------------------------------
        
        
        function names = simulate_LH_point_With_Noise()
            
            % Load some external data
            [v, f, lh_coords, lh_vertices] = pitt.diss.SimAnalysis.load_init_data();
            [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = pitt.diss.Sim.load_vertex_info();
            path_idxs_2d = pitt.diss.SimulateABWithDist.loadABPath();
            
            A = 3673;
            B = 534;
            
            % Append A to start simulations from A->B (not including B)
            path_idxs_2d = [A; path_idxs_2d];
            
            % Run this many at each Path->B pair
            simulations_per_AB_pair = 10;
            
            
            % Simulation params:
            numAverages = 1;
            end_time    = 60000;
            
            names = [];
            % Loop over each path point
            for i = 1:length(path_idxs_2d)
                
                % Label is path point + B
                label_idxs = [path_idxs_2d(i), B]
                
                % Loop over number of simulations at each pair
                for j = 1:simulations_per_AB_pair
                    
                    sim_name = sprintf( 'AB_slide_lh_only_tar_%04i_orig_%04i_sim_%02i', B, path_idxs_2d(i), j);
                    pitt.diss.SimAnalysis.new_lh_point_simulation( v, f, lh_coords, lh_vertices, sim_name, label_idxs, numAverages, end_time );
                    names = [names; sim_name];
                end
                
            end
            
            
            
        end
        
        
        function plv_write_and_launch_analyze( names )
            
            Freqs = [6 15];
            
            for i = 1:size(names,1)
                
                pitt.diss.SimulateABWithDist.write_lh_point_simulation( names(i,:) );
                pitt.diss.SimulateABWithDist.plv_analysis_lh_point_simulation( names(i,:), Freqs );
                
            end
            
        end
        
        
        function write_lh_point_simulation(sim_name)
            
            
            base_path = '~/simu/';
            
            cd( [base_path, sim_name] );
            
            pitt.Depend.MNEadd;
            lh_meg_data = mne_read_stc_file('fwdsubj-sim-cortex-lh.stc');
            
            data        = lh_meg_data.data;
            
            NDatSplits  = 100;
            
            proj_path   = 'plvabdist';
            subj_path   = 'absim';
            data_type   = [sim_name,'_raw'];
            
            % Legion job name
            output_path  = sprintf( '%s/%s/%s', proj_path, subj_path, data_type ) ;
            
            % Construct the job data
            job_data    = data;
            
            fprintf( 'Writing' );
            legion.stream.Util.write_split_data( job_data, output_path, NDatSplits );
            
        end
        
        
        function plv_analysis_lh_point_simulation( sim_name, freqs )
            
            
            master = legion.stream.Master(2);
            
            grid_job        = ['job',sim_name,'_plv_calc'];
            for i = freqs;
               grid_job     = [grid_job,sprintf('_%i',(i)) ]; 
            end
            
            proj_path       = 'plvabdist';
            subj_path       = 'absim';
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
        
        % names is a character array (it is converted to cell array first)
        % Launches onto SGE a reintegration scheme for every file
        function name_cell_arr = distr_reintegration( names )
            
            name_cell_arr = cell(0);
            
            for i = 1:size(names,1)
                
                proj_path       = 'plvabdist';
                subj_path       = 'absim';
                out_data_type   = [names(i,:),'_plv'];
                
                base_path       = '/synapse/logs/schmidtb';
                
                name_cell_arr{i}   = sprintf( '%s/%s/%s/%s', base_path, proj_path, subj_path, out_data_type );
                
            end
            
            pitt.diss.SimAnalysis.plv_post_analysis( name_cell_arr', 'job_sim_ab_with_dist');
            
        end
        
        
        
        
        %-------------------------------------------------------------------
        %       Path Calculation in PLV space
        %-------------------------------------------------------------------
        
        % path_idxs_2d = pitt.diss.SimulateABWithDist.loadABPath();
        function path_A_B_idxs = loadABPath()
            
            data = load( '~/simu/A_B_path_idxs.mat' );
            path_A_B_idxs = data.path_idxs_2d;
            
        end
        
        % Return the idxs which form the shortest path from origin to
        % target
        % Calls shortestDistance()
        function path_idxs_2d = shorestPath( distances, origin_idx, target_idx, radius )
            
            path_idxs_2d = [];
            
            cont = 1;
            
            currIdx = origin_idx;
            counter = 0;
            while( cont )
                
                next_idx = pitt.diss.SimulateABWithDist.shortestDistance( distances, currIdx, target_idx, radius );
                
                if( next_idx == target_idx )
                    cont = 0; % Break loop
                    break;
                end
                
                if( next_idx == currIdx )
                    error('Stuck'); 
                end
                
                path_idxs_2d = [path_idxs_2d; next_idx];
                currIdx = next_idx;
                
                counter = counter+1;
                
                if( counter > 4098 )
                    path_idxs_2d
                    error('Loop exceed total points');
                    break;
                end
            end
            
            
        end
        
        % Perform step-wise shortest distance calculation
        %{
        Assume regular grid. Given two points or origin and target, find
        next idx point which is shortest step from origin to target
        %}
        function next_idx = shortestDistance( distances, origin_idx, target_idx, radius )
            
            % Narrow local search
            [S,Isort] = sort( distances( origin_idx,: ) );
            search_idxs = Isort(1:radius);
            
            % Find closest local point to target
            [V,Imin] = min(distances( search_idxs, target_idx ) );
            
            % Reconstruct original idx
            next_idx = Isort(Imin);
            
        end
        
        
        function disp_path_vertices( h, path_idxs_2d )
            
            vd = zeros( 4098,1 );
            vd( path_idxs_2d ) = 1;
            pitt.exp.simu.GraphAnalysis.overlayPLVData( h, vd );
            
        end
        
    end
    
end