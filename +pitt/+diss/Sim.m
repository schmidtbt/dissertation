classdef Sim < handle
   
    properties (Constant)
        
        mne_bin_path = '/synapse/home/schmidtb/mne/bin/';
        
    end
    
    
    methods (Static)
        
        
        %____________________________________
        %       Simulation Pathway
        %____________________________________
        
        % Simulate equation based labels through forward model onto sensors
        function exec = forward_simulate( simdef, fwd_name )
            
            if( nargin == 1 )
                fwd_name = 'fwd';
            end
            
            cd( simdef.working_dir );
            
            exec = sprintf('unix(''%smne_simu --fwd %s --label %s --meg --out %s%s --sfreq %s --tmin %s --tmax %s --senscov %s --nave %s --raw --timecourse %s'');', ...
                pitt.diss.Sim.mne_bin_path, simdef.fwdfile, simdef.labelfile, fwd_name, simdef.simname, num2str(simdef.sfreq),num2str(simdef.start_time), num2str(simdef.end_time), simdef.covfile, num2str(simdef.numaverages), simdef.timefile);
            disp( exec );
            fprintf( '\nCreating output forward simulation file: %s%s \n', fwd_name, simdef.simname);
            eval( exec );
        end
        
        % Combine multiple sensor-space simulations
        %{
        fwd-sim_name should be a *.fif file
        %}
        function combine_fwd_simulation( simdef, fwd_name1, fwd_name2, fwd_sim_name )
            cd( simdef.working_dir );
            
            % Read data from different foward simulations
            [data1, raw1] = pitt.diss.Sim.read_fwd_data( simdef, [fwd_name1, simdef.simname] );
            [data2, raw2] = pitt.diss.Sim.read_fwd_data( simdef, [fwd_name2, simdef.simname] );
            
            combined_data = data1 + data2;
            
            % Now write data to simdef.simname *.fif file
            pitt.diss.Sim.write_sensor_file( simdef, combined_data, raw1, fwd_sim_name );
            
        end
        
        % Pre-process event timings for inverse simulation
        function exec = process_event_timings( simdef, fwd_sim_name )
            cd( simdef.working_dir );
            
            % Generate average file with correct timings
            pitt.diss.Sim.write_fake_ave_file( simdef )
            
            exec = sprintf('unix(''%smne_process_raw --raw %s --filteroff --ave %s --saveavetag %s_ave'');',pitt.diss.Sim.mne_bin_path, fwd_sim_name, simdef.avefile, simdef.simname);
            disp( exec );
            eval( exec );
            
            
        end
        
        % Project sensor-space calculations onto cortical surface over time
        function exec = inverse_simulate( simdef, fwd_sim_name, inbackground )
            
            if( nargin == 2 )
                inbackground = 0;
            end
            
            cd( simdef.working_dir );
            
            lst = find( fwd_sim_name == '.' );
            fwd_sim_name = fwd_sim_name( 1:lst-1 );
            
            exec = sprintf('%smne_make_movie --inv %s --meas %s_ave.fif --signed --set 1 --picknormalcomp --subject %s --stc %s.stc',pitt.diss.Sim.mne_bin_path, simdef.invfile, [fwd_sim_name,simdef.subject_name,'-sim.fif'], ...
                simdef.subject_name, [fwd_sim_name, '-cortex']);
            disp( exec );
            
            if( inbackground )
                
                fid = fopen('sim.tcsh', 'w' );
                
                fprintf( fid, '#Execute simulation\n' );
                %fprintf( fid, '. ~/bin/setup_mne\n' );
                fprintf( fid, 'hostname\n', exec );
                fprintf( fid, 'cd %s\n', pwd );
                fprintf( fid, 'sh ./bgsim.sh\n' );
                
                fclose( fid );
                
                system( 'chmod u+x sim.sh' );
                
                fid = fopen('bgsim.sh', 'w' );
                
                fprintf( fid, '#!/bin/sh\n' );
                fprintf( fid, '#Matlab blows. Workaround\n' );
                fprintf( fid, 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/synapse/home/schmidtb/mne/lib\nexport MATLAB_ROOT=/usr/local/MATLAB/R2011a\nexport MNE_ROOT=/synapse/home/schmidtb/mne\nexport PATH=$PATH:/synapse/home/schmidtb/mne/bin\nexport SUBJECTS_DIR=~/subjects_dir/\nexport PATH=$PATH:$MNE_ROOT/lib\n' );
                fprintf( fid, '%s\n', exec );
                
                fclose( fid );
                
                system( 'chmod u+x sim.tcsh' );
                system( 'chmod u+x bgsim.sh' );
                %[stat, res] = system( './bgsim.sh' );
                [stat, res] = system( 'qsub -wd `pwd` sim.tcsh' );
                res
                
            end
            
            %eval( exec );
            
        end
        
        
        
        %____________________________________
        %       UTLIITIES
        %____________________________________
        
        function [data, raw] = read_fwd_data( simdef, name )
            cd( simdef.working_dir );
            
            exec = sprintf('raw = fiff_setup_read_raw(''%s'');', name);
            disp( exec );
            eval( exec );
            disp( '' );
            data = fiff_read_raw_segment(raw);
        end
        
        
        function write_sensor_file( simdef, data, raw, newname )
            cd( simdef.working_dir );
            
            picks           = 1:raw.info.nchan;
            [outfid,cals]   = fiff_start_writing_raw( newname, raw.info,picks );
            
            fprintf( '\nWriting new sensor file: %s\n', newname );
            fiff_write_raw_buffer( outfid, data, cals );
            
        end
        
        function write_fake_ave_file( simdef )
            
            fid = fopen( simdef.avefile , 'w' );
            
            fprintf( fid, 'average {\n\tname\t\t"fake ave"\n\teventfile\t%s\n\tcategory {\n\t\tname\t"fake"\n\t\tevent\t1\n\t\ttmin\t%i\n\t\ttmax\t%.3f\n}\n', simdef.evefile, simdef.start_time, (simdef.end_time-1)/1000 );
            
        end
        
        
        %____________________________________
        %       LABEL BUILDING
        %____________________________________
        function [lh_vertices, rh_vertices, lh_coords, rh_coords, distances] = load_vertex_info()
            
            base_path = '~/simu/';
            
            lh_vertices = load( [base_path,'lh_vertex.mat']);
            rh_vertices = load( [base_path,'rh_vertex.mat']);
            lh_coords = load( [base_path,'lh_coords.mat']);
            rh_coords = load( [base_path,'rh_coords.mat']);
            distances = load( [base_path,'distances.mat']);
            
            lh_vertices = double(lh_vertices.lh_vertex);
            rh_vertices = double(rh_vertices.rh_vertex);
            lh_coords = lh_coords.lh_coords;
            rh_coords = rh_coords.rh_coords;
            distances = distances.distances;
            
            %vertex_idxs = [lh_vertices;rh_vertices];
        end
        
        %{
        
        label_name should end with -rh.label or -lh.label
        
        vertex_data is [mne_idx, x, y, z]
        
        %}
        function write_new_label( simdef, label_name, vertex_data )
            cd( simdef.working_dir );
            
            fid = fopen( label_name, 'w' );
            fprintf( fid, '#generated by matlab, cause MNE is too bloated\n' );
            fprintf( fid, '%i\n', size( vertex_data, 1 ) );
            for i = 1:size( vertex_data,1 )
                fprintf( fid, '%i %f %f %f 0\n', vertex_data(i,1), vertex_data(i,2), vertex_data(i,3), vertex_data(i,4) );
            end
            fprintf( fid, '\n' );
            fclose( fid );
            
        end
        
    end
    
    
    
    
end