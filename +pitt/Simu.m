classdef Simu < handle % pitt.Simu
    
    properties (Constant)
        path_mnebin = '/synapse/home/schmidtb/mne/bin/';
        path_base = '/synapse/home/schmidtb/simu/base/';
        
        start_time = 0; %start time for simulation (in ms)
        end_time = 600000; %end time for simulation (in ms)
        sfreq = 150;
    end
    
    % Execpts the necessary simulation setup files to be located in
    % path_base
    properties
        
        subjectname;
        
        numaverages = 4;
        
        fwdfile = [pitt.Simu.path_base, 'fwd.fif']; % Forward simulation file
        invfile = [pitt.Simu.path_base, 'inv.fif']; % Inverse simulation file
        covfile = [pitt.Simu.path_base, 'cov.fif']; % Covariance file
        evefile = [pitt.Simu.path_base, 'fake_events.eve']; % Artificial events
        avefile = [pitt.Simu.path_base, 'fake_ave.ave']; % Artificial average file
        
        subject_name;
        simname;
        
        path_sim_dir; % Path to this simulation data run
        
        label1file;
        label2file;
        
        time1file;
        time2file;
        
        lhdata;
        rhdata;
        
        dipole_loc_label1
        label1_data
        label1_avg
        
        dipole_loc_label2
        label2_data
        label2_avg
        
    end
    
    methods
        
        function obj = Simu( subject_name, path_sim_dir, label1, label2, time1file, time2file )
            obj.subject_name = subject_name;
            obj.simname = [obj.subject_name, '-sim.fif'];
            obj.path_sim_dir = path_sim_dir;
            %Labels should reside in the path_sim_dir
            obj.label1file = [pitt.Simu.path_base, label1];
            obj.label2file = [pitt.Simu.path_base, label2];
            obj.time1file = [obj.path_sim_dir,'/', time1file];
            obj.time2file = [obj.path_sim_dir,'/', time2file];
            pitt.Depend.MNEadd;
        end
        
        function simulate( obj )
            obj.enter_dir();
            obj.simu_data();
            obj.combine_data();
            obj.process_movie();
        end
        
        function enter_dir( obj )
            eval(sprintf('cd %s', obj.path_sim_dir));
        end
        
        % Use forward solution to simulate the time files into the brain.
        % Outputs two files with 1* and 2* prepended to them.
        function simu_data( obj )
            obj.enter_dir();
            
            % Simulate first data file
            eval(sprintf('unix(''%smne_simu --fwd %s --label %s --meg --out 1%s --sfreq %s --tmin %s --tmax %s --senscov %s --nave %s --raw --timecourse %s'');', ...
                obj.path_mnebin, obj.fwdfile, obj.label1file, obj.simname, num2str(pitt.Simu.sfreq),num2str(pitt.Simu.start_time), num2str(pitt.Simu.end_time), obj.covfile, num2str(obj.numaverages), obj.time1file));
            
            % Simulate second data file
            eval(sprintf('unix(''%smne_simu --fwd %s --label %s --meg --out 2%s --sfreq %s --tmin %s --tmax %s --senscov %s --nave %s --raw --timecourse %s'');', ...
                obj.path_mnebin, obj.fwdfile, obj.label2file, obj.simname, num2str(pitt.Simu.sfreq),num2str(pitt.Simu.start_time), num2str(pitt.Simu.end_time), obj.covfile, num2str(obj.numaverages), obj.time2file));
            
        end
        
        % Since simu can only simulate one timecourse at a time, we add the
        % two to create a single brain with both time-courses.
        function combine_data( obj )
            
            % First, read in the saved data from simu_data()
            eval(sprintf('raw1 = fiff_setup_read_raw(''1%s'');',obj.simname));
            data1 = fiff_read_raw_segment(raw1);
            
            eval(sprintf('raw2 = fiff_setup_read_raw(''2%s'');',obj.simname));
            data2 = fiff_read_raw_segment(raw2);
            
            % Combine the two data types
            data_sum=data1+data2;
            
            picks=1:raw1.info.nchan;
            [outfid,cals] = fiff_start_writing_raw(obj.simname,raw1.info,picks);
            fiff_write_raw_buffer(outfid,data_sum,cals);
            
            picks=1:raw1.info.nchan;
            [outfid,cals] = fiff_start_writing_raw(obj.simname,raw1.info,picks);
            fiff_write_raw_buffer(outfid,data_sum,cals);
            
        end
        
        function process_movie( obj )
            tic;
            eval(sprintf('unix(''%smne_process_raw --raw %s --filteroff --ave %s --saveavetag %s_ave'');',obj.path_mnebin, obj.simname, obj.avefile, obj.simname));
            
            eval(sprintf('unix(''%smne_make_movie --inv %s --meas %s_ave.fif --signed --set 1 --picknormalcomp --subject %s --stc %s.stc'');',obj.path_mnebin, obj.invfile, [obj.subject_name,'-sim', obj.subject_name,'-sim.fif'], ...
                obj.subject_name, obj.simname));
            toc;
        end
        
        function post_simulate( obj )
            obj.read_stc_file();
            obj.read_labels( '/synapse/home/schmidtb/subjects_dir/Opt065/label/' );
        end
        
        function read_stc_file( obj )
            
            obj.enter_dir();
            
            disp('Reading stc data');
            disp(sprintf('\tReading LH data'));
            eval(sprintf('lh_meg_data=mne_read_stc_file(''%s-lh.stc'');',obj.simname));
            %disp(sprintf('\tReading RH data'));
            obj.lhdata = lh_meg_data;
            disp(sprintf('\t***Skipping RH data****' ) );
            %eval(sprintf('rh_meg_data=mne_read_stc_file(''%s-rh.stc'');',obj.simname));
            %obj.rhdata = rh_meg_data;
            disp(sprintf('done reading stc data'));
            
        end
        
        function read_labels( obj, labeldir )
            cd( labeldir );
            % Label 1
            [obj.dipole_loc_label1, obj.label1_data, obj.label1_avg] = obj.read_single_label( obj.label1file );
            [obj.dipole_loc_label2, obj.label2_data, obj.label2_avg] = obj.read_single_label( obj.label2file );
        end
        
        function [dipole_loc_label, label_data, label_avg] = read_single_label( obj, label )
            
            if ~isempty(strfind(label,'lh'))
                vertex=obj.lhdata.vertices;
                vert_hem='lh';
            elseif ~isempty(strfind(label,'rh'))
                vertex=obj.rhdata.vertices;
                vert_hem='rh';
            else
                error('Label file did not contain lh or rh');
            end
            disp(sprintf('Reading inverse label: %s', label ));
            eval(sprintf('[labelloc,x,y,z,val] = inverse_read_label(''%s'');',label));
            
            n = 0;
            for k = 1:length(labelloc)
                if ~isempty(find(vertex==labelloc(k))) %#ok<EFIND>
                    n = n+1;
                    dipole_loc_label(n)=find(vertex==labelloc(k));
                end
            end
            
            if ~isempty(strfind(label,'lh'))
                label_data=obj.lhdata.data(dipole_loc_label,:);
            elseif ~isempty(strfind(label,'rh'))
                label_data=obj.rhdata.data(dipole_loc_label,:);
            end
            
            if size(label_data,1)>1
                label_avg=squeeze(mean(rectify_MEG_in_ROI(label_data)));
            else
                label_avg=label_data;
            end
            
        end
        
        
    end
    
    
    methods( Static )
       
        function lagsin_multiple_1_avg
            tic;
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_2', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        function lagsin_multiple_1_avg_2
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_3', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        
        function lagsin_multiple_1_avg_3
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_4', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        
        function lagsin_multiple_1_avg_4
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_5', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        
        
        
        function lagsin_multiple_patched_2_avg1
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_2_avg1/', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        function lagsin_multiple_patched_3_avg1
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_3_avg1/', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        function lagsin_multiple_patched_4_avg1
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_4_avg1', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1;
            sim.simulate();
            toc;
        end
        
        
        function lagsin_multiple_patched_2_avg3
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_2_avg3', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 3;
            sim.simulate();
            toc;
        end
        
        function lagsin_multiple_patched_3_avg3
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_3_avg3', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 3;
            sim.simulate();
            toc;
        end
        
        function lagsin_multiple_patched_4_avg3
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_4_avg3', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 3;
            sim.simulate();
            toc;
        end
        
        
        function lagsin_multiple_patched_4_avg1000
            tic; 
            pitt.Depend.MNEadd();
            sim = pitt.Simu( 'subj', '~/simu/lagsin_multiple_patched_4_avg1000', 'S1-lh.label', 'DLPF-lh.label', 'sinewave_timecourse1.txt', 'sinewave_timecourse2.txt' );
            sim.numaverages = 1000;
            sim.simulate();
            toc;
        end
        
        
    end
    
    
    
    
end