classdef Jobs < handle
   % Define specific jobs run on the system. and maintain their code as a
   % log of what was done
    methods( Static )
        
        function master = lagsin_multiple_patched_4_avg1_5lags_labels_and_random_100pts()
            
            % Timepoints
            tpts        = [1000:1100];
            
            % Lags
            NLags       = 5;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1';
            
            % Number of cores to run on
            NCores      = 20;
            
            % Legion job name
            jobName     = 'lagsin_multiple_patched_4_avg1_5lags_labels_and_random_100pts';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Extract label data
            [l_data, l_vertices] = pitt.exp.simu.Processing.extractLabelDataAllAll( sim );
            
            % Extract random vertices
            [r_data, r_vertices] = pitt.exp.simu.Processing.getRandomVertices( sim, size( l_data,1 ) );
            
            % Construct the job data
            job_data    = [l_data; r_data];
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.simu.Processing.gregressKernel( NLags );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            fprintf('Now saving original data' );
            save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = lagsin_multiple_patched_4_avg1_5lags_labels_and_rdm_100pts_lag()
            
            % Timepoints
            tpts        = [3000:3100];
            
            % Lags
            NLags       = 5;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1';
            
            % Number of cores to run on
            NCores      = 20;
            
            % Legion job name
            jobName     = 'lagsin_multiple_patched_4_avg1_5lags_labels_and_rdm_100pts_lag';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Extract label data
            [l_data, l_vertices] = pitt.exp.simu.Processing.extractLabelDataAllAll( sim );
            
            % Extract random vertices
            [r_data, r_vertices] = pitt.exp.simu.Processing.getRandomVertices( sim, size( l_data,1 ) );
            
            % Construct the job data
            job_data    = [l_data; r_data];
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.simu.Processing.gregressKernel( NLags );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            fprintf('Now saving original data' );
            save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = lagsin_multiple_patched_4_avg1_5lags_labels_and_rdm_100pts_lag2()
            
            % Timepoints
            tpts        = [5000:5100];
            
            % Lags
            NLags       = 5;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1';
            
            % Number of cores to run on
            NCores      = 20;
            
            % Legion job name
            jobName     = 'lagsin_multiple_patched_4_avg1_5lags_labels_and_rdm_100pts_lag2';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Extract label data
            [l_data, l_vertices] = pitt.exp.simu.Processing.extractLabelDataAllAll( sim );
            
            % Extract random vertices
            [r_data, r_vertices] = pitt.exp.simu.Processing.getRandomVertices( sim, size( l_data,1 ) );
            
            % Construct the job data
            job_data    = [l_data; r_data];
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.simu.Processing.gregressKernel( NLags );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            fprintf('Now saving original data' );
            save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = lagsin_multiple_patched_4_avg1_5lags_all_all_100pts()
            
            % Timepoints
            tpts        = [1000:5100];
            
            % Lags
            NLags       = 5;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1';
            
            % Number of cores to run on
            NCores      = 60;
            
            % Legion job name
            jobName     = 'lagsin_multiple_patched_4_avg1_5lags_all_all_100pts_lag2';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Extract label data
            [l_data, l_vertices] = pitt.exp.simu.Processing.extractLabelDataAllAll( sim );
            
            % Extract random vertices
            [r_data, r_vertices] = pitt.exp.simu.Processing.getRandomVertices( sim, size( l_data,1 ) );
            
            % Construct the job data
            job_data    = sim.lhdata.data;
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.simu.Processing.gregressKernel( NLags );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
        function master = lagsin_multiple_patched_4_avg1000_5lags_all_all_100pts()
            
            % Timepoints
            tpts        = [1000:5100];
            
            % Lags
            NLags       = 5;
            
            % Simulation location
            sim_data    = 'lagsin_multiple_patched_4_avg1000';
            
            % Number of cores to run on
            NCores      = 60;
            
            % Legion job name
            jobName     = 'lagsin_multiple_patched_4_avg1000_5lags_all_all_100pts';
            
            
            % generate our sim object with the simulated data
            sim         = pitt.exp.simu.Util.loadSimData( sim_data );
            
            % Extract label data
            [l_data, l_vertices] = pitt.exp.simu.Processing.extractLabelDataAllAll( sim );
            
            % Extract random vertices
            [r_data, r_vertices] = pitt.exp.simu.Processing.getRandomVertices( sim, size( l_data,1 ) );
            
            % Construct the job data
            job_data    = sim.lhdata.data;
            job_data    = job_data( :, tpts );
            
            
            kernel      = pitt.exp.simu.Processing.gregressKernel( NLags );
            master      = legion.DistrMaster( kernel, job_data, NCores, jobName );
            
            % Save how we got this data
            %fprintf('Now saving original data' );
            %save( sprintf('/synapse/logs/schmidtb/%s/extracted_data.mat', jobName), 'l_data', 'l_vertices', 'r_data', 'r_vertices' ); 
            
            % Execute the job
            master.run();
            
        end
        
    end
    
end