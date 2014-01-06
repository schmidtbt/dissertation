classdef Cluster < handle
    
    properties (Constant )
        % This path is on the remote system to matlab (the compute-0-7 part
        % needs to be replaced based on the host (done in program)
        matlabExecution = 'nohup /synapse/pkg/R2011a/compute-0-7/bin/matlab -nosplash -nodesktop -nojvm -r ';
        % Log files will be dumped into this directory
        clusterLogs = '~/matlab/cluster_logs/';
    end
    
    properties
        logLocation = '~/matlab/cluster_logs/';
        hostname = 'compute-0-7';
        jobName = 'default';
    end
    
    
    methods
        function obj = Cluster( hostname, jobName )
            if( nargin == 0 )
                obj.hostname = pitt.Cluster.roundRobin;
                obj.setJobName( 'default' );
            else
                obj.hostname = hostname;
                obj.setJobName( jobName );
            end
        end
        
        function setJobName( obj, jobName )
            obj.jobName = jobName;
            mkdir( [obj.logLocation, obj.jobName ] );
            obj.logLocation = [pitt.Cluster.clusterLogs, obj.jobName, '/' ];
        end
        
        function setHostname( obj, host )
            obj.hostname = host;
        end
        
        function remoteExecute( obj, execString )
            
            executable = strrep( pitt.Cluster.matlabExecution, 'compute-0-7', obj.hostname );
            
            % Add a date-time to the log names
            timenow = datestr(now, 'mmddyyyy-HH-MM-SS-FFF');
            
            % Setup the remote and local log filenames
            localOutput = ['loc_',obj.jobName, '-', timenow];
            remoteOutput = [obj.hostname, '--', obj.jobName,  '-', timenow];
            
            % Escape the " character for shell execution
            executionString = strrep( execString, '"', '\"' );
            
            % Append a quit command to the matlab command
            executionString = ['disp([''Started At:'', datestr(now, ''mmddyyyy-HH-MM-SS-FFF'')]);',executionString, ';disp([''Finished At:'', datestr(now, ''mmddyyyy-HH-MM-SS-FFF'')]);quit;'];
            
            % Assemble the ssh command, use -f to launch a background ssh
            % command before execution. 
            string = ['ssh -f ', obj.hostname, ' "', executable, ' \"' executionString, '\" >> ', obj.logLocation, remoteOutput, '.out 2>> ', obj.logLocation, remoteOutput, '.err < /dev/null &" >> ', obj.logLocation, localOutput, '.out 2>> ', obj.logLocation, localOutput, '.err < /dev/null '];
            fprintf('Executing Output:\t %s\n', remoteOutput);
            fprintf('On Host:\t\t %s\n', obj.hostname);
            fprintf('Now executing command: \t %s\n', string );
            
            % Execute the command
            system( string );
            
            fid = fopen( [obj.logLocation, 'job_spec'], 'w' );
            fprintf(fid, '\nOutput:\t %s\n', remoteOutput);
            fprintf(fid, 'On Host:\t\t %s\n', obj.hostname);
            fprintf(fid, 'String Command:\t\t %s\n', execString);
            fprintf(fid, 'Execute At:\t\t %s\n\n\n', timenow);
            %fprintf(fid, 'Now executing command: \t %s\n', string );
            fclose( fid );
            
        end
    end
    
    methods (Static)
        
        function newLogDir( newLogDir )
            mkdir( newLogDir );
        end
        
        function kill( hostname )
            cmd = ['ssh ', hostname, ' "killall MATLAB" &'];
            disp('Executing:');
            disp( cmd );
            system( cmd ); 
        end
        
        %{
        hostname = which system to enter (ssh keyless)
            If empty, [], will use roundrobin choice for which slave
        executionString = matlab commands to execute on remote machine (a
        quit is appended to this string)
        outputLocation = (Optional) name of job to save output and error
        logs for remote (rem_*) and local (loc_*)
        %}
        function runOnHost( hostname, executionString, outputLocation )
            
            if( isempty(hostname) )
                hostname = pitt.Cluster.roundRobin;
            end
            
            if( nargin ==2 )
                outputLocation = 'foo';
            end
            
            % Replace the hostdependent path of matlab
            executable = strrep( pitt.Cluster.matlabExecution, 'compute-0-7', hostname );
            
            % Add a date-time to the log names
            timenow = datestr(now, 'mmddyyyy-HH-MM-SS-FFF');
            
            % Setup the remote and local log filenames
            localOutput = ['loc_',outputLocation, '-', timenow];
            remoteOutput = ['rem_',outputLocation, '-', hostname, '-', timenow];
            
            % Escape the " character for shell execution
            executionString = strrep( executionString, '"', '\"' );
            
            % Append a quit command to the matlab command
            executionString = [executionString, ';quit;'];
            
            % Assemble the ssh command, use -f to launch a background ssh
            % command before execution. 
            string = ['ssh -f ', hostname, ' "', executable, ' \"' executionString, '\" >> ', pitt.Cluster.clusterLogs, remoteOutput, '.out 2>> ', pitt.Cluster.clusterLogs, remoteOutput, '.err < /dev/null &" >> ', pitt.Cluster.clusterLogs, localOutput, '.out 2>> ', pitt.Cluster.clusterLogs, localOutput, '.err < /dev/null '];
            fprintf('Executing Output:\t %s\n', remoteOutput);
            fprintf('On Host:\t\t %s\n', hostname);
            fprintf('Now executing command: \t %s\n', string );
            
            % Execute the command
            system( string );
            
        end
        
        
        function hostString = hostSwitch( idx )
            hostString = ['compute-0-',num2str(idx)];
        end
        
        function hostString = roundRobin
            hostString = pitt.Cluster.hostSwitch( randi(7,1) );
        end
    end
    
    
end