classdef StationaryScales < handle
    
    properties
        data;
    end
    
    methods
        
        % Data = numLocations x numTimepoints
        function obj = StationaryScales( data )
            obj.data = data;
        end
        
    end
    
    methods (Access=private)
       
    end
    
    methods (Static)
        
        % 0 Output indicates stationarity
        function output = runADF( data, numScales )
            s = pitt.jobs.Scale( data, pitt.exp.StationaryScales.genADFPipe() );
            warning off all
            output = s.runByWindow( numScales, 0 );
            warning on all
        end
        
        % 1 Output indicates stationarity
        function output = runKPSS( data, numScales )
            s = pitt.jobs.Scale( data, pitt.exp.StationaryScales.genKPSSPipe() );
            warning off all
            output = s.runByWindow( numScales, 0 );
            warning on all
        end
        
        function pipe = genADFPipe()
            
            pipe = pitt.jobs.Pipe; 
            pipe.add( @cca_detrend, 'X');
            pipe.add( @cca_rm_temporalmean, 'X');
            pipe.add( @cca_check_cov_stat, 'X', 10 );
            
        end
        
        function pipe = genKPSSPipe()
            
            pipe = pitt.jobs.Pipe; 
            pipe.add( @cca_detrend, 'X');
            pipe.add( @cca_rm_temporalmean, 'X');
            pipe.add( @cca_kpss, 'X' );
            
        end
        
        
        % Save data then re-load it later. Executed locally
        function distrRunADF( data, jobName, numscales )
            grid = pitt.grid.Grid( jobName );
            save( ['/synapse/logs/schmidtb/',jobName,'/data.mat'], 'data' );
            grid.setExecutable( [ 'pitt.exp.StationaryScales.remoteRunADF( ', num2str(numscales),' );' ] );
            grid.executePBS();
        end
        
        % Executed on remote host to begin processing there (via load
        % command of previously saved data
        function remoteRunADF( numscales )
            
            disp('Adding dependencies');
            pitt.Depend.CCAadd;
            which cca_adf
            disp('loading');
            load('data.mat');
            disp('executing');
            output = pitt.exp.StationaryScales.runADF( data, numscales );
            disp('saving');
            save( 'output.mat', 'output', '-v7.3'  );
            
        end
        
        
        % Save data then re-load it later. Executed locally
        function distrRunKPSS( data, jobName, numscales )
            grid = pitt.grid.Grid( jobName );
            save( ['/synapse/logs/schmidtb/',jobName,'/data.mat'], 'data' );
            grid.setExecutable( [ 'pitt.exp.StationaryScales.remoteRunKPSS( ', num2str(numscales),' );' ] );
            grid.executePBS();
        end
        
        % Executed on remote host to begin processing there (via load
        % command of previously saved data
        function remoteRunKPSS( numscales )
            
            disp('Adding dependencies');
            pitt.Depend.CCAadd;
            disp('loading');
            load('data.mat');
            disp('executing');
            output = pitt.exp.StationaryScales.runKPSS( data, numscales );
            disp('saving');
            save( 'output.mat', 'output', '-v7.3'  );
            
        end
        
        
    end
    
    
end