classdef JobRunner < handle
    
    properties
        useParallel = false;
    end
    
    methods
        function value = runJob(obj, value)
            if( length(value) == 1 )
                %Singular Job to run
                obj.internalRun( value )
            else
                % Multiple Job Files
                if( obj.useParallel )
                    disp('RUNNING IN PARALLEL');
                    matlabpool;
                    parfor i  = 1:length( value )
                        value(i) = obj.internalRun( value(i) );
                    end
                    matlabpool close;
                else
                    for i  = 1:length( value )
                        obj.internalRun( value(i) )
                    end
                end
            end
        end
        
        
        function parallel(obj, bool)
            obj.useParallel = bool;
        end
        
        
        function complete( obj, jobName )
            disp('----------------');
            disp(['Completed: ', jobName ]);
            disp('----------------');
        end
        
        function fail( obj, jobName )
            disp('----------------');
            warning(['FAILED: ', jobName ]);
            disp('----------------');
        end
    end
    
    methods( Access = private )
        function aJob = internalRun( obj, aJob )
            if( aJob.execute )
                obj.complete( aJob.jobName );
            else
                obj.fail( aJob.jobName );
            end
        end
    end
end