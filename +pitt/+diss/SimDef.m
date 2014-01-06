classdef SimDef < handle
   
    properties (Constant)
       
        path_base = '/synapse/home/schmidtb/simu/base/';
        
    end
    
    properties
        
        numaverages = 4;
        sfreq       = 150;
        start_time  = 0;
        end_time    = 600000;
        
        
        subject_name    = 'subj';
        simname         = 'subj-sim.fif';
        
        fwdfile = 'fwd.fif'; % Forward simulation file
        invfile = 'inv.fif'; % Inverse simulation file
        covfile = 'cov.fif'; % Covariance file
        evefile = 'fake_events.eve'; % Artificial events
        avefile = 'fake_ave.ave'; % Artificial average file
        
        working_dir = './';
        
        timefile = 'timecourse.txt';
        
        labelfile = 'label.label';
        
        
    end
    
    methods
        
        
    end
    
    
    
end