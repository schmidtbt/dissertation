classdef Scripts < handle
    
    methods (Static)
        
        
        % --------------------------------
        % Pre processing
        % --------------------------------
        
        function scripts()
            %error('This function is not meant to be run externally. Copy/Pase code for proper execution');
            % --------------------------------
            % Colin
            % --------------------------------
            optColin = pitt.cerebro.Reg.genUserData( 'Colin' );
            
            optColin = pitt.cerebro.DispBrain.outputSurface( optColin, 'pial' );
            
            
            % ----------------------------------------------------------------
            %                       RAW DATA SPACE LOADING
            % ----------------------------------------------------------------
            
            % --------------------------------
            % Opt001
            % --------------------------------
            opt001 = pitt.cerebro.Reg.genUserData( 'Opt001' );
            opt001 = pitt.cerebro.Reg.appendObservedData( opt001, 'rest2', '~/data/r21dataNM/Opt001_resting2_ds4_filt_1_50-ave-Opt001-oct-6-src-fwd.fif' );
            %opt001 = pitt.cerebro.Reg.appendObservedData( opt001, 'rest2', '~/data/r21dataNM/Opt001_resting1_projon_ds4_filt_0_50_tSSS-ave-Opt001-oct-6-src-meg-inv.fif' );
            %opt001 = pitt.cerebro.Reg.appendObservedData( opt001, 'rest2', '~/fwd.fif' );
            opt001 = pitt.cerebro.Reg.appendObservedData( opt001, 'empty', '~/data/r21dataNM/Opt001_resting1_projon_ds4_filt_0_50-ave-Opt001-oct-6-src-fwd.fif' );
            opt001 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt001, 'empty' );
            
            diffL = mean(opt001.empty.surface.L.v) - mean(optColin.pial.surface.L.v);
            diffR = mean(opt001.empty.surface.R.v) - mean(optColin.pial.surface.R.v);
            
            vL =  opt001.rest2.surface.L.v - repmat(diffL, size( opt001.empty.surface.L.v,1 ),1 );
            vR =  opt001.rest2.surface.R.v - repmat(diffR, size( opt001.empty.surface.R.v,1 ),1 );
            
            kL = dsearchn( vL, optColin.pial.surface.L.v );
            kR = dsearchn( vR, optColin.pial.surface.R.v );
            
            opt001.empty.reg.colin.L = kL;
            opt001.empty.reg.colin.R = kR;
            
            
            % --------------------------------
            % Opt052
            % --------------------------------
            opt052 = pitt.cerebro.Reg.genUserData( 'Opt052' );
            opt052 = pitt.cerebro.Reg.appendObservedData( opt052, 'rest2', '~/data/r21dataNM/Opt052_resting2_projon_ds4_filt_0_50_tSSS-ave-Opt052-oct-6-src-fwd.fif' );
            opt052 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt052, 'rest2' );
            opt052 = pitt.cerebro.Reg.regiserLoResToHiResInNewSpace( opt052, optColin, 'rest2', 'colin' );
            
            
            % --------------------------------
            % Opt060
            % --------------------------------
            opt060 = pitt.cerebro.Reg.genUserData( 'Opt060' );
            opt060 = pitt.cerebro.Reg.appendObservedData( opt060, 'rest2', '~/data/r21dataNM/Opt060_resting2_projon_ds4_filt_0_50_tSSS-ave-Opt060-oct-6-src-fwd.fif' );
            opt060 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt060, 'rest2' );
            opt060 = pitt.cerebro.Reg.regiserLoResToHiResInNewSpace( opt060, optColin, 'rest2', 'colin' );
            
            
            % --------------------------------
            % Opt062
            % --------------------------------
            opt062 = pitt.cerebro.Reg.genUserData( 'Opt062' );
            opt062 = pitt.cerebro.Reg.appendObservedData( opt062, 'rest2', '~/data/r21dataNM/Opt062_resting2_projon_ds4_filt_0_50_tSSS-ave-Opt062-oct-6-src-fwd.fif' );
            opt062 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt062, 'rest2' );
            opt062 = pitt.cerebro.Reg.regiserLoResToHiResInNewSpace( opt062, optColin, 'rest2', 'colin' );
            
            % --------------------------------
            % Opt065
            % --------------------------------
            opt065 = pitt.cerebro.Reg.genUserData( 'Opt065' );
            opt065 = pitt.cerebro.Reg.appendObservedData( opt065, 'rest2', '~/data/r21dataNM/Opt065_resting2_projon_ds4_filt_0_50_tSSS-ave-Opt065-oct-6-src-fwd.fif' );
            opt065 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt065, 'rest2' );
            opt065 = pitt.cerebro.Reg.regiserLoResToHiResInNewSpace( opt065, optColin, 'rest2', 'colin' );
            
            
            % --------------------------------
            % Opt067
            % --------------------------------
            opt067 = pitt.cerebro.Reg.genUserData( 'Opt067' );
            opt067 = pitt.cerebro.Reg.appendObservedData( opt067, 'rest2', '~/data/r21dataNM/Opt067_resting2_projon_ds4_filt_0_50_tSSS-ave-Opt067-oct-6-src-fwd.fif' );
            opt067 = pitt.cerebro.Reg.smoothLoResToHiResInSphere( opt067, 'rest2' );
            opt067 = pitt.cerebro.Reg.regiserLoResToHiResInNewSpace( opt067, optColin, 'rest2', 'colin' );
            
            
            opts = struct;
            opts.optColin = optColin;
            opts.opt001 = opt001;
            opts.opt052 = opt052;
            opts.opt060 = opt060;
            opts.opt062 = opt062;
            opts.opt065 = opt065;
            opts.opt067 = opt067;
            
        end
        
        % For preprocessing many frequencies into averaged bands
        function averageFreqBandsInDir( size )
            
            theta = [5:8];
            alpha = [9:13];
            betaL = [14:20];
            betaH = [21:30];
            
            total = zeros( size, size, length(theta));
            fprintf( 'theta\n' );
            for i = theta
                fprintf( 'loading: %s\n', sprintf('freq_%03i',i) );
                total(:,:,i) = pitt.cerebro.DataLoader.loadOptRest( sprintf('freq_%03i',i) );
            end
            output = squeeze( mean( total,3 ) );
            fprintf( 'Saving\n' );
            save(['freq_', 'theta'], 'output' );
            
            
            total = zeros( size, size, length(alpha));
            fprintf( 'alpha\n' );
            for i = alpha
                fprintf( 'loading: %s\n', sprintf('freq_%03i',i) );
                total(:,:,i) = pitt.cerebro.DataLoader.loadOptRest( sprintf('freq_%03i',i) );
            end
            output = squeeze( mean( total,3 ) );
            fprintf( 'Saving\n' );
            save(['freq_', 'alpha'], 'output' );
            
            
            total = zeros( size, size, length(betaL));
            fprintf( 'betaL\n' );
            for i = betaL
                fprintf( 'loading: %s\n', sprintf('freq_%03i',i) );
                total(:,:,i) = pitt.cerebro.DataLoader.loadOptRest( sprintf('freq_%03i',i) );
            end
            output = squeeze( mean( total,3 ) );
            fprintf( 'Saving\n' );
            save(['freq_', 'betaL'], 'output' );
                
            
            total = zeros( size, size, length(betaH));
            fprintf( 'betaH\n' );
            for i = betaH
                fprintf( 'loading: %s\n', sprintf('freq_%03i',i) );
                total(:,:,i) = pitt.cerebro.DataLoader.loadOptRest( sprintf('freq_%03i',i) );
            end
            output = squeeze( mean( total,3 ) );
            fprintf( 'Saving\n' );
            save(['freq_', 'betaH'], 'output' );
            
        end
        
        
        
        % --------------------------------
        %  Primary Processing
        % --------------------------------
        
        function x = processBand( band, opts )
            
            x       = struct;
            x.band  = band;
            
            % --------------------------------
            % Freq Data
            % --------------------------------
            fprintf( 'Loading \n' );
            x.freqopt001er = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt001/empty_done/freq_',band] );
            x.freqopt052r2 = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt052/rest2_done/freq_',band] );
            x.freqopt060r2 = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt060/rest2_done/freq_',band] );
            x.freqopt062r2 = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt062/rest2_done/freq_',band] );
            x.freqopt065r2 = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt065/rest2_done/freq_',band] );
            x.freqopt067r2 = pitt.cerebro.DataLoader.loadOptRest( ['~/data/all/opt067/rest2_done/freq_',band] );
            
            
            % --------------------------------
            % EVC Calculation in Self-Space
            % --------------------------------
            fprintf( 'EVC' );
            x.evcopt052r2 = pitt.cerebro.EVC.eigenVectorCentrality( x.freqopt052r2 ); fprintf( '...done' );
            x.evcopt060r2 = pitt.cerebro.EVC.eigenVectorCentrality( x.freqopt060r2 ); fprintf( '...done' );
            x.evcopt062r2 = pitt.cerebro.EVC.eigenVectorCentrality( x.freqopt062r2 ); fprintf( '...done' );
            x.evcopt065r2 = pitt.cerebro.EVC.eigenVectorCentrality( x.freqopt065r2 ); fprintf( '...done' );
            x.evcopt067r2 = pitt.cerebro.EVC.eigenVectorCentrality( x.freqopt067r2 ); fprintf( '...done' );
            fprintf('\n');
            
            % --------------------------------
            % EVC In Colin Space
            % --------------------------------
            fprintf( 'EVC Colin \n' );
            x.evc52 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opts.opt052.rest2, x.evcopt052r2 );
            x.evc60 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opts.opt060.rest2, x.evcopt060r2 );
            x.evc62 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opts.opt062.rest2, x.evcopt062r2 );
            x.evc65 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opts.opt065.rest2, x.evcopt065r2 );
            x.evc67 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opts.opt067.rest2, x.evcopt067r2 );
            
            
            % --------------------------------
            % Remove the Poor Inverse model reconstruction
            % --------------------------------
            fprintf( 'Ave EVC \n' );
            mevc = mean([x.evc52;x.evc60;x.evc62;x.evc65;x.evc67]);
            mevc( find(mevc<.011) ) = 0;
            
            x.mevc = mevc;
            
            fprintf( 'Preproc EVC' );
            [x.freqopt052r2m, x.evcopt052r2m, x.evc52m] = pitt.cerebro.Scripts.preprocessWithMeanEVC( opts.opt052, mevc, x.freqopt052r2 ); fprintf( '...done' );
            [x.freqopt060r2m, x.evcopt060r2m, x.evc60m] = pitt.cerebro.Scripts.preprocessWithMeanEVC( opts.opt060, mevc, x.freqopt060r2 ); fprintf( '...done' );
            [x.freqopt062r2m, x.evcopt062r2m, x.evc62m] = pitt.cerebro.Scripts.preprocessWithMeanEVC( opts.opt062, mevc, x.freqopt062r2 ); fprintf( '...done' );
            [x.freqopt065r2m, x.evcopt065r2m, x.evc65m] = pitt.cerebro.Scripts.preprocessWithMeanEVC( opts.opt065, mevc, x.freqopt065r2 ); fprintf( '...done' );
            [x.freqopt067r2m, x.evcopt067r2m, x.evc67m] = pitt.cerebro.Scripts.preprocessWithMeanEVC( opts.opt067, mevc, x.freqopt067r2 ); fprintf( '...done' );
            fprintf('\n');
            
            x.mfixevc = mean([x.evc52m;x.evc60m;x.evc62m;x.evc65m;x.evc67m]);
            
            % Create clusters
            fprintf( 'Calc Idxs \n' );
            x.idx = pitt.cerebro.Cluster.kClusterOnSurface( opts.optColin, 'sphere', 40, x.mfixevc.*100'.^3 );
            
            % Create subject specific clustering on self-space
            fprintf( 'LoRes Space \n' );
            x.opt1idx  = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt001, 'empty', x.idx );
            x.opt52idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt052, 'rest2', x.idx );
            x.opt60idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt060, 'rest2', x.idx );
            x.opt62idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt062, 'rest2', x.idx );
            x.opt65idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt065, 'rest2', x.idx );
            x.opt67idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt067, 'rest2', x.idx );
            
            
            % Obtain surface data corresponding to cluster-cluster freq data
            fprintf( 'Orig Data \n' );
            x.opt1loCdata  = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt1idx,  x.freqopt001er );
            x.opt52loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt52idx, x.freqopt052r2m );
            x.opt60loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt60idx, x.freqopt060r2m );
            x.opt62loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt62idx, x.freqopt062r2m );
            x.opt65loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt65idx, x.freqopt065r2m );
            x.opt67loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt67idx, x.freqopt067r2m );
            
            % Generate cell data for rest vs empty conditions
            x.ercell = cell(0);
            x.ercell{1} = x.opt1loCdata;
            x.restcell = cell(0);
            
            x.restcell{1} = x.opt60loCdata;
            x.restcell{2} = x.opt62loCdata;
            x.restcell{3} = x.opt65loCdata;
            x.restcell{4} = x.opt67loCdata;
            x.restcell{5} = x.opt52loCdata;
            
            % Run out permutation testing on cluster-cluster significance
            fprintf( 'Stats \n' );
            [x.tstats,x.output,x.tpop] = pitt.cerebro.Stats.permute_calc_stats( x.restcell, x.ercell );
            
            fprintf( 'Done \n' );
        end
        
        function x = reexec( x, opts )
            
            fprintf( 'LoRes Space \n' );
            x.opt1idx  = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt001, 'empty', x.idx );
            x.opt52idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt052, 'rest2', x.idx );
            x.opt60idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt060, 'rest2', x.idx );
            x.opt62idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt062, 'rest2', x.idx );
            x.opt65idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt065, 'rest2', x.idx );
            x.opt67idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opts.opt067, 'rest2', x.idx );
            
            
            % Obtain surface data corresponding to cluster-cluster freq data
            fprintf( 'Orig Data \n' );
            x.opt1loCdata  = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt1idx,  x.freqopt001er );
            x.opt52loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt52idx, x.freqopt052r2m );
            x.opt60loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt60idx, x.freqopt060r2m );
            x.opt62loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt62idx, x.freqopt062r2m );
            x.opt65loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt65idx, x.freqopt065r2m );
            x.opt67loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( x.opt67idx, x.freqopt067r2m );
            
            % Generate cell data for rest vs empty conditions
            x.ercell = cell(0);
            x.ercell{1} = x.opt1loCdata;
            x.restcell = cell(0);
            
            x.restcell{1} = x.opt60loCdata;
            x.restcell{2} = x.opt62loCdata;
            x.restcell{3} = x.opt65loCdata;
            x.restcell{4} = x.opt67loCdata;
            x.restcell{5} = x.opt52loCdata;
            
            % Run out permutation testing on cluster-cluster significance
            fprintf( 'Stats \n' );
            [x.tstats,x.output,x.tpop] = pitt.cerebro.Stats.permute_calc_stats( x.restcell, x.ercell );
            
        end
        
        function [freqm, evcLom, evcHim] = preprocessWithMeanEVC( opt, mevc, freq )
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt, 'rest2', mevc );
            lst = find( vd > 0 );
            freqm = freq;
            freqm(lst,lst) = freqm(lst,lst)./2;
            evcLom = pitt.cerebro.EVC.eigenVectorCentrality( freqm );
            evcHim = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt.rest2, evcLom );
        end
        
        
        
        % --------------------------------
        %  Outputs
        % --------------------------------
        
        function procTstats( dobj, prepend_name, optColin )
            
            st = pitt.cerebro.Stats.hemisphereTstats( pitt.diss.PLVAnalysis.full_map( dobj.tstats ), [.01, .01, .01], 40 );
            
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'inflated', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, mean(pitt.cerebro.Stats.clusterHemi( st, 'LR', 40 )) ) );
            print(1,[prepend_name,'_LR_inflated_cent'], '-djpeg' );
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'sphere', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, mean(pitt.cerebro.Stats.clusterHemi( st, 'LR', 40 )) ) );
            print(2,[prepend_name,'_LR_sphere_cent'], '-djpeg' );
            
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'inflated', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, mean(pitt.cerebro.Stats.clusterHemi( st, 'I', 40 )) ) );
            print(1,[prepend_name,'_I_inflated_cent'], '-djpeg' );
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'sphere', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, mean(pitt.cerebro.Stats.clusterHemi( st, 'I', 40 )) ) );
            print(2,[prepend_name,'_I_sphere_cent'], '-djpeg' );
            
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'inflated', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, max(pitt.cerebro.Stats.clusterHemi( st, 'LR', 40 )) ) );
            print(1,[prepend_name,'_LR_inflated_pos'], '-djpeg' );
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'sphere', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, max(pitt.cerebro.Stats.clusterHemi( st, 'LR', 40 )) ) );
            print(2,[prepend_name,'_LR_sphere_pos'], '-djpeg' );
            
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'inflated', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, max(pitt.cerebro.Stats.clusterHemi( st, 'I', 40 )) ) );
            print(1,[prepend_name,'_I_inflated_pos'], '-djpeg' );
            pitt.cerebro.DispBrain.displayRawDataOntoSurfaceImgs( optColin, 'sphere', pitt.cerebro.DispBrain.upsampleClusterIdxs( dobj.idx, max(pitt.cerebro.Stats.clusterHemi( st, 'I', 40 )) ) );
            print(2,[prepend_name,'_I_sphere_pos'], '-djpeg' );
            
        end
        
        
        % --------------------------------
        %  Deprecated
        % --------------------------------
        
        function processTheta
            
            % --------------------------------
            % Freq Data
            % --------------------------------
            freq6opt052r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt052/rest2/freq_006' );
            freq6opt060r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt060/rest2/freq_006' );
            freq6opt062r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt062/rest2/freq_006' );
            freq6opt065r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt065/rest2/freq_006' );
            freq6opt067r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt067/rest2/freq_006' );
            
            
            % --------------------------------
            % EVC Calculation in Self-Space
            % --------------------------------
            evc6opt052r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt052r2 );
            evc6opt060r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt060r2 );
            evc6opt062r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt062r2 );
            evc6opt065r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt065r2 );
            evc6opt067r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt067r2 );
            
            
            % --------------------------------
            % EVC In Colin Space
            % --------------------------------
            evc52 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt052.rest2, evc6opt052r2 );
            evc60 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt060.rest2, evc6opt060r2 );
            evc62 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt062.rest2, evc6opt062r2 );
            evc65 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt065.rest2, evc6opt065r2 );
            evc67 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt067.rest2, evc6opt067r2 );
            
            
            % --------------------------------
            % Remove the Poor Inverse model reconstruction
            % --------------------------------
            mevc = mean([evc52;evc60;evc62;evc65;evc67]);
            mevc( find(mevc<.011) ) = 0;
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt052, 'rest2', mevc );
            lst = find( vd > 0 );
            freq6opt052r2m = freq6opt052r2;
            freq6opt052r2m(lst,lst) = freq6opt052r2m(lst,lst)./2;
            evc6opt052r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt052r2m );
            evc52m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt052.rest2, evc6opt052r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt060, 'rest2', mevc );
            lst = find( vd > 0 );
            freq6opt060r2m = freq6opt060r2;
            freq6opt060r2m(lst,lst) = freq6opt060r2m(lst,lst)./2;
            evc6opt060r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt060r2m );
            evc60m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt060.rest2, evc6opt060r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt062, 'rest2', mevc );
            lst = find( vd > 0 );
            freq6opt062r2m = freq6opt062r2;
            freq6opt062r2m(lst,lst) = freq6opt062r2m(lst,lst)./2;
            evc6opt062r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt062r2m );
            evc62m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt062.rest2, evc6opt062r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt065, 'rest2', mevc );
            lst = find( vd > 0 );
            freq6opt065r2m = freq6opt065r2;
            freq6opt065r2m(lst,lst) = freq6opt065r2m(lst,lst)./2;
            evc6opt065r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt065r2m );
            evc65m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt065.rest2, evc6opt065r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt067, 'rest2', mevc );
            lst = find( vd > 0 );
            freq6opt067r2m = freq6opt067r2;
            freq6opt067r2m(lst,lst) = freq6opt067r2m(lst,lst)./2;
            evc6opt067r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq6opt067r2m );
            evc67m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt067.rest2, evc6opt067r2m );
            
            % Create clusters
            idx = pitt.cerebro.Cluster.kClusterOnSurface( optColin, 'sphere', 40, (mean([evc52m;evc60m;evc62m;evc65m;evc67m]).*100)'.^3 );
            
            % Create subject specific clustering on self-space
            opt1idx  = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt001, 'rest2', idx );
            opt52idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt052, 'rest2', idx );
            opt60idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt060, 'rest2', idx );
            opt62idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt062, 'rest2', idx );
            opt65idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt065, 'rest2', idx );
            opt67idx = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt067, 'rest2', idx );
            
            
            % Obtain surface data corresponding to cluster-cluster freq data
            opt1loCdata  = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt1idx, freq6opt001er );
            opt52loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt52idx, freq6opt052r2m );
            opt60loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt60idx, freq6opt060r2m );
            opt62loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt62idx, freq6opt062r2m );
            opt65loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt65idx, freq6opt065r2m );
            opt67loCdata = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt67idx, freq6opt067r2m );
            
            
            
            % Generate cell data for rest vs empty conditions
            ercell = cell(0);
            ercell{1} = opt1loCdata;
            restcell = cell(0);
            
            restcell{1} = opt60loCdata;
            restcell{2} = opt62loCdata;
            restcell{3} = opt65loCdata;
            restcell{4} = opt67loCdata;
            restcell{5} = opt52loCdata;
            
            % Run out permutation testing on cluster-cluster significance
            [tstats,output] = pitt.cerebro.Stats.permute_calc_stats( restcell, ercell );
            
            % Display signif clusters
            figure; imagesc( output )
            
        end
        
        function processAlpha
            
            % --------------------------------
            % Freq Data
            % --------------------------------
            freq10opt001er = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt001/empty/freq_010' );
            freq10opt052r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt052/rest2/freq_010' );
            freq10opt060r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt060/rest2/freq_010' );
            freq10opt062r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt062/rest2/freq_010' );
            freq10opt065r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt065/rest2/freq_010' );
            freq10opt067r2 = pitt.cerebro.DataLoader.loadOptRest( '~/data/opt067/rest2/freq_010' );
            
            
            % --------------------------------
            % EVC Calculation in Self-Space
            % --------------------------------
            evc10opt052r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt052r2 );
            evc10opt060r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt060r2 );
            evc10opt062r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt062r2 );
            evc10opt065r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt065r2 );
            evc10opt067r2 = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt067r2 );
            
            
            % --------------------------------
            % EVC In Colin Space
            % --------------------------------
            evc52 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt052.rest2, evc10opt052r2 );
            evc60 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt060.rest2, evc10opt060r2 );
            evc62 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt062.rest2, evc10opt062r2 );
            evc65 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt065.rest2, evc10opt065r2 );
            evc67 = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt067.rest2, evc10opt067r2 );
            
            
            % --------------------------------
            % Remove the Poor Inverse model reconstruction
            % --------------------------------
            mevc = mean([evc52;evc60;evc62;evc65;evc67]);
            mevc( find(mevc<.011) ) = 0;
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt052, 'rest2', mevc );
            lst = find( vd > 0 );
            freq10opt052r2m = freq10opt052r2;
            freq10opt052r2m(lst,lst) = freq10opt052r2m(lst,lst)./2;
            evc10opt052r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt052r2m );
            evc52m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt052.rest2, evc10opt052r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt060, 'rest2', mevc );
            lst = find( vd > 0 );
            freq10opt060r2m = freq10opt060r2;
            freq10opt060r2m(lst,lst) = freq10opt060r2m(lst,lst)./2;
            evc10opt060r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt060r2m );
            evc60m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt060.rest2, evc10opt060r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt062, 'rest2', mevc );
            lst = find( vd > 0 );
            freq10opt062r2m = freq10opt062r2;
            freq10opt062r2m(lst,lst) = freq10opt062r2m(lst,lst)./2;
            evc10opt062r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt062r2m );
            evc62m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt062.rest2, evc10opt062r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt065, 'rest2', mevc );
            lst = find( vd > 0 );
            freq10opt065r2m = freq10opt065r2;
            freq10opt065r2m(lst,lst) = freq10opt065r2m(lst,lst)./2;
            evc10opt065r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt065r2m );
            evc65m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt065.rest2, evc10opt065r2m );
            
            vd = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLow( opt067, 'rest2', mevc );
            lst = find( vd > 0 );
            freq10opt067r2m = freq10opt067r2;
            freq10opt067r2m(lst,lst) = freq10opt067r2m(lst,lst)./2;
            evc10opt067r2m = pitt.cerebro.EVC.eigenVectorCentrality( freq10opt067r2m );
            evc67m = pitt.cerebro.DispBrain.upsampleToHiResRegisteredDualHemi( opt067.rest2, evc10opt067r2m );
            
            % Create clusters
            idx10 = pitt.cerebro.Cluster.kClusterOnSurface( optColin, 'sphere', 40, (mean([evc52m;evc60m;evc62m;evc65m;evc67m]).*100)'.^3 );
            
            % Create subject specific clustering on self-space
            opt1idx10  = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt001, 'rest2', idx10 );
            opt52idx10 = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt052, 'rest2', idx10 );
            opt60idx10 = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt060, 'rest2', idx10 );
            opt62idx10 = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt062, 'rest2', idx10 );
            opt65idx10 = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt065, 'rest2', idx10 );
            opt67idx10 = pitt.cerebro.DispBrain.downsampleHiResRegisteredToLowByMode( opt067, 'rest2', idx10 );
            
            
            % Obtain surface data corresponding to cluster-cluster freq data
            opt1loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt1idx10, freq10opt001er );
            opt52loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt52idx10, freq10opt052r2m );
            opt60loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt60idx10, freq10opt060r2m );
            opt62loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt62idx10, freq10opt062r2m );
            opt65loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt65idx10, freq10opt065r2m );
            opt67loCdata10 = pitt.cerebro.DispBrain.loResPopSampleFromHiResCluster( opt67idx10, freq10opt067r2m );
            
            % Generate cell data for rest vs empty conditions
            ercell = cell(0);
            ercell{1} = opt1loCdata10;
            restcell = cell(0);
            
            restcell{1} = opt60loCdata10;
            restcell{2} = opt62loCdata10;
            restcell{3} = opt65loCdata10;
            restcell{4} = opt67loCdata10;
            restcell{5} = opt52loCdata10;
            
            % Run out permutation testing on cluster-cluster significance
            [tstats,output] = pitt.cerebro.Stats.permute_calc_stats( restcell, ercell );
            
            % Display signif clusters
            figure; imagesc( output )
            
        end
        
        
        
    end
    
end

