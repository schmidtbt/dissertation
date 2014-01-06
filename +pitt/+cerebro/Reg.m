classdef Reg < handle
    
    methods (Static)
        
        function opt = genUserData( optNameUpper )
            
            opt = struct();
            
            opt.name = optNameUpper;
            opt.namelower = lower( optNameUpper );
            
            [opt.inflated.surface.L.v,opt.inflated.surface.L.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/lh.inflated'] );
            [opt.inflated.surface.R.v,opt.inflated.surface.R.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/rh.inflated'] );
            
            [opt.sphere_reg.surface.L.v,opt.sphere_reg.surface.L.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/lh.sphere.reg'] );
            [opt.sphere_reg.surface.R.v,opt.sphere_reg.surface.R.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/rh.sphere.reg'] );
            
            [opt.sphere.surface.L.v,opt.sphere.surface.L.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/lh.sphere'] );
            [opt.sphere.surface.R.v,opt.sphere.surface.R.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/rh.sphere'] );
            
            [opt.pial.surface.L.v,opt.pial.surface.L.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/lh.pial'] );
            [opt.pial.surface.R.v,opt.pial.surface.R.f] = mne_read_surface( ['~/data/sMRI/',optNameUpper,'/surf/rh.pial'] );
            
        end
        
        % We allow only a single trial data per subject (we cannot store multiple observed data fields
        function opt = appendObservedData( opt, named_field, fwd_fif_file )
            
            % We assume that the index numbers correspond to those from the FreeSurfer registration in genuserdata()
            % Otherwise we're screwed regardless and it's all pointless.
            
            s               = mne_read_source_spaces( fwd_fif_file );
            
            % We extract the used vertices from the fwd files data b/c mne is too lazy to do it for us
            usedL           = s(1).rr( find( s(1).inuse ), : );
            usedR           = s(2).rr( find( s(2).inuse ), : );
            
            % Setup the new vertex system
            opt.(named_field).surface.L.v = usedL;
            opt.(named_field).surface.R.v = usedR;
            
            % Record original vertex numbers used:
            opt.(named_field).orig.L.v    = s(1).vertno;
            opt.(named_field).orig.R.v    = s(2).vertno;
            
            % Now we have to find the new corresponding face indices to go with our smaller mesh (thanks again MNE)
            Lu   = find( s(1).inuse );
            tris = s(1).use_tris;
            trisConv = tris;
            for i = 1:length(Lu)
                lst = find( tris == Lu(i));
                if( length(lst) == 0 ); error('here'); end;
                trisConv(lst) = i;
            end
            
            % Hack for bad downsampled v,f pairs - remove MNE vertex labels from faces which did not appear in the new
            % list
            lst = find( trisConv > length(Lu) );
            trisConv(lst) = 1;
            if( length(lst) > 0 ); warning('Implemented'); end;
            
            opt.(named_field).surface.L.f = trisConv;
            
            
            
            Lu   = find( s(2).inuse );
            tris = s(2).use_tris;
            trisConv = tris;
            for i = 1:length(Lu)
                lst = find( tris == Lu(i));
                trisConv(lst) = i;
            end
            
            % Hack for bad downsampled v,f pairs
            lst = find( trisConv > length(Lu) );
            trisConv(lst) = 1;
            
            
            opt.(named_field).surface.R.f = trisConv;
            
        end
        
        % Smooth loRes sampling onto hiRes
        function opt = smoothLoResToHiResInSphere( opt, loRes_named_field )
            
            vL = opt.(loRes_named_field).surface.L.v;
            vR = opt.(loRes_named_field).surface.R.v;
            
            vHL = opt.sphere.surface.L.v;
            vHR = opt.sphere.surface.R.v;
            
            kL = dsearchn( vL, vHL );
            kR = dsearchn( vR, vHR );
            
            opt.(loRes_named_field).hiRes.L = kL;
            opt.(loRes_named_field).hiRes.R = kR;
            
        end
        
        function optLo = regiserLoResToHiResInNewSpace( optLo, optCommon, loRes_named_field, reg_named_field )
            
            % Low dim vertex space
            %vL      = optLo.(loRes_named_field).surface.L.v;
            %vR      = optLo.(loRes_named_field).surface.R.v;
            
            vIdxL   =  optLo.(loRes_named_field).orig.L.v(1:3:end,:);
            vIdxR   =  optLo.(loRes_named_field).orig.R.v(1:3:end,:);
            
            vLoRegL = optLo.sphere_reg.surface.L.v;
            vLoRegR = optLo.sphere_reg.surface.R.v;
            
            % Hack
            lstL = find( vIdxL > size( vLoRegL,1 ) );
            lstR = find( vIdxR > size( vLoRegR,1 ) );
            vIdxL(lstL) = randi( size( vLoRegL,1 ), length(lstL),1 );
            vIdxR(lstR) = randi( size( vLoRegR,1 ), length(lstR),1 );
            
            
            vRegL   = vLoRegL( vIdxL,: );
            vRegR   = vLoRegR( vIdxR,: );
            
            % Common registered space
            vHLcm   = optCommon.sphere_reg.surface.L.v;
            vHRcm   = optCommon.sphere_reg.surface.R.v;
            
            % Search for common space
            kL      = dsearchn( vRegL, vHLcm );
            kR      = dsearchn( vRegR, vHRcm );
            
            optLo.(loRes_named_field).reg.(reg_named_field).L = kL;
            optLo.(loRes_named_field).reg.(reg_named_field).R = kR;
            
        end
        
        
    end
    
end

