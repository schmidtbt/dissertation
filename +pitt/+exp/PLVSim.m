classdef PLVSim < handle
   
    methods ( Static )
        
        function sim_plv = simulate_plv_random
           
            N = 10;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            
        end
        
        %{
        
        Basic plv locking
        A is sum of B and C at different frequences
        
        sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic
           
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 5;
            F2 = 25;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            A = sind( F1*pi*t ) + sind( F2*pi*t );
            B = sind( F1*pi*t );
            C = sind( F2*pi*t );
            
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        function sim_plv = simulate_plv_abc_basic_2
           
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            A = sind( F1*pi*t ) + sind( F2*pi*t );
            B = sind( F1*pi*t );
            C = sind( F2*pi*t );
            
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_3;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_3
           
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            A = sind( F1*pi*t ) + sind( F2*pi*t );
            B = sind( F1*pi*t );
            C = sind( F2*pi*t );
            
            
            %sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_rand_with_sin;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_rand_with_sin
           
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            %A = sind( F1*pi*t ) + sind( F2*pi*t );
            B = sind( F1*pi*t );
            %C = sind( F2*pi*t );
            
            
            %sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            %sim_plv( Cidx, : ) = C;
            
            
        end
        
        % Fixed offsets -- complete phase locked between all nodes
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_4;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_4
           
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            A = sind( F1*pi*t + 45 ) + sind( F1*pi*t + 90 );
            B = sind( F1*pi*t + 45 );
            C = sind( F1*pi*t + 90 );
            
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_5;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_5
            
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            Factor1 = 10;
            Factor2 = 10;
            
            Factor1 = 200;
            Factor2 = 200;
            
            F1phase = Factor1.*randn( 1,S );
            F2phase = Factor2.*randn( 1,S );
            
            %F1phase = zeros( 1,S );
            %F2phase = zeros( 1,S );
            
            %Off1 = 45;
            %Off2 = 900;
            
            Off1 = 0;
            Off2 = 0;
            
            A = sind( F1*pi*t + F1phase + Off1 ) + sind( F1*pi*t + F2phase + Off2 );
            B = sind( F1*pi*t + F1phase + Off1 );
            C = sind( F1*pi*t + F2phase + Off2 );
            
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_6;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_6
            
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            Factor1 = 10;
            Factor2 = 10;
            
            %Factor1 = 200;
            %Factor2 = 200;
            
            F1phase = Factor1.*randn( 1,S );
            F2phase = Factor2.*randn( 1,S );
            
            %F1phase = zeros( 1,S );
            %F2phase = zeros( 1,S );
            
            %Off1 = 45;
            %Off2 = 180;
            
            Off1 = 0;
            Off2 = 0;
            
            SNRA = 3;
            SNRB = 3;
            SNRC = 3;
            
            randA = SNRA .* randn(1, S );
            randB = SNRB .* randn(1, S );
            randC = SNRC .* randn(1, S );
            
            
            A = sind( F1*pi*t + F1phase + Off1 ) + sind( F1*pi*t + F2phase + Off2 ) + randA;
            B = sind( F1*pi*t + F1phase + Off1 ) + randB;
            C = sind( F1*pi*t + F2phase + Off2 ) + randC;
            
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        
        %{
        Example of butterworth filtering for PLV calculation
        
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_7;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_7
            
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            
            [B,A] = butter( 5, [5/45, 10/45] );
            
            y = filter( B,A,randn( 1,S ) );
            
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            Factor1 = 10;
            Factor2 = 10;
            
            %Factor1 = 200;
            %Factor2 = 200;
            
            F1phase = Factor1.*randn( 1,S );
            F2phase = Factor2.*randn( 1,S );
            
            %F1phase = zeros( 1,S );
            %F2phase = zeros( 1,S );
            
            %Off1 = 45;
            %Off2 = 180;
            
            Off1 = 0;
            Off2 = 0;
            
            SNRA = 1;
            SNRB = 1;
            SNRC = 3;
            
            randA = SNRA .* randn(1, S );
            randB = SNRB .* randn(1, S );
            randC = SNRC .* randn(1, S );
            
            
            %A = sind( F1*pi*t + F1phase + Off1 ) + sind( F1*pi*t + F2phase + Off2 ) + randA;
            %B = sind( F1*pi*t + F1phase + Off1 ) + randB;
            %C = sind( F1*pi*t + F2phase + Off2 ) + randC;
            
            
            C = y + randA;
            B = y + randB;
            A = B + C + randC;
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        %{
            sim_plv = pitt.exp.PLVSim.simulate_plv_abc_basic_8;
            plv = pitt.exp.PLVSim.calculate_plv( sim_plv );
            figure(5); pitt.exp.PLVSim.disp( sim_plv, plv );
        %}
        function sim_plv = simulate_plv_abc_basic_8
            
            N = 4;
            S = 10000;
            
            for i = 1:N;
                
                sim_plv(i,:) = randn( 1, S );
                
            end
            
            F1 = 10;
            F2 = 21;
            
            
            [B,A] = butter( 5, [5/45, 10/45] );
            
            y1 = filter( B,A,randn( 1,S ) );
            y2 = filter( B,A,randn( 1,S ) );
            
            Aidx = 1;
            Bidx = 2;
            Cidx = 3;
            
            t = 1:S;
            
            SNRA = 1.6;
            SNRB = .2;
            SNRC = .2;
            
            randA = SNRA .* randn(1, S );
            randB = SNRB .* randn(1, S );
            randC = SNRC .* randn(1, S );
            
            %C = y1;% + randA;
            C = y1 + randC;
            %B = y2;% + randB;
            B = y2 + randB;
            A = B + C + randA;
            
            sim_plv( Aidx, : ) = A;
            sim_plv( Bidx, : ) = B;
            sim_plv( Cidx, : ) = C;
            
            
        end
        
        function plv = calculate_plv( data, freq )
            
            num_vertex = size( data,1 );
            
            
            Fs = 150;
            Freqs = 5:40;
            
            plv = -1 .* ones( length( Freqs ), num_vertex, num_vertex );
            
            Fs = 150;
            Freqs = 5:40;
            
            kernel = pitt.exp.plv.Kernels.plvkernel(Fs, Freqs);
            
            
            for i = 1:num_vertex
                
                for j = 1:num_vertex
                    
                    if( j>i ); continue; end;
                    
                    kernel.initial( [data(i,:); data(j,:)] );
                    output = kernel.execute;
                    
                    for k =Freqs
                        plv(k-Freqs(1)+1,i,j) = output{k};
                    end
                    
                end
                
            end
            
            lst = find( plv == -1 );            
            plv( lst ) = 0;
            
            disp( 'finished calculating' );
        end
        
        function [AB,AC,BC] = extract_plv( plv )
            
            [f, x, y] = size( plv );
            
            AB = -1 .* ones( 1,f );
            AC = -1 .* ones( 1,f );
            BC = -1 .* ones( 1,f );
            
            for i = 1:f
                
                AB(i) = plv( i,2,1 );
                AC(i) = plv( i,3,1 );
                BC(i) = plv( i,3,2 );
                
            end
            
        end
        
        function disp( sim_plv, plv )
            
            [AB,AC,BC] = pitt.exp.PLVSim.extract_plv( plv );
            
            subplot( 4,3,[1 2 3] );
            plot( sim_plv(1:3,1:100)' );
            legend('A','B','C');
            
            subplot( 4,3,[4 5 6] );
            plot( 5:40, [AB; AC; BC] );
            legend('AB','AC','BC');
            
            subplot( 4,3,7 );
            pwelch( sim_plv( 1,: ) ); title( 'A' );
            
            
            subplot( 4,3,8 );
            pwelch( sim_plv( 2,: ) ); title( 'B' );
            
            
            subplot( 4,3,9 );
            pwelch( sim_plv( 3,: ) ); title( 'C' );
            
            subplot( 4,3,10 );
            mscohere( sim_plv( 1,: ), sim_plv( 2,: ) ); title( 'AB' );
            subplot( 4,3,11 );
            mscohere( sim_plv( 1,: ), sim_plv( 3,: ) ); title( 'AC' );
            subplot( 4,3,12 );
            mscohere( sim_plv( 2,: ), sim_plv( 3,: ) ); title( 'BC' );
            
        end
        
    end
    
    
end