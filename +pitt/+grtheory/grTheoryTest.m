% Test propitt.grtheory.gram for pitt.grtheory.grTheory functions
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

clear all

% select test:
% ntest=1 - pitt.grtheory.grBase
% ntest=2 - pitt.grtheory.grCoBase
% ntest=3 - pitt.grtheory.grCoCycleBasis
% ntest=4 - pitt.grtheory.grColEdge
% ntest=5 - pitt.grtheory.grColVer
% ntest=6 - pitt.grtheory.grComp
% ntest=7 - pitt.grtheory.grCycleBasis
% ntest=8 - pitt.grtheory.grDecOrd
% ntest=9 - pitt.grtheory.grDistances
% ntest=10 - pitt.grtheory.grEccentricity
% ntest=11 - pitt.grtheory.grIsEulerian
% ntest=12 - pitt.grtheory.grIsomorph
% ntest=13 - pitt.grtheory.grMaxComSu
% ntest=14 - pitt.grtheory.grMaxFlows
% ntest=15 - pitt.grtheory.grMaxMatch
% ntest=16 - pitt.grtheory.grMaxStabSet
% ntest=17 - pitt.grtheory.grMinAbsEdgeSet
% ntest=18 - pitt.grtheory.grMinAbsVerSet
% ntest=19 - pitt.grtheory.grMinCutSet
% ntest=20 - pitt.grtheory.grMinEdgeCover
% ntest=21 - pitt.grtheory.grMinSpanTree
% ntest=22 - pitt.grtheory.grMinVerCover
% ntest=23 - pitt.grtheory.grPERT
% ntest=24 - pitt.grtheory.grPlot
% ntest=25 - pitt.grtheory.grShortPath
% ntest=26 - pitt.grtheory.grShortVerPath
% ntest=27 - pitt.grtheory.grTranClos
% ntest=28 - pitt.grtheory.grTravSale

ntest=26;

switch ntest % selected test
  case 1, % pitt.grtheory.grBase test
    disp('The pitt.grtheory.grBase test')
    V=[0 4;4 4;0 0;4 0;8 4;8 0;12 4;12 0;...
      -2 8;-4 4;0 -4;-4 -4;-4 0];
    E=[1 2;3 4;2 4;2 5;4 6;5 6;5 7;6 8;8 7;...
       1 9;9 10;10 1;3 11;11 12;12 13;13 3;3 12;13 11];
    pitt.grtheory.grPlot(V,E,'d','%d','');
    title('\bfThe initial dipitt.grtheory.graph')
    BG=pitt.grtheory.grBase(E);
    disp('The bases of dipitt.grtheory.graph:')
    disp(' N    vertexes')
    for k1=1:size(BG,1),
      fprintf('%2.0f    ',k1)
      fprintf('%d  ',BG(k1,:))
      fprintf('\n')
    end
  case 2, % pitt.grtheory.grCoBase test
    disp('The pitt.grtheory.grCoBase test')
    V=[0 4;4 4;0 0;4 0;8 4;8 0;12 4;12 0;...
      -2 8;-4 4;0 -4;-4 -4;-4 0];
    E=[2 1;3 4;2 4;2 5;4 6;5 6;5 7;6 8;8 7;...
       1 9;9 10;10 1;3 11;11 12;12 13;13 3;3 12;13 11];
    pitt.grtheory.grPlot(V,E,'d','%d','');
    title('\bfThe initial dipitt.grtheory.graph')
    CBG=pitt.grtheory.grCoBase(E);
    disp('The contrabasis of dipitt.grtheory.graph:')
    disp(' N    vertexes')
    for k1=1:size(CBG,1),
      fprintf('%2.0f    ',k1)
      fprintf('%d  ',CBG(k1,:))
      fprintf('\n')
    end
  case 3, % pitt.grtheory.grCoCycleBasis test
    disp('The pitt.grtheory.grCoCycleBasis test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    E=[E [1:size(E,1)]']; % edges with numbers
    pitt.grtheory.grPlot(V,E,'g','','%d'); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    CoCycles=pitt.grtheory.grCoCycleBasis(E); % all independences cut-sets
    for k1=1:size(CoCycles,2),
      pitt.grtheory.grPlot(V,E(find(~CoCycles(:,k1)),:),'g','','%d'); % one cocycle
      title(['\bfCocycle N' num2str(k1)]);
    end
  case 4, % pitt.grtheory.grColEdge test
    disp('The pitt.grtheory.grColEdge test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V,E,'g','','%d'); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    mCol=pitt.grtheory.grColEdge(E); % the color problem for edges
    fprintf('The colors of edges\n N edge    N col\n');
    fprintf('  %2.0f        %2.0f\n',[1:length(mCol);mCol']);
    pitt.grtheory.grPlot(V,[E,mCol],'g','','%d'); % plot the colored pitt.grtheory.graph
    title('\bfThe pitt.grtheory.graph with colored edges')
  case 5, % pitt.grtheory.grColVer test
    disp('The pitt.grtheory.grColVer test')
    t=[0:4]';
    V=[[5*sin(2*pi*t/5) 5*cos(2*pi*t/5)];...
       [4*sin(2*pi*(t-0.5)/5) 4*cos(2*pi*(t-0.5)/5)];...
       [2*sin(2*pi*(t-0.5)/5) 2*cos(2*pi*(t-0.5)/5)];[0 0]];
    E=[1 7;7 2;2 8;8 3;3 9;9 4;4 10;10 5;5 6;6 1;...
      1 10;2 6;3 7;4 8;5 9;1 12;2 13;3 14;4 15;5 11;...
      6 14;7 15;8 11;9 12;10 13;...
      11 12;12 13;13 14;14 15;15 11;1 16;2 16;3 16;4 16;5 16];
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    nCol=pitt.grtheory.grColVer(E); % the color problem for vertexes
    fprintf('The colors of vertexes\n N ver     N col\n');
    fprintf('  %2.0f        %2.0f\n',[1:length(nCol);nCol']);
    pitt.grtheory.grPlot([V,nCol],E,'g','%d',''); % plot the colored pitt.grtheory.graph
    title('\bfThe pitt.grtheory.graph with colored vertexes')
  case 6, % pitt.grtheory.grComp test
    disp('The pitt.grtheory.grComp test')
    V=[0 4;4 4;0 0;4 0;8 4;8 0;12 4;12 0;...
      -2 8;-4 4;0 -4;-4 -4;-4 0];
    E=[1 2;3 4;2 4;2 5;4 6;5 6;5 7;6 8;8 7;...
       1 9;9 10;10 1;3 11;11 12;12 13;13 3;3 12;13 11];
    ncV=pitt.grtheory.grComp(E);
    pitt.grtheory.grPlot([V ncV],E,'g','%d','');
    title(['\bfThis pitt.grtheory.graph have ' num2str(max(ncV)) ' component(s)'])
    E=[2 4;2 5;4 6;5 6;5 7;6 8;8 7;...
       1 9;9 10;10 1;3 11;11 12;12 13;13 3;3 12;13 11];
    ncV=pitt.grtheory.grComp(E);
    pitt.grtheory.grPlot([V ncV],E,'g','%d','');
    title(['\bfThis pitt.grtheory.graph have ' num2str(max(ncV)) ' component(s)'])
    E=[2 4;2 5;4 6;5 6;5 7;6 8;8 7;...
       1 9;9 10;10 1;3 11;11 12;3 12];
    ncV=pitt.grtheory.grComp(E,size(V,1));
    pitt.grtheory.grPlot([V ncV],E,'g','%d','');
    title(['\bfThis pitt.grtheory.graph have ' num2str(max(ncV)) ' component(s)'])
  case 7, % pitt.grtheory.grCycleBasis test
    disp('The pitt.grtheory.grCycleBasis test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    Cycles=pitt.grtheory.grCycleBasis(E); % all independences cycles
    for k1=1:size(Cycles,2),
      pitt.grtheory.grPlot(V,E(find(Cycles(:,k1)),:),'g','%d',''); % one cycle
      title(['\bfCycle N' num2str(k1)]);
    end
  case 8, % pitt.grtheory.grDecOrd test
    disp('The pitt.grtheory.grDecOrd test')
    V=[0 4;1 4;2 4;3 4;4 4;0 3;1 3;2 3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2;3 2;4 3;5 4;6 1;2 7;8 2;3 8;9 4;9 5;10 5;7 6;8 7;...
       8 9;10 9;11 6;7 12;13 8;14 9;15 10;12 11;13 12;13 14;...
       14 13;15 14;16 11;12 17;13 18;20 15;17 16;17 18;18 17;...
       19 18;19 20;21 16;17 22;18 22;22 18;18 23;19 24;20 25;...
       21 22;22 21;23 24;24 23;24 25];
    pitt.grtheory.grPlot(V,E,'d');
    title('\bfThe initial dipitt.grtheory.graph')
    [Dec,Ord]=pitt.grtheory.grDecOrd(E); % solution
    disp('The classes of mutually connected vertexes:')
    disp(' N    vertexes')
    for k1=1:size(Dec,2),
      fprintf('%2.0f    ',k1)
      fprintf('%d  ',find(Dec(:,k1)))
      fprintf('\n')
    end
    fprintf('The partial ordering of the classes:\n  ')
    fprintf('%3.0f',1:size(Ord,2))
    fprintf('\n')
    for k1=1:size(Ord,1), % the matrix of partial ordering
      fprintf('%2.0f ',k1)
      fprintf(' %1.0f ',Ord(k1,:))
      fprintf('\n')
    end
    V1=V;
    for k1=1:size(Dec,2),
      V1(find(Dec(:,k1)),3)=k1; % weight = number of the class
    end
    pitt.grtheory.grPlot(V1,E,'d','%d','');
    title('\bfThe classes of mutually connected vertexes')
  case 9, % pitt.grtheory.grDispances test
    disp('The pitt.grtheory.grDispances test')
    V=[0 4;1 4;2 4;3 4;4 4;0 3;1 3;2 3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2 1;3 2 2;4 3 3;5 4 4;6 1 5;2 7 6;8 2 7;3 8 1;...
       9 4 9;9 5 8;10 5 7;7 6 6;8 7 5;8 9 4;10 9 3;11 6 2;...
       7 12 1;13 8 2;14 9 3;15 10 4;12 11 5;13 12 6;13 14 7;...
       14 13 8;5 5 10;15 14 9;16 11 8;12 17 7;13 18 6;...
       20 15 5;17 16 4;17 18 3;18 17 2;19 18 1;19 20 2;...
       5 5 8; 21 16 3;17 22 4;18 22 5;22 18 6;18 23 7;...
       19 24 8;20 25 9;21 22 8;22 21 7;23 24 6;10 10 8;...
       24 23 5;24 25 4];
    s=1; % one vertex
    t=25; % other vertex
    pitt.grtheory.grPlot(V,E); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    [dSP,sp]=pitt.grtheory.grDistances(E(:,1:2),s,t); % the nonweighted task
    fprintf('The steps number between all vertexes:\n   ')
    fprintf('%4.0f',1:size(dSP,2))
    fprintf('\n')
    for k1=1:size(dSP,1), % the matrix of distances
      fprintf('%3.0f ',k1)
      fprintf(' %2.0f ',dSP(k1,:))
      fprintf('\n')
    end
    pitt.grtheory.grPlot(V(:,1:2),[sp(1:end-1);sp(2:end)]','d','%d','')
    title(['\bfThe shortest way between vertex ' ...
      num2str(s) ' and vertex ' num2str(t)])
    [dSP,sp]=pitt.grtheory.grDistances(E,1,25); % the weighted task
    fprintf('The distances between all vertexes:\n   ')
    fprintf('%4.0f',1:size(dSP,2))
    fprintf('\n')
    for k1=1:size(dSP,1), % the matrix of distances
      fprintf('%3.0f ',k1)
      fprintf(' %2.0f ',dSP(k1,:))
      fprintf('\n')
    end
    pitt.grtheory.grPlot(V(:,1:2),[sp(1:end-1);sp(2:end)]','d','%d','')
    title(['\bfThe way with minimal weight between vertex ' ...
      num2str(s) ' and vertex ' num2str(t)])
  case 10, % pitt.grtheory.grEccentricity test
    disp('The pitt.grtheory.grEccentricity test')
    V=[0 4;1 4;2 4;3 4;4 4;0 3;1 3;2 3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2 1;3 2 2;4 3 3;5 4 4;6 1 5;2 7 6;8 2 7;3 8 1;...
       9 4 9;9 5 8;10 5 7;7 6 6;8 7 5;8 9 4;10 9 3;11 6 2;...
       7 12 1;13 8 2;14 9 3;15 10 4;12 11 5;13 12 6;13 14 7;...
       14 13 8;5 5 10;15 14 9;16 11 8;12 17 7;13 18 6;...
       20 15 5;17 16 4;17 18 3;18 17 2;19 18 1;19 20 2;...
       5 5 8; 21 16 3;17 22 4;18 22 5;22 18 6;18 23 7;...
       19 24 8;20 25 9;21 22 8;22 21 7;23 24 6;10 10 8;...
       24 23 5;24 25 4];
    pitt.grtheory.grPlot(V,E); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    [Ec,Rad,Diam,Cv,Pv]=pitt.grtheory.grEccentricity(E(:,1:2)); % the nonweighted task
    fprintf('The nonweighted eccentricities of all vertexes:\n N ver  Ecc\n')
    fprintf('%4.0f  %4.0f\n',[1:size(V,1);Ec])
    fprintf('The radius of pitt.grtheory.graph Rad=%d\n',Rad)
    fprintf('The diameter of pitt.grtheory.graph Diam=%d\n',Diam)
    fprintf('The center vertexes is:')
    fprintf('  %d',Cv)
    fprintf('\nThe periphery vertexes is:')
    fprintf('  %d',Pv)
    fprintf('\n')
    [Ec,Rad,Diam,Cv,Pv]=pitt.grtheory.grEccentricity(E); % the weighted task
    fprintf('The weighted eccentricities of all vertexes:\n N ver  Ecc\n')
    fprintf('%4.0f  %4.0f\n',[1:size(V,1);Ec])
    fprintf('The radius of pitt.grtheory.graph Rad=%d\n',Rad)
    fprintf('The diameter of pitt.grtheory.graph Diam=%d\n',Diam)
    fprintf('The center vertexes is:')
    fprintf('  %d',Cv)
    fprintf('\nThe periphery vertexes is:')
    fprintf('  %d',Pv)
    fprintf('\n')
  case 11, % pitt.grtheory.grIsEulerian test
    disp('The pitt.grtheory.grIsEulerian test')
    V=[0 4;1 4;2 4;3 4;4 4;0 3;1 3;2 3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2;3 2;4 3;5 4;6 1;2 7;8 2;3 8;9 4;9 5;10 5;7 6;8 7;...
       8 9;10 9;11 6;7 12;13 8;14 9;15 10;12 11;13 12;13 14;...
       14 13;15 14;16 11;12 17;13 18;20 15;17 16;17 18;18 17;...
       19 18;19 20;21 16;17 22;18 22;22 18;18 23;19 24;20 25;...
       21 22;22 21;23 24;24 23;24 25];
    [eu,cEu]=pitt.grtheory.grIsEulerian(E);
    switch eu,
      case 1,
        st='';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      case 0.5,
        st='semi-';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      otherwise,
        st='not ';
        E=[E(:,1:2), [1:size(E,1)]'];
    end
    pitt.grtheory.grPlot(V,E,'g','','%d');
    title(['\bf This pitt.grtheory.graph is ' st 'Eulerian'])
    V=[0 4;1 4;2 4;3 3.7;4 4;0 3;1 3.3;2 3.3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2;3 2;4 3;5 4;3 5;6 1;2 7;8 2;3 8;9 4;9 5;...
       10 4;10 5;7 6;8 7;6 9;20 15;11 16;21 18;19 23;...
       8 9;10 9;11 6;7 12;14 9;15 10;12 11;13 12;13 14;...
       14 13;15 14;16 11;12 17;13 18;20 15;17 16;17 18;...
       19 18;19 20;21 16;17 22;18 22;18 23;19 24;20 25;...
       21 22;22 21;23 24;24 23;24 25];
    [eu,cEu]=pitt.grtheory.grIsEulerian(E);
    switch eu,
      case 1,
        st='';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      case 0.5,
        st='semi-';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      otherwise,
        st='not ';
        E=[E(:,1:2), [1:size(E,1)]'];
    end
    pitt.grtheory.grPlot(V,[E,[1:size(E,1)]'],'g','');
    title(['\bf This pitt.grtheory.graph is ' st 'Eulerian'])
    V=[0 4;1 4;2 4;3 3.7;4 4;0 3;1 3.3;2 3.3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2;3 2;4 3;5 4;3 5;6 1;2 7;8 2;3 8;9 4;9 5;...
       10 4;10 5;7 6;8 7;6 9;20 15;11 16;21 18;19 23;...
       8 9;10 9;11 6;7 12;14 9;15 10;12 11;13 12;13 14;...
       14 13;15 14;16 11;12 17;13 18;20 15;17 16;17 18;...
       19 18;19 20;21 16;17 22;18 22;18 23;19 24;20 25;...
       21 22;22 21;23 24;24 25];
    [eu,cEu]=pitt.grtheory.grIsEulerian(E);
    switch eu,
      case 1,
        st='';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      case 0.5,
        st='semi-';
        E=[E(cEu,1:2), [1:size(E,1)]'];
      otherwise,
        st='not ';
        E=[E(:,1:2), [1:size(E,1)]'];
    end
    pitt.grtheory.grPlot(V,[E,[1:size(E,1)]'],'g','');
    title(['\bf This pitt.grtheory.graph is ' st 'Eulerian'])
  case 12, % pitt.grtheory.grIsomorph test
    disp('The pitt.grtheory.grIsomorph test')
    V1=[0 0;1 0;2 0;0 1;1 1;2 1];
    E1=[1 4;1 5;1 6;2 4;2 5;2 6;3 4;3 5;3 6];
    V2=[cos([1:6]*pi/3);sin([1:6]*pi/3)]';
    E2=[1 2;2 3;3 4;4 5;5 6;6 1;1 4;2 5;3 6];
    pitt.grtheory.grPlot(V1,E1,'g','%d','');
    title('\bfThe pitt.grtheory.graph 1')
    pitt.grtheory.grPlot(V2,E2,'g','%d','');
    title('\bfThe pitt.grtheory.graph 2')
    [IsIsomorph,Permut]=pitt.grtheory.grIsomorph(E1,E2);
    if IsIsomorph,
      disp('Their pitt.grtheory.graphs are isomorphic:')
      fprintf('Number of vertexes:\n  V1   V2\n')
      fprintf(' %3.0f  %3.0f\n',[[1:max(E1(:))];Permut'])
    else
      disp('Their pitt.grtheory.graphs are nonisomorphic')
    end;
  case 13, % pitt.grtheory.grMaxComSu test
    disp('The pitt.grtheory.grMaxComSu test')
    V=[0 0 2;1 1 3;1 0 3;1 -1 4;2 1 1;2 0 2;2 -1 3;3 1 4;...
       3 0 5;3 -1 1;4 0 5]; % vertexes coordinates and weights
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V(:,1:2),E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed vertexes')
    nMS=pitt.grtheory.grMaxComSu(E); % the maximal complete supitt.grtheory.graph
    fprintf('Number of vertexes on the maximal complete supitt.grtheory.graph = %d\n',...
      length(nMS));
    disp('In a maximal complete supitt.grtheory.graph is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
    nMS=pitt.grtheory.grMaxComSu(E,V(:,3)); % the weightd maximal complete supitt.grtheory.graph
    fprintf(['Number of vertexes on the weighed maximal complete supitt.grtheory.graph '...
      '= %d\n'],length(nMS));
    disp('In a weighed maximal complete supitt.grtheory.graph is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
  case 14, % pitt.grtheory.grMaxFlows test
    disp('The pitt.grtheory.grMaxFlows test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;...
       2 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 2;...
       9 10 2;8 11 5;9 11 4;10 11 4]; % arrows and weights
    s=1; % the network source
    t=11; % the network sink
    fprintf('The source of the net s=%d\nThe sink of the net t=%d\n',s,t)
    pitt.grtheory.grPlot(V,E,'d','','%d'); % the initial dipitt.grtheory.graph
    title('\bfThe dipitt.grtheory.graph of the net')
    [v,mf]=pitt.grtheory.grMaxFlows(E,s,t); % the maximal flow
    disp('The solution of the maximal flows problem')
    disp('  N arrow       flow')
    fprintf('   %2.0f      %12.8f\n',[[1:length(v)];v'])
    fprintf('The maximal flow =%12.8f\n',mf)
    pitt.grtheory.grPlot(V,[E(:,1:2),v],'d','','%6.4f'); % plot the dipitt.grtheory.graph
    title('\bfThe flows on the arrows')
  case 15, % pitt.grtheory.grMaxMatch test
    disp('The pitt.grtheory.grMaxMatch test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;...
       2 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 2;...
       9 10 2;8 11 5;9 11 4;10 11 4]; % arrows and weights
    pitt.grtheory.grPlot(V,E,'g','','%d'); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    nMM=pitt.grtheory.grMaxMatch(E(:,1:2)); % the maximal matching
    fprintf('Number of edges on the maximal matching = %d\n',...
      length(nMM));
    disp('In a maximal matching is the edges with numbers:');
    fprintf('%d  ',nMM);
    fprintf('\nThe total weight = %d\n',sum(E(nMM,3)));
    pitt.grtheory.grPlot(V,E(nMM,:),'g','','%d'); % the maximal matching
    title('\bfThe maximal matching')
    nMM=pitt.grtheory.grMaxMatch(E); % the weighed maximal matching
    fprintf('Number of edges on the weighed maximal matching = %d\n',...
      length(nMM));
    disp('In a weighed maximal matching is the edges with numbers:');
    fprintf('%d  ',nMM);
    fprintf('\nThe total weight = %d\n',sum(E(nMM,3)));
    pitt.grtheory.grPlot(V,E(nMM,:),'g','','%d'); % the weighed maximal matching
    title('\bfThe weighed maximal matching')
  case 16, % pitt.grtheory.grMaxStabSet test
    disp('The pitt.grtheory.grMaxStabSet test')
    V=[0 0 2;1 1 3;1 0 3;1 -1 4;2 1 1;2 0 2;2 -1 3;3 1 4;...
       3 0 5;3 -1 1;4 0 5]; % vertexes coordinates and weights
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V(:,1:2),E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    nMS=pitt.grtheory.grMaxStabSet(E); % the maximal stable set
    fprintf('Number of vertexes on the maximal stable set = %d\n',...
      length(nMS));
    disp('In a maximal stable set is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed vertexes')
    nMS=pitt.grtheory.grMaxStabSet(E,V(:,3)); % the weightd maximal stable set
    fprintf(['Number of vertexes on the weighed maximal stable set '...
      '= %d\n'],length(nMS));
    disp('In a weighed maximal stable set is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
  case 17, % pitt.grtheory.grMinAbsEdgeSet test
    disp('The pitt.grtheory.grMinAbsEdgeSet test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;...
       2 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 4;...
       9 10 5;8 11 5;9 11 4;10 11 4]; % arrows and weights
    pitt.grtheory.grPlot(V,E,'g','','%d'); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    nMS=pitt.grtheory.grMinAbsEdgeSet(E(:,1:2)); % the minimal absorbant set of edges
    fprintf('Number of edges on the minimal absorbant set = %d\n',...
      length(nMS));
    disp('In a minimal absorbant set is the edges with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(E(nMS,3)));
    pitt.grtheory.grPlot(V,E(nMS,:),'g','','%d'); % the minimal absorbant set of edges
    title('\bfThe minimal absorbant set of edges')
    nMS=pitt.grtheory.grMinAbsEdgeSet(E); % the minimal weighed absorbant set of edges
    fprintf('Number of edges on the minimal weighed absorbant set = %d\n',...
      length(nMS));
    disp('In a minimal weighed absorbant set is the edges with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(E(nMS,3)));
    pitt.grtheory.grPlot(V,E(nMS,:),'g','','%d'); % the minimal weighed absorbant set of edges
    title('\bfThe minimal weighed absorbant set of edges')
  case 18, % pitt.grtheory.grMinAbsVerSet test
    disp('The pitt.grtheory.grMinAbsVerSet test')
    V=[0 0 2;1 1 3;1 0 3;1 -1 4;2 1 1;2 0 2;2 -1 3;3 1 4;...
       3 0 5;3 -1 1;4 0 5]; % vertexes coordinates and weights
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V(:,1:2),E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph')
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed vertexes')
    nMS=pitt.grtheory.grMinAbsVerSet(E); % the minimal absorbant set of vertexes
    fprintf('Number of vertexes on the minimal absorbant set = %d\n',...
      length(nMS));
    disp('In a minimal absorbant set is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
    nMS=pitt.grtheory.grMinAbsVerSet(E,V(:,3)); % the weightd minimal absorbant set of vertexes
    fprintf(['Number of vertexes on the weighed minimal absorbant set '...
      '= %d\n'],length(nMS));
    disp('In a weighed minimal absorbant set is the vertexes with numbers:');
    fprintf('%d  ',nMS);
    fprintf('\nThe total weight = %d\n',sum(V(nMS,3)));
  case 19, % pitt.grtheory.grMinCutSet test
    disp('The pitt.grtheory.grMinCutSet test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;...
       2 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 2;...
       9 10 2;8 11 5;9 11 4;10 11 4]; % arrows and weights
    s=1; % the network source
    t=11; % the network sink
    fprintf('The source of the net s=%d\nThe sink of the net t=%d\n',s,t)
    pitt.grtheory.grPlot(V,E,'d','','%d'); % the initial dipitt.grtheory.graph
    title('\bfThe dipitt.grtheory.graph of the net')
    [nMCS,mf]=pitt.grtheory.grMinCutSet(E,s,t); % the minimal cut-set
    fprintf('The first minimal cut-set include arrows:');
    fprintf('  %d',nMCS);
    fprintf(['\nThe maximal flow through '...
      'each minimal cut-set = %12.6f\n'],mf)
    pitt.grtheory.grPlot(V,E(setdiff(1:size(E,1),nMCS),:),'d','','%d');
    title('\bfThe dipitt.grtheory.graph without first minimal cut-set')
  case 20, % pitt.grtheory.grMinEdgeCover test
    disp('The pitt.grtheory.grMinEdgeCover test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates and weights
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;2 6 2;3 6 5;...
       3 7 2;4 7 3;5 6 1;6 7 1;5 8 5;6 8 2;6 9 3;7 9 2;...
       7 10 3;8 9 2;9 10 2;8 11 5;9 11 4;10 11 4]; % edges and weights
    pitt.grtheory.grPlot(V,E,'g',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    nMC=pitt.grtheory.grMinEdgeCover(E(:,1:2)); % the minimal edge covering
    fprintf('Number of edges on the minimal edge covering = %d\n',...
      length(nMC));
    disp('In a minimal edge cover is the edges with numbers:');
    fprintf('%d  ',nMC);
    fprintf('\nThe total weight = %d\n',sum(E(nMC,3)));
    pitt.grtheory.grPlot(V,E(nMC,:),'g',''); % the minimal edge covering
    title('\bfThe minimal edge covering')
    nMC=pitt.grtheory.grMinEdgeCover(E); % the weighed minimal edge covering
    fprintf('Number of edges on the weighed minimal edge covering = %d\n',...
      length(nMC));
    disp('In a weighed minimal edge cover is the edges with numbers:');
    fprintf('%d  ',nMC);
    fprintf('\nThe total weight = %d\n',sum(E(nMC,3)));
    pitt.grtheory.grPlot(V,E(nMC,:),'g',''); % the weighed minimal edge covering
    title('\bfThe weighed minimal edge covering')
  case 21, % pitt.grtheory.grMinSpanTree test
    disp('The pitt.grtheory.grMinSpanTree test')
    V=[0 4;1 4;2 4;3 4;4 4;0 3;1 3;2 3;3 3;4 3;...
       0 2;1 2;2 2;3 2;4 2;0 1;1 1;2 1;3 1;4 1;...
       0 0;1 0;2 0;3 0;4 0];
    E=[1 2 1;3 2 2;4 3 3;5 4 4;6 1 5;2 7 6;8 2 7;3 8 8;...
       9 4 9;9 5 8;10 5 7;7 6 6;8 7 5;8 9 4;10 9 3;11 6 2;...
       7 12 1;13 8 2;14 9 3;15 10 4;12 11 5;13 12 6;13 14 7;...
       14 13 8;5 5 10;15 14 9;16 11 8;12 17 7;13 18 6;...
       20 15 5;17 16 4;17 18 3;18 17 2;19 18 1;19 20 2;...
       5 5 8; 21 16 3;17 22 4;18 22 5;22 18 6;18 23 7;...
       19 24 8;20 25 9;21 22 8;22 21 7;23 24 6;10 10 8;...
       24 23 5;24 25 4];
    pitt.grtheory.grPlot(V,E); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed edges')
    nMST=pitt.grtheory.grMinSpanTree(E(:,1:2)); % the spanning tree
    fprintf('Number of edges on the spanning tree = %d\n',length(nMST));
    fprintf('The total weight = %d\n',sum(E(nMST,3)));
    pitt.grtheory.grPlot(V,E(nMST,:)); % the spanning tree
    title('\bfThe spanning tree')
    nMST=pitt.grtheory.grMinSpanTree(E); % the minimal spanning tree
    fprintf('Number of edges on the minimal spanning tree = %d\n',...
      length(nMST));
    fprintf('The total weight = %d\n',sum(E(nMST,3)));
    pitt.grtheory.grPlot(V,E(nMST,:)); % the minimal spanning tree
    title('\bfThe minimal spanning tree')
  case 22, % pitt.grtheory.grMinVerCover test
    disp('The pitt.grtheory.grMinVerCover test')
    V=[0 0 2;1 1 3;1 0 3;1 -1 4;2 1 1;2 0 2;2 -1 3;3 1 4;...
       3 0 7;3 -1 1;4 0 5]; % vertexes coordinates and weights
    E=[1 2;1 3;1 4;2 3;3 4;2 5;2 6;3 6;3 7;4 7;6 5;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;9 10;8 11;9 11;10 11]; % edges
    pitt.grtheory.grPlot(V,E,'g','%d',''); % the initial pitt.grtheory.graph
    title('\bfThe initial pitt.grtheory.graph with weighed vertexes')
    nMC=pitt.grtheory.grMinVerCover(E); % the minimal vertex cover
    fprintf('Number of vertexes on the minimal vertex cover = %d\n',...
      length(nMC));
    disp('In a minimal vertex cover is the vertexes with numbers:');
    fprintf('%d  ',nMC);
    fprintf('\nThe total weight = %d\n',sum(V(nMC,3)));
    pitt.grtheory.grPlot(V(nMC,:)); % the solution of the MinVerCover problem
    title('\bfThe minimal vertex cover')
    nMC=pitt.grtheory.grMinVerCover(E,V(:,3)); % the weightd minimal vertex cover
    fprintf(['Number of vertexes on the weighed minimal vertex cover '...
      '= %d\n'],length(nMC));
    disp('In a weighed minimal vertex cover is the vertexes with numbers:');
    fprintf('%d  ',nMC);
    fprintf('\nThe total weight = %d\n',sum(V(nMC,3)));
    pitt.grtheory.grPlot(V(nMC,:)); % the solution of the weighed MinVerCover problem
    title('\bfThe weighed minimal vertex cover')
  case 23, % pitt.grtheory.grPERT test
    disp('The pitt.grtheory.grPERT test')
    V=[1 1;0 0;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       4 0;3 -1;3 0]; % events coordinates
    E=[2 1 5;2 3 5;2 4 5;1 3 2;3 4 2;1 5 3;...
       1 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 11 3;7 11 2;7 10 3;8 11 2;...
       11 10 2;8 9 5;11 9 4;10 9 4]; % works and their times
    pitt.grtheory.grPlot(V,E,'d','%d','%d');
    title('\bfThe schema of project')
    [CrP,Ts,Td]=pitt.grtheory.grPERT(E);
    pitt.grtheory.grPlot([V Ts'],[CrP(1:end-1);CrP(2:end)]','d','%d','');
    title('\bfThe critical path and start times for events')
    pitt.grtheory.grPlot([V Ts'],[E(:,1:2) Td],'d','%d','%d')
    title('\bfThe start times for events and delay times for works')
  case 24, % pitt.grtheory.grPlot test
    disp('The pitt.grtheory.grPlot test')
    V=[0 0 2;1 1 3;1 0 3;1 -1 4;2 1 1;2 0 2;2 -1 3;3 1 4;...
       3 0 5;3 -1 1;4 0 5]; % vertexes coordinates and weights
    E=[1 2 5;1 1 2;1 1 5;2 2 3;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;2 6 2;3 6 5;3 7 2;...
       4 7 3;5 6 1;6 7 1;5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 2;...
       9 10 2;8 11 5;9 11 4;10 11 4;1 2 8;1 3 4;1 3 5;1 3 6]; % edges (arrows) and weights
    pitt.grtheory.grPlot(V(:,1:2),E,'d');
    title('\bfThe dipitt.grtheory.graph with weighed multiple arrows and loops')
    pitt.grtheory.grPlot(V,E(:,1:2),[],'%d','');
    title('\bfThe pitt.grtheory.graph with weighed vertexes without edges numeration')
    pitt.grtheory.grPlot(V(:,1:2));
    title('\bfThe disconnected pitt.grtheory.graph')
    pitt.grtheory.grPlot([],fullfact([5 5]),'d','','',0.8)
    title('\bfThe directed clique\rm \itK\rm_5')
  case 25, % pitt.grtheory.grShortPath test
    disp('The pitt.grtheory.grShortPath test')
    V=[0 0;1 1;1 0;1 -1;2 1;2 0;2 -1;3 1;...
       3 0;3 -1;4 0]; % vertexes coordinates
    E=[1 2 5;1 3 5;1 4 5;2 3 2;3 4 2;2 5 3;...
       2 6 2;3 6 5;3 7 2;4 7 3;5 6 1;6 7 1;...
       5 8 5;6 8 2;6 9 3;7 9 2;7 10 3;8 9 2;...
       9 10 2;8 11 5;9 11 4;10 11 4]; % arrows and weights
    s=1; % the network source
    t=11; % the network sink
    fprintf('The source of the net s=%d\nThe sink of the net t=%d\n',s,t)
    pitt.grtheory.grPlot(V(:,1:2),E,'d','','%d');
    title('\bfThe dipitt.grtheory.graph with weighed edges')
    [dSP,sp]=pitt.grtheory.grShortPath(E,s,t);
    disp('The shortest paths between all vertexes:');
    fprintf('    %2.0f',1:size(dSP,2));
    fprintf('\n');
    for k1=1:size(dSP,1),
      fprintf('%2.0f',k1)
      fprintf('%6.2f',dSP(k1,:))
      fprintf('\n')
    end
    pitt.grtheory.grPlot(V(:,1:2),[sp(1:end-1);sp(2:end)]','d','%d','')
    title(['\bfThe shortest path from vertex ' ...
      num2str(s) ' to vertex ' num2str(t)])
  case 26, % pitt.grtheory.grShortVerPath test
    disp('The pitt.grtheory.grShortVerPath test')
    V=[-8 0; -8 1; -7 0; -7 1; -6 0; -6 1; -6 2; -5 2; -4 1; -3 2; 
       -3 1; -2 0; -2 1; -1 2;  0 1;  0 2];
    E=[1 3; 2 4; 3 5; 3 6; 4 6; 4 7; 5 9; 6 9; 7 8; 8 9; 8 10; 9 11; 
       11 12; 11 13; 10 14; 13 14; 13 15; 14 16];
    V(:,3)=[5;9;1;2;2;4;9;9;1;4;6;5;7;7;3;5];
    pitt.grtheory.grPlot(V,E,'d','','',3);
    title('\bfThe dipitt.grtheory.graph with weighted vertices')
    [dMWP,ssp]=pitt.grtheory.grShortVerPath(E,V(:,3)); % solution
    E1=repmat(dMWP(2:end-1),2,1);
    E2=[dMWP(1); E1(:); dMWP(end)];
    E3=reshape(E2,2,length(E2)/2)';
    disp('The Minimal Vertexes Path:')
    fprintf('%d   ',dMWP)
    fprintf('\nThe total weight = %d\n',ssp)
    pitt.grtheory.grPlot(V,E3,'d','','',3)
    title('\bfThe Minimal Vertexes Weight Path')
  case 27, % pitt.grtheory.grTranClos test
    disp('The pitt.grtheory.grTranClos test')
    V=[0 0;1 1;1.2 0.2;1 -1;2 1.2;2.2 0.2;2 -1.2;3 1;...
       3.2 -0.2;3 -1;4 0]; % vertexes coordinates
    E=[1 2;1 3;1 4;2 3;3 4;2 5;...
       2 6;3 6;3 7;4 7;5 6;6 7;...
       5 8;6 8;6 9;7 9;7 10;8 9;...
       9 10;8 11;9 11;10 11]; % arrows
    pitt.grtheory.grPlot(V,E,'d','%d','');
    title('\bfThe initial dipitt.grtheory.graph')
    Ecl=pitt.grtheory.grTranClos(E); % solution
    pitt.grtheory.grPlot(V,Ecl,'d','%d','');
    title('\bfThe transitive closure for the dipitt.grtheory.graph E')
  case 28, % pitt.grtheory.grTravSale test
    disp('The pitt.grtheory.grTravSale test')
    C=[0 3 7 4 6 4;4 0 3 7 8 5;6 9 0 3 2 1;...
       8 6 3 0 9 8;3 7 4 6 0 4;4 5 8 7 2 0];
    disp('The distances between cities:')
    fprintf('  %18.0f',1:size(C,2))
    fprintf('\n')
    for k1=1:size(C,1),
      fprintf('%2.0f',k1)
      fprintf('%20.12f',C(k1,:))
      fprintf('\n')
    end
    [pTS,fmin]=pitt.grtheory.grTravSale(C);
    disp('The order of cities:')
    fprintf('%d   ',pTS)
    fprintf('\nThe minimal way =%3.0f\n',fmin)
  otherwise,
    error('Select the test')
end