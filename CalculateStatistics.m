clear;
RestLength = 10;
data = {};
for L= [20,50,80] 
    for Ninit = 200 
        for Rcross = 1 
            for kappaR1 = [1,10,100] 
                for kappaR2 = 0
                    for ensemble= 0 
                        for Cross_angle = 0
                            for Rdecross = 1 
                                for kappa = [0.1] 
                                    for InitialSize = 100
                                        for InitialAlign = 0
                                            for Edecross = 100
                                                for geo = {'cubecylinder'}
                                                    for Max_distance = 20 
                                                        for Bkappa = [1,2,10] 
                                                            for radiuscylinder = [30,50,70,90] 
                                                                pattern1=['./data/',...
                                                                    'Trajectory_',num2str(ensemble),...
                                                                    '_Ninit_',num2str(Ninit),...
                                                                    '_Linit_',num2str(L),...
                                                                    '_Rcross_',num2str(Rcross),...
                                                                    '_Rdecross_',num2str(Rdecross),...
                                                                    '_Edecross_',num2str(Edecross),...
                                                                    '_kappa_',num2str(kappa),...
                                                                    '_kappaR1_',num2str(kappaR1),...
                                                                    '_kappaR2_',num2str(kappaR2),...
                                                                    '_Crossangle_',num2str(Cross_angle),...
                                                                    '_geometry_',geo{1},...
                                                                    '_InitialSize_',num2str(InitialSize),...
                                                                    '_InitialAlign_',num2str(InitialAlign),...
                                                                    '_MaxDistance_',num2str(Max_distance),...
                                                                    '_RadiusCylinder_',num2str(radiuscylinder),...
                                                                    '_Bkappa_',num2str(Bkappa),...
                                                                    '_date_*.txt'];
                                                                  %  '*.txt']
                                                                
                                                                pattern2=['./data/',...
                                                                    'Connection_',num2str(ensemble),...
                                                                    '_Ninit_',num2str(Ninit),...
                                                                    '_Linit_',num2str(L),...
                                                                    '_Rcross_',num2str(Rcross),...
                                                                    '_Rdecross_',num2str(Rdecross),...
                                                                    '_Edecross_',num2str(Edecross),...
                                                                    '_kappa_',num2str(kappa),...
                                                                    '_kappaR1_',num2str(kappaR1),...
                                                                    '_kappaR2_',num2str(kappaR2),...
                                                                    '_Crossangle_',num2str(Cross_angle),...
                                                                    '_geometry_',geo{1},...
                                                                    '_InitialSize_',num2str(InitialSize),...
                                                                    '_InitialAlign_',num2str(InitialAlign),...
                                                                    '_MaxDistance_',num2str(Max_distance),...
                                                                    '_RadiusCylinder_',num2str(radiuscylinder),...
                                                                    '_Bkappa_',num2str(Bkappa),...
                                                                    '_date_*.txt'];
                                                                  %  '*.txt'];
                                                                
                                                                fileLists = dir(pattern1);
                                                                fileLists2 = dir(pattern2);
                                                                
                                                                for kk = 1:numel(fileLists)
                                                                    filename = [fileLists(kk).folder,'/',fileLists(kk).name];
                                                                    filename2 = [fileLists2(kk).folder,'/',fileLists2(kk).name];
                                                                    traj = load(filename);
                                                                    conn = load(filename2);
                                                                    TT = unique(traj(:,1));
                                                                    jj = 0;
                                                                    NematicOrder = zeros(3,length(TT));
                                                                    LocalNematicOrder = zeros(3,length(TT));
                                                                    NumClusters = zeros(1,length(TT));
                                                                    NumNeighbours = zeros(1, length(TT));
                                                                    DistanceToCenter = zeros(1, length(TT));
                                                                    NearestDistanceToCenter = 2*radiuscylinder*ones(1, length(TT));
                                                                    AngularDistance = zeros(Ninit,length(TT));
                                                                    SpinDistance = zeros(Ninit,length(TT));
                                                                    for tt=TT'
                                                                        jj = jj+1;
                                                                        index1 = find(abs(traj(:,1)-tt)<1e-3);
                                                                        if ~isempty(conn)
                                                                            index2 = find(abs(conn(:,1)-tt)<1e-3);
                                                                        else
                                                                            index2 = [];
                                                                        end
                                                                        
                                                                        N = traj(index1,6:8);
                                                                        Cen = traj(index1,3:5);
                                                                        Len = traj(index1,2);
                                                                        AngularDistance(:,jj) = atan2(Cen(:,2),Cen(:,1));
                                                                        SpinDistance(:,jj) = atan2(N(:,2),N(:,1));
                                                                        S = N'*N/length(index1)-1/3*eye(3);
                                                                        NematicOrder(:,jj) = 1.5*eig(S);
                                                                        if isempty(index2)
                                                                            k1 = (1:Ninit)';
                                                                            k2 = (1:Ninit)';
                                                                        else
                                                                            k1 = [conn(index2,2);(1:Ninit)'];
                                                                            k2 = [conn(index2,4);(1:Ninit)'];
                                                                        end
                                                                        aj = accumarray([k1(:),k2(:)],1);
                                                                        A = aj + aj' - 2*eye(Ninit);
                                                                        Gr = graph(A);
                                                                        NumNeighbours(jj) = mean(Gr.degree);
                                                                        bins = conncomp(Gr);
                                                                        NumClusters(jj) = max(bins);
                                                                        for ii = 1:Ninit %max(bins)
                                                                            % Nc = N(bins == ii,:);
                                                                            Nc = N([ii;neighbors(Gr,ii)],:);
                                                                            Sc = Nc'*Nc/size(Nc,1)-1/3*eye(3);
                                                                            LocalNematicOrder(:,jj) = LocalNematicOrder(:,jj) + eig(1.5*Sc);%size(Nc,1)*eig(1.5*Sc);
                                                                            p1 = Cen(ii,:) + 0.5*Len(ii)*N(ii,:);
                                                                            p2 = Cen(ii,:) - 0.5*Len(ii)*N(ii,:);
                                                                            dis = DistBetween2Segment(p1, p2, [0,0,0], [0,0,-1000]);
                                                                            DistanceToCenter(jj) = DistanceToCenter(jj) + dis;
                                                                            NearestDistanceToCenter(jj) = min([dis,NearestDistanceToCenter(jj)]);
                                                                                
                                                                        end
                                                                        LocalNematicOrder(:,jj) = LocalNematicOrder(:,jj)/Ninit; % weighted by the number of filaments in a cluster
                                                                        DistanceToCenter(jj) = DistanceToCenter(jj)/Ninit;
                                                                    end
                                                                    
                                                                    if ~isempty(conn)
                                                                        ed = [-1;TT+0.1*(TT(2)-TT(1))];
                                                                        [Num,edges] = histcounts(conn(conn(:,15)>0,1),ed);
                                                                        index = find(conn(:,15)>0);
                                                                        subs = uint16(conn(index,1)/(TT(2)-TT(1))+1);
                                                                        %subs = subs - subs(1) + 1;
                                                                        theta11 = conn(index,7);
                                                                        theta12 = conn(index,9);
                                                                        theta21 = conn(index,11);
                                                                        theta22 = conn(index,13);
                                                                        Len = conn(index,14);
                                                                        Etor = .5*kappaR1.*(theta11.^2+theta21.^2)+...
                                                                            .5*kappaR2.*(theta12.^2+theta22.^2);
                                                                        Etra = .5*kappa.*(Len-RestLength).^2;
                                                                        if isempty(subs)
                                                                            subs = length(TT);
                                                                            Etor = 0;
                                                                            Etra = 0;
                                                                        else
                                                                            if max(subs) < length(TT)
                                                                                subs = [subs;length(TT)];
                                                                                Etor = [Etor;0];
                                                                                Etra = [Etra;0];
                                                                            end
                                                                        end
                                                                        EtorPerBond = accumarray(subs,Etor)./(Num'+eps);
                                                                        EtraPerBond = accumarray(subs,Etra)./(Num'+eps);
                                                                    else
                                                                        Num = zeros(1,length(TT));
                                                                        EtorPerBond = zeros(length(TT),1);
                                                                        EtraPerBond = zeros(length(TT),1);
                                                                    end
                                                                    
                                                                    Tstart = 1;
                                                                    S1 = max(NematicOrder,[],1);
                                                                    S2 = min(NematicOrder,[],1);
                                                                    Local_S1 = max(LocalNematicOrder, [], 1);
                                                                    Local_S2 = min(LocalNematicOrder, [], 1);
                                                                    mean_S1 = mean(S1(TT'>Tstart));
                                                                    std_S1 = std(S1(TT'> Tstart));
                                                                    mean_S2 = mean(S2(TT'> Tstart));
                                                                    std_S2 = std(S2(TT'>Tstart));
                                                                    mean_Local_S1 = mean(Local_S1(TT' > Tstart));
                                                                    std_Local_S1 = std(Local_S1(TT' > Tstart));
                                                                    mean_Local_S2 = mean(Local_S2(TT' > Tstart));
                                                                    std_Local_S2 = std(Local_S2(TT' > Tstart));
                                                                    
                                                                    mean_NumBonds = mean(Num(TT'>Tstart));
                                                                    std_NumBonds = std(Num(TT'>Tstart));
                                                                    
                                                                    mean_NumClusters = mean(NumClusters(TT' > Tstart));
                                                                    std_NumClusters = std(NumClusters(TT' > Tstart));
                                                                    
                                                                    mean_NumNeighbours = mean(NumNeighbours(TT' > Tstart));
                                                                    std_NumNeighbours = std(NumNeighbours(TT' > Tstart));
                                                                    
                                                                    mean_Etra = mean(EtraPerBond(TT'>Tstart));
                                                                    std_Etra = std(EtraPerBond(TT'>Tstart));
                                                                    mean_Etor = mean(EtorPerBond(TT'>Tstart));
                                                                    std_Etor = std(EtorPerBond(TT'>Tstart));
                                                                    
                                                                    mean_Dist = mean(DistanceToCenter(TT'>Tstart));
                                                                    std_Dist = std(DistanceToCenter(TT'>Tstart));
                                                                    
                                                                    mean_NearestDist = mean(NearestDistanceToCenter(TT' > Tstart));
                                                                    std_NearestDist = std(NearestDistanceToCenter(TT' > Tstart));
                                                                    
                                                                    AngularDistance = unwrap(AngularDistance,[],2);
                                                                    AngularVelocity = diff(AngularDistance,1,2)/(TT(2)-TT(1));
                                                                    Ang = mean(AngularVelocity,1);
                                                                    mean_AngularVelocity = mean(Ang(TT(1:end-1)' > Tstart));
                                                                    std_AngularVelocity = std(Ang(TT(1:end-1)' > Tstart));
                                                                    
                                                                    SpinDistance = unwrap(SpinDistance,[],2);
                                                                    SpinVelocity = diff(SpinDistance,1,2)/(TT(2)-TT(1));
                                                                    Spin = mean(SpinVelocity,1);
                                                                    mean_SpinVelocity = mean(Spin(TT(1:end-1)' > Tstart));
                                                                    std_SpinVelocity = std(Spin(TT(1:end-1)' > Tstart));
                                                                    
                                                                    el = {ensemble,Ninit,L,Rcross,Rdecross,...
                                                                        Edecross,kappa,kappaR1,kappaR2,...
                                                                        Cross_angle,geo{1},InitialSize,InitialAlign,...
                                                                        Max_distance,...
                                                                        radiuscylinder,...
                                                                        Bkappa,...
                                                                        mean_S1,std_S1,mean_S2,std_S2,...
                                                                        mean_Local_S1,std_Local_S1,mean_Local_S2,std_Local_S2,...
                                                                        mean_NumBonds,std_NumBonds,...
                                                                        mean_NumClusters,std_NumClusters,...
                                                                        mean_NumNeighbours, std_NumNeighbours,...
                                                                        mean_Etra,std_Etra,...
                                                                        mean_Etor,std_Etor,...
                                                                        mean_Dist,std_Dist,...
                                                                        mean_NearestDist, std_NearestDist,...
                                                                        mean_AngularVelocity, std_AngularVelocity,...
                                                                        mean_SpinVelocity, std_SpinVelocity
                                                                        };
                                                                    data = [data;el];
                                                                    
                                                                    
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
data
T = cell2table(data,...
    'VariableNames',...
    {'Ensemble', 'Ninit', 'Linit', 'Rcross', 'Rdecross',...
      'Edecross','kappa','kappaR1','kappaR2',...
      'Crossangle','geo','InitialSize','InitialAlign',...
      'MaxDistance',...
      'CylinderRadius',...
      'Bkappa',...
      'mean_S1','std_S1','mean_S2','std_S2',...
      'mean_Local_S1','std_Local_S1','mean_Local_S2','std_Local_S2',...
      'mean_NumBonds','std_NumBonds',...  
      'mean_NumClusters','std_NumClusters',...
      'mean_NumNeighbours', 'std_NumNeighbours',...
      'mean_Etra','std_Etra',...
      'mean_Etor','std_Etor',...
      'mean_Dist','std_Dist',...
      'mean_NearestDist','std_NearestDist',...
      'mean_AngularVelocity', 'std_AngularVelocity',...
      'mean_SpinVelocity', 'std_SpinVelocity'});
  
writetable(T,['./figures/Statistics_',date,'.csv']);
