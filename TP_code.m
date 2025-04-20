clear all; clc;

% Pour La longitude $
    g=9.8;
    Ro1max=g;
    Ro1min=0.21*g; 

    A1Long=[0 1 0 0; 0 0 -Ro1max 0; 0 0 0 1; 0 0 0 0];
    A2Long=[0 1 0 0; 0 0 -Ro1min 0; 0 0 0 1; 0 0 0 0];
    B1Long=[0;1;0;1];
    C1Long=[1;0; 0; 1];
    %% Definition de la matrice deinit positive
        Q1Long=sdpvar(4,4,'symmetric');

    %% Calcul des gains
        F1Long=sdpvar(1,4);
        F2Long=sdpvar(1,4);

     %% Definition de la condition de lemme congruence
     
        %%%Remarque : Pour chaque pair (Ai, Bi) on a 2 
        % équations LMI. Dans notre cas on aura deux équations LMI car B_1=B_2

        LM1=Q1Long*A1Long'+A1Long*Q1Long-B1Long*F1Long-F1Long'*B1Long';
        LM2=Q1Long*A2Long'+A2Long*Q1Long-B1Long*F2Long-F2Long'*B1Long';
        LM5=Q1Long*A1Long'+A1Long*Q1Long-B1Long*F2Long-F2Long'*B1Long';
        LM6=Q1Long*A2Long'+A2Long*Q1Long-B1Long*F1Long-F1Long'*B1Long';

        %% Definition de la condition de lemme congruence pour L
     
        %%%Remarque : Pour chaque pair (Ai, Bi) on a 2 
        % équations LMI. Dans notre cas on aura deux équations LMI car B_1=B_2
        F3Long=sdpvar(1,4);
        F4Long=sdpvar(1,4);
        LM3=A1Long'*Q1Long+ Q1Long*A1Long-C1Long'*F3Long'-F3Long*C1Long;
        LM4=A2Long'*Q1Long+ Q1Long*A2Long-C1Long'*F4Long'-F4Long*C1Long;
     %% Resolution du problème d'optimisation

        ProbemLong=[Q1Long>=(1e-5)*eye(4),LM1<=(-1e-5)*eye(4),LM2<=(-1e-5)*eye(4),LM5<=(-1e-5)*eye(4),LM6<=(-1e-5)*eye(4)];
        SolutionLong=solvesdp(ProbemLong);
        Q1Long = double(Q1Long);
        P1Long = inv(Q1Long);
        K1Long = double(F1Long)*P1Long;
        K2Long = double(F2Long)*P1Long;


      %% Resolution du problème d'optimisation L

% %         ProbemLong2=[Q1Long>=(1e-5)*eye(4),LM3<=(-1e-5)*eye(4),LM4<=(-1e-5)*eye(4)];
% %         SolutionLong2=solvesdp(ProbemLong);
% %         Q1Long = double(Q1Long);
% %         L1Long = P1Long*double(F3Long);
% %         L2Long = P1Long*double(F4Long);

     %% Etude de la stabilité.
     disp(['La valeur du K1Long est:' num2str(K1Long)]);
     disp(['La valeur du K2Long est:' num2str(K2Long)]);
     %%valeurs Propres
     modLong1=eig(A1Long-B1Long*K1Long);
     modLong2=eig(A2Long-B1Long*K2Long);
     disp('La valeur du mode 1 est:')
      disp(modLong1);
     disp('La valeur du mode 2 est:')
     disp(modLong2);


        

     % Pour La latitude $
    g=9.8;
    Ro2max=g;
    Romin=0.21*g; 

    A1Lat=[0 1 0 0; 0 0 Ro2max 0; 0 0 0 1; 0 0 0 0];
    A2Lat=[0 1 0 0; 0 0 Ro1min 0; 0 0 0 1; 0 0 0 0];
    B1Lat=[0;1;0;1];
    
    %% Definition de la matrice deinit positive
        Q1Lat=sdpvar(4,4,'symmetric');

    %% Calcul des gains
        F1Lat=sdpvar(1,4);
        F2Lat=sdpvar(1,4);

     %% Definition de la condition de lemme congruence
     
        %%%Remarque : Pour chaque pair (Ai, Bi) on a 2 
        % équations LMI. Dans notre cas on aura deux équations LMI car B_1=B_2

        LM1_Lat=Q1Lat*A1Lat'+A1Lat*Q1Lat-B1Lat*F1Lat-F1Lat'*B1Lat';
        LM2_Lat=Q1Lat*A2Lat'+A2Lat*Q1Lat-B1Lat*F2Lat-F2Lat'*B1Lat';
        LM3_Lat=Q1Lat*A1Lat'+A1Lat*Q1Lat-B1Lat*F2Lat-F2Lat'*B1Lat';
        LM4_Lat=Q1Lat*A2Lat'+A2Lat*Q1Lat-B1Lat*F1Lat-F1Lat'*B1Lat';

     %% Resolution du problème d'optimisation

        ProbemLat=[Q1Lat>=(1e-8)*eye(4),LM1_Lat<=(-1e-8)*eye(4),LM2_Lat<=(-1e-8)*eye(4)];
        SolutionLat=solvesdp(ProbemLat);
        Q1Lat = double(Q1Lat);
        P1Lat = inv(Q1Lat);
        K1Lat = double(F1Lat)*P1Lat;
        K2Lat = double(F2Lat)*P1Lat;

     %% Etude de la stabilité.
     disp(['La valeur du K1Lat est:' num2str(K1Lat)]);
     disp(['La valeur du K2Lat est:' num2str(K2Lat)]);
     %%valeurs Propres
     modLat1=eig(A1Lat-B1Lat*K1Lat);
     modLat2=eig(A2Lat-B1Lat*K2Lat);
     disp('La valeur du mode 1 est:')
      disp(modLat1);
     disp('La valeur du mode 2 est:')
     disp(modLat2);


        
        

        
