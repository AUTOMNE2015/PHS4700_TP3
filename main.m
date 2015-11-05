function main
    hold on
    fprintf('start\n');
    sol = zeros(8);
    sol = Devoir3([0 0 0], [6.85, 0.0, 6.85], 0.66 );
    celldisp(sol)
    pointsBalle = sol{4};
    
    
    x1 = pointsBalle(:, 1);
    y1 = pointsBalle(:, 2);
    z1 = pointsBalle(:, 3);
    scatter3(x1,y1,z1);

    pointsBoite = sol{5};
    x2 = pointsBoite(:, 1);
    y2 = pointsBoite(:, 2);
    z2 = pointsBoite(:, 3);
    scatter3(x2,y2,z2);
   % option 1
   
   % option 1 end

   % option 2
   
   % option 2 end

   % option 3

   % option 3 end
end

function y = dt()
    y = 0.1;
end

function y = Pos0Balle()
    y = [0 0 2];
end


function y = RayonBalle()
     y = 0.0335;
end

function y = MasseBalle()
    y = 0.058;
end

function y = AireBalle()
    y = pi*RayonBalle()*RayonBalle();
end

function y = AireBoite()
    y = RayonBoite()*RayonBoite() + HauteurBoite()*HauteurBoite()
end


function y = Pos0Boite()
    y = [3 0 10];
end

function y = RayonBoite()
    y = 0.05;
end

function y = HauteurBoite()
    y = 0.15;
end

function y = MasseBoite()
    y = 0.075;
end


function y = FrottementAir()
    y = 0.1;
end

function y = Epsilon()
    y = 0.5
end

%wboitei vitesse angulaire initiale de la boˆ?te.
% vballei vitesse initiale du centre de masse de la balle.
% tballe temps tl ou la balle est lanc´ee.
% [Dev vbaf vbof rba rbo tc]=Devoir3(wboitei,vballei,tballe)
function y = Devoir3(wboitei,vballei,tballe)
    deltaT = dt();
    qSolBalle = zeros(8);
    qSolBoite = zeros(7);
    posBalle = Pos0Balle();
    posBoite = Pos0Boite();
    
    qSolBalle(1,:) = [0 vballei(1) vballei(2) vballei(3) posBalle(1) posBalle(2) posBalle(3) 0];
    qSolBoite(1,:) = [0 0 0 0 posBoite(1) posBoite(2) posBoite(3)];
    
    i = 1; %nb iterations
    result = 0; %is there collision
    while( i < 10 && result == 0 )
        if(qSolBoite(i,1) >= tballe)
            % Calculer la balle avec precision et imposer son Deltat a la
            % boite
            qSolBalle(i+1,:) = SEDRK4c(qSolBalle(i,2:7), qSolBalle(i,1),  deltaT, 0.25, @fonctionGballe);
            qSolBoite(i+1,:) = SEDRK4cImprecis(qSolBoite(i,2:7), qSolBoite(i,1), qSolBalle(i+1, 8), 0.25, @fonctionGboite);
        else
            % Calculer boite seulement, la balle ne bouge pas au debut
            qSolBoite(i+1,:) = SEDRK4c(qSolBoite(i,2:7), qSolBoite(i,1), deltaT, 0.25, @fonctionGboite);
            qSolBalle(i+1,:) = [qSolBoite(i+1,1) vballei(1) vballei(2) vballei(3) posBalle(1) posBalle(2) posBalle(3) deltaT];
        end 
        i = i + 1;
        result = CollisionDetect();
    end
    
    %qsol(1) est le temps
    % result-1 = 0 si sol, sinon = 1 si collision
    if(result-1 == 1)
        %calcul vitesse finaux apres collision
    end
    
    %qfinal = [ 0 0 0 0 0 0];
   % qfinal = cell(1,6);
    vbaf = [0 0 0; 0 0 0];
    vbof = [0 0 0; 0 0 0];
    rba = qSolBalle(:,5:7);
    rbo = qSolBoite(:,5:7);
    ti = qSolBalle(:,1);
    tc = 0;
    qfinal = {result-1 vbaf vbof rba rbo tc};
%     qfinal(2) = {vbaf};
%     qfinal(3) = {vbof};
%     qfinal(4) = {rba};
%     qfinal(5) = {rbo};
%     qfinal(6) = {tc};
    
    %qfinal = [result-1 vbaf vbof rba rbo tc];
    y = qfinal; 
end

function y = CollisionDetect()
    y = 0;
end
function [t2 qs1 qs2 qs3 qs4 qs5 qs6]= SEDRK4cImprecis (q0 ,t0 ,Deltat ,Err , fonctiong )
    % Solution ED dq/dt= fonctiong (q)
    % Methode de Runge - Kutta d’ ordre 4
    % Contrôle d’ erreur et interpolation Richardson
    % qs : vecteur final (v,r)
    % q0 : vecteur initial (v,r)
    % Deltat : intervalle de temps
    % Err : Precision requise ( erreur relative )
    % fonctiong : membre de droite de ED.
    % maxit : nombre maximum itérations

    % Solution sur intervalle complet
    qs= SEDRK4t0 (q0 ,t0 ,Deltat , fonctiong );
    t2 = t0+ Deltat;
    qs1 = qs(1);
    qs2 = qs(2);
    qs3 = qs(3);
    qs4 = qs(4);
    qs5 = qs(5);
    qs6 = qs(6);
    
end
function [t2 qs1 qs2 qs3 qs4 qs5 qs6 Deltat2]= SEDRK4c (q0 ,t0 ,Deltat ,Err , fonctiong )
    % Solution ED dq/dt= fonctiong (q)
    % Methode de Runge - Kutta d’ ordre 4
    % Contrôle d’ erreur et interpolation Richardson
    % qs : vecteur final (v,r)
    % q0 : vecteur initial (v,r)
    % Deltat : intervalle de temps
    % Err : Precision requise ( erreur relative )
    % fonctiong : membre de droite de ED.
    % maxit : nombre maximum itérations
    maxit =100;
    ndiv =2;
    nbsubi =1;
    count =1;
    nbel= length (q0 );
    % Solution sur intervalle complet
    qs= SEDRK4t0 (q0 ,t0 ,Deltat , fonctiong );
    % Debut itération
    while (nbsubi < maxit)
        nbsubi = nbsubi * ndiv;
        q2=q0;
        t2=t0;
        Deltat2 = Deltat / nbsubi ;
        % Solution sur pas plus fin
        for ietape =1: nbsubi
            count = count +1;
            q2= SEDRK4t0 (q2 ,t2 , Deltat2 , fonctiong );
            t2=t2+ Deltat2;
        end
        avgsol =( q2 (2: nbel )+ qs (2: nbel ))/2.;
        ErrSol =( q2 (2: nbel )-qs (2: nbel ));
        % Evaluation erreur et la solution moyenne
        MaxErr =max(abs ( ErrSol ./ avgsol ));
        qs=q2;
        if MaxErr < Err
            % Extrapolation de Richardson
            qs (2: nbel )= qs (2: nbel )+ ErrSol /15.;
            break ;
        end
        if nbsubi >= maxit
            fprintf('abort');
            break ;
        end
    end
    qs1 = qs(1);
    qs2 = qs(2);
    qs3 = qs(3);
    qs4 = qs(4);
    qs5 = qs(5);
    qs6 = qs(6);
    
end

function qs= SEDRK4t0 (q0 ,t0 ,Deltat , fonctiong )
%
% Solution ED dq/dt= fonctiong (q)
% Methode de Runge - Kutta d’ ordre 4
% qs : vecteur final [q(tf )]
% q0 : vecteur initial [q(ti )]
% Deltat : intervalle de temps
% fonctiong : membre de droite de ED.
% Ceci est un m- file de matlab
% qui retourne [dq/dt(ti )]
%
    k1= fonctiong (q0 )* Deltat ;
    k2= fonctiong (q0+k1 /2)* Deltat ;
    k3= fonctiong (q0+k2 /2)* Deltat ;
    k4=  fonctiong (q0+k3)* Deltat ;
    qs=q0 +( k1 +2* k2 +2* k3+k4 )/6;
end

function y = fonctionGballe(q0)
    ax = -FrottementAir()*AireBalle()*q0(1);
    ay = -FrottementAir()*AireBalle()*q0(2);
    az = -FrottementAir()*AireBalle()*q0(3) - 9.8;
    y = [ax ay az q0(1) q0(2) q0(3)];
end

function y = fonctionGboite(q0)
    ax = -FrottementAir()*AireBoite()*q0(1);
    ay = -FrottementAir()*AireBoite()*q0(2);
    az = -FrottementAir()*AireBoite()*q0(3) - 9.8;
    y = [ax ay az q0(1) q0(2) q0(3)];
end


function qs= SEDEuler (q0 ,Deltat , fonctiong )
    % Solution ED dq/dt= fonctiong (q,t)
    % Methode de Euler
    % qs : vecteur final [tf q(tf )]
    % q0 : vecteur initial [ti q(ti )]
    % Deltat : intervalle de temps
    % fonctiong : membre de droite de ED.
    % Ceci est un m- file de matlab
    % qui retourne [1 dq/dt(ti )]
    qs=q0+ fonctiong* Deltat + [Deltat 0 0 0 0 0 0];
    maxIt = 100;
    nbIt = 1;
    % Extrapolation de Richardson
    while (nbIt < maxIt)
        % Calculer avec plus de precision
        nbIt = nbIt*2;
        DeltatPetit = Deltat/nbIt;
        q2=q0+ fonctiong* DeltatPetit + [DeltatPetit 0 0 0 0 0 0];
        % Determiner le taux d'erreur
        avgsol=(q2+qs)/2.;
        ErrSol=(q2-qs);
        MaxErr=max(abs(ErrSol./avgsol));
        qs = q2;
        % Si le taux d'erreur est acceptable, quitter la boucle
        if(MaxErr < 0.05)
            qs = qs + ErrSol/15.;
            break;
        end
    end
end

function y = fonctionG(q0, t0)
    norme = norm([q0(2) q0(3) q0(4)]);
    ax = -(((pi*db()*db())/8)*pAir()*Cv()*norme*q0(2))/mb();
    ay = -(((pi*db()*db())/8)*pAir()*Cv()*norme*q0(3))/mb();
    az = -(((pi*db()*db())/8)*pAir()*Cv()*norme*q0(4))/mb() - 9.8;
    y = [0 ax ay az q0(2) q0(3) q0(4)];
end

