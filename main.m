function main
    hold on
    fprintf('start\n');
    sol = zeros(8);
    
    % option 1
%     sol = Devoir3([0 0 0], [6.85, 0.0, 6.85], 0.66);
%     celldisp(sol)
%     pointsBalle = sol{4};
%     
%     
%     x1 = pointsBalle(:, 1);
%     y1 = pointsBalle(:, 2);
%     z1 = pointsBalle(:, 3);
%     scatter3(x1,y1,z1);
% 
%     pointsBoite = sol{5};
%     x2 = pointsBoite(:, 1);
%     y2 = pointsBoite(:, 2);
%     z2 = pointsBoite(:, 3);
%     scatter3(x2,y2,z2);
    

   
   % option 1 end

   % option 2
%     sol = Devoir3([0 2.3 0], [6.85, 0.0, 6.85], 0.66);
%     celldisp(sol)
%     pointsBalle = sol{4};
%     
%     
%     x1 = pointsBalle(:, 1);
%     y1 = pointsBalle(:, 2);
%     z1 = pointsBalle(:, 3);
%     scatter3(x1,y1,z1);
% 
%     pointsBoite = sol{5};
%     x2 = pointsBoite(:, 1);
%     y2 = pointsBoite(:, 2);
%     z2 = pointsBoite(:, 3);
%     scatter3(x2,y2,z2);
   % option 2 end

   % option 3
    sol = Devoir3([0 0 0], [28, 0.5, 10], 1.1);
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
    y = RayonMinBoite()*RayonMinBoite() + HauteurBoite()*HauteurBoite()
end

function y = PosRelativeCoinBoite()
    r = RayonMinBoite()
    h = HauteurBoite()
    y = [0 r h/2;
        r/sqrt(2) r/sqrt(2) h/2;
        r 0 h/2;
        r/sqrt(2) -r/sqrt(2) h/2;
        0 -r h/2;
        -r/sqrt(2) -r/sqrt(2) h/2;
        -r 0 h/2;
        -r/sqrt(2) r/sqrt(2) h/2;
        0 r -h/2;
        r/sqrt(2) r/sqrt(2) -h/2;
        r 0 -h/2;
        r/sqrt(2) -r/sqrt(2) -h/2;
        0 -r -h/2;
        -r/sqrt(2) r/sqrt(2) -h/2;
        -r 0 -h/2;
        -r/sqrt(2) -r/sqrt(2) -h/2];
end

function y = Pos0Boite()
    y = [3 0 10];
end

function y = RayonMinBoite()
    y = 0.05;
end

function y = RayonMaxBoite()
    y = sqrt(HauteurBoite()*HauteurBoite()/4 + RayonMinBoite()*RayonMinBoite());
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
    y = 0.5;
end

function y = TauxErreur()
    y = 0.1;
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
    while( i < 100 && result == 0 )
        if(qSolBoite(i,1) >= tballe)
            % Calculer la balle avec precision et imposer son Deltat a la
            % boite
            qSolBalle(i+1,:) = SEDRK4c(qSolBalle(i,2:7), qSolBalle(i,1),  deltaT, TauxErreur(), @fonctionGballe);
            qSolBoite(i+1,:) = SEDRK4cImprecis(qSolBoite(i,2:7), qSolBoite(i,1), qSolBalle(i+1, 8), TauxErreur(), @fonctionGboite);
        else
            % Calculer boite seulement, la balle ne bouge pas au debut
            qSolBoite(i+1,:) = SEDRK4cImprecis(qSolBoite(i,2:7), qSolBoite(i,1), deltaT, TauxErreur(), @fonctionGboite);
            qSolBalle(i+1,:) = [qSolBoite(i+1,1) vballei(1) vballei(2) vballei(3) posBalle(1) posBalle(2) posBalle(3) deltaT];
        end 
        i = i + 1;
        result = DetectionCollision(qSolBoite(i,5:7), qSolBalle(i,5:7),wboitei, qSolBalle(i, 1));
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
    tc = ti(size(ti));
    qfinal = {result-1 vbaf vbof rba rbo tc};
    y = qfinal;
%     qfinal(2) = {vbaf};
%     qfinal(3) = {vbof};
%     qfinal(4) = {rba};
%     qfinal(5) = {rbo};
%     qfinal(6) = {tc};
%     
%     qfinal = [result-1 vbaf vbof rba rbo tc];

end

function y = DetectionCollision(posBoite, posBalle, wboitei, temps)
    %rayon balle (Constante)
    %rayon boite max (distance coin j et centre de masse) (Constante)
    %position balle (posBalle)
    %position boite (posBoite)
    
    if(norm(posBoite-posBalle) > RayonBalle() + RayonMaxBoite())
        if(posBalle(3) - RayonBalle() <= 0)
            y = 1;
        else
            y = 0;
        end
        %0
    else
        y = DetectionCollisionPlansDivision(posBoite, posBalle, wboitei, temps);  %(retourne 0 ou 2)
    end
    
    %0 = pas de collision
    %1 = collision sol
    %2 = collision balle-boite
end

function y = DetectionCollisionPlansDivision(posBoite, posBalle, wboitei, temps)
    %position des coins de la boite.
    %plan de la boite : [i + 8, i + 1, i]. 8 plans + 2.
    posCoinBoiteRotate = RotaterVecteur(wboitei, temps);
    y=0;
    
    for i = 1:8
       planCourant = [posCoinBoiteRotate(i + 8) posCoinBoiteRotate(i + 1) posCoinBoiteRotate(i)];
       vecteurNormal = cross(planCourant(1),planCourant(2));
       vecteurNormalUnitaire = vecteurNormal/norm(vecteurNormal);
       distance = dot(vecteurNormalUnitaire, (posBalle - (posCoinBoiteRotate(i + 8)+posBoite)));
       if(distance > RayonBalle)
           y = -2;
           break;
       end
    end
        

    %0 = pas de collision
    %2 - collision balle-boite
    y = y + 2;
end

function y = RotaterVecteur(temps,wboitei)
    angle = temps*wboitei;
    quatRot = angle2quat(angle(1), angle(2), angle(3));
    pointNonRotate = PosRelativeCoinBoite();
    
    posCoinBoiteRotate = quatrotate(quatRot, pointNonRotate);
    
    y = posCoinBoiteRotate;
    
    
%     x2 = posCoinBoiteRotate(:, 1);
%     y2 = posCoinBoiteRotate(:, 2);
%     z2 = posCoinBoiteRotate(:, 3);
%     scatter3(x2,y2,z2);
end

function y = CalculerCollision()
    y = 0;
end

function y = SEDRK4cImprecis (q0 ,t0 ,Deltat ,Err , fonctiong )
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
    y = [t2 qs(1) qs(2) qs(3) qs(4) qs(5) qs(6)];
    
end
function y= SEDRK4c (q0 ,t0 ,Deltat ,Err , fonctiong )
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
    y = [t2 qs(1) qs(2) qs(3) qs(4) qs(5) qs(6) Deltat];
    
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
%     k1= fonctiong (q0 )* Deltat ;
%     k2= fonctiong (q0+k1 /2)* Deltat ;
%     k3= fonctiong (q0+k2 /2)* Deltat ;
%     k4=  fonctiong (q0+k3)* Deltat ;
%     qs=q0 +( k1 +2* k2 +2* k3+k4 )/6;
fprintf('q0 avant calcul:\n');
fprintf(mat2str(q0));
qs = q0 + fonctiong (q0 )* Deltat;
fprintf('\nqs apres calcul:\n');
fprintf(mat2str(qs));
end

function y = fonctionGballe(q0)
    ax = -FrottementAir()*AireBalle()*q0(1)/MasseBalle();
    ay = -FrottementAir()*AireBalle()*q0(2)/MasseBalle();
    az = -FrottementAir()*AireBalle()*q0(3)/MasseBalle()  - 9.8;
    y = [ax ay az q0(1) q0(2) q0(3)];
end

function y = fonctionGboite(q0)
    ax = -FrottementAir()*AireBoite()*q0(1)/MasseBoite();
    ay = -FrottementAir()*AireBoite()*q0(2)/MasseBoite();
    az = -FrottementAir()*AireBoite()*q0(3)/MasseBoite()  - 9.8;
    fprintf('fonctionGBoite\n');
    
    y = [ax ay az q0(1) q0(2) q0(3)];
    %printf(mat2str(y));
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


