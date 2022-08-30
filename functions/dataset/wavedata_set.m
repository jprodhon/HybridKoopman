function sim_sets =  wavedata_set(N_SIM,parameter)
    global L NX NT
    N = max(size(parameter));
    lx = [0:L/NX:L]'; lx = lx(1:NX);
    sim_sets = zeros(1+N,N_SIM,NX,NT);
    for k=1:N_SIM
        larg = 0.2+rand/10; ampl = 0.9+rand/10; deb = 0;
        UVp  = 20*ampl * sin((lx-deb)*pi/larg).^2.*(lx<larg+deb).*(lx>deb);
        UV = 0*UVp;
        sim_sets(1,k,:,:) = dynamique_wave(UV,UVp,1);
        for i=1:N
            sim_sets(1+i,k,:,:) = dynamique_wave(UV,UVp,parameter(i));
        end
    end
end

function sol = dynamique_wave(UV,UVp,parameter)
    global T NT L NX
    kappa = 1; 
    dt = T/(NT-1); 
    dx = L/(NX-1); 
    MM  = kappa * speye(NX);
    KK = -gallery('tridiag',ones(NX-1,1),-2*ones(NX,1),ones(NX-1,1))/dx^2; 
    FF = zeros(NX,1);
    sol_esp_temps1 = UV; 
    for iter = 1:NT-1
        [UV,UVp]=Newmark(parameter*KK,0*KK,MM ,UV, UVp,dt,FF,[1, 0; NX, 0]);
        sol_esp_temps1 = [sol_esp_temps1, UV];
    end
    sol = sol_esp_temps1;
end


function [tempd,tempv]=Newmark(K,C,M,d0,v0,dt,f,pdisp)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for dynamic solution of second-order
%  FE equations considering boundary conditions.
% INPUT:
%    K : matrice rigidité, dim(K)= nd x nd
%    C : amortissemeNT visqueux,  dim(C)= nd x nd
%    M : matrice de masse, dim(M)= nd x nd
%    d0 : déplacemeNT initial d(0), dim(f)= nd x 1
%    v0 : vitesse initiale v(0), dim(f)= nd x 1
%    ip : [dt tottime alfa delta [nsnap nhist t(i) ...  dof(i) ... ]] :
%         iNTegration and output parameters
%    f : load vector, dim(f)= n x nstep+1,
%        If dim(f)= n x 1 it is supposed that the values
%        are kept constaNT during time iNTegration.
%    pdisp : matrice des CDL, dim(pdisp)= nbc x nstep+2,
%            where nbc = numbre de conditions aux limites
%            (constaNT ou fonction du temps)
%            la 1ère  colonne coNTieNT les DDL où soNT imposés les déplacemeNTs
%             les colonnes suivaNTes coNTienneNT l'histoire du chargemeNT
% %            If dim(pdisp)= nbc x 2  cela veut dire que le chargemeNT est
%            maiNTenu  constaNT pendaNT l'iNTégration
% OUTPUT:
%    D :   solution en déplacemeNT, coNTenaNT deux vecteurs [d0 d]
%        nhist est les DDL à sortir ici tous
%         dim(D)=nhist x 2
%    V :   vitesse = dérivée du déplacemeNT [v0 v]
%          dim(V)=nhist x 2
%    A :   Accélération [A0 A]. dim(A)=nhist x 2
%-------------------------------------------------------------
    alfa=0.25;    
    delta=0.50;
    b1 = dt*dt*0.5*(1-2*alfa);  
    b2 = (1-delta)*dt;
    b3 = delta*dt;              
    b4 = alfa*dt*dt;
    nd = size(K,1);
    nr = size(pdisp,1);
    a0=M\(f-C*v0-K*d0);
    tempd=zeros(nd,1);    
    tempv=zeros(nd,1);   
    fdof=[1:nd]';
    pdof = pdisp(:,1); 
    pd   = pdisp(:,2);
    pv   = (pd - d0(pdof))/dt;
    fdof(pdof)=[];
    Keff = M(fdof,fdof)+b3*C(fdof,fdof)+b4*K(fdof,fdof);
    dold=d0(fdof);    vold=v0(fdof);    aold=a0(fdof);
    dpred = dold + dt*vold + b1*aold;
    vpred = vold + b2*aold;
    pdeff = C(fdof,pdof)*pv + K(fdof,pdof)*pd;
    reff = f(fdof)-C(fdof,fdof)*vpred-K(fdof,fdof)*dpred-pdeff;
    anew = Keff\reff;  
    dnew = dpred+b4*anew;  
    vnew = vpred+b3*anew;
    tempd(pdof)=pd;  
    tempv(pdof)=pv;
    tempd(fdof)=dnew;         
    tempv(fdof)=vnew;         
end