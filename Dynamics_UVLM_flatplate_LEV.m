% linear dynamics function
% X = [eta; etad]
function [dXdt] = Dynamics_UVLM_flatplate_LEV(t,X,tag,nm,ns,mmod,Uref,Phi,vec_b,vec_dy,Kinematics)
global fix_p_uvlm sch_p_uvlm clvec cmvec a0vec Qvec

M = fix_p_uvlm.M;
N = fix_p_uvlm.N;
c = fix_p_uvlm.c;
%b = fix_p_uvlm.b;
%S = fix_p_uvlm.S;
pvt = fix_p_uvlm.pvt;
%xc = fix_p_uvlm.xc;
yc = fix_p_uvlm.yc;
%zc varies every step
%x = fix_p_uvlm.x;
y = fix_p_uvlm.y;
%z varies every step
%Sp = fix_p_uvlm.Sp;
%cp = fix_p_uvlm.cp;
%bp = fix_p_uvlm.bp;
%xgrid = fix_p_uvlm.xgrid;
ygrid = fix_p_uvlm.ygrid;
%zgrid varies every step
rho = fix_p_uvlm.rho;
%dt = fix_p_uvlm.dt;
%expect_vort = fix_p_uvlm.expect_vort;

eta = X(1:nm,1);
etad = X(nm+1:2*nm,1);

if tag == 1
    Qvec = zeros(nm,1); % vector of generalised forces
    pdyn = 0.5*rho*Uref^2; % dynamic pressure
        
    [x,z,xdot,zdot,xgrid,zgrid,xc,zc,xcdot,zcdot,ni,taui,tauj,Lalpha] = UpdateUVLM_FlatPlate;
    sch_p_uvlm.Lalpha = Lalpha;
    
    %% Run UVLM
    Kinematics.X0 = -Uref*t;
    [coeff,sch_p_uvlm] = UVLM_FlatPlateLEV_coupled(x,y,z,xdot,zdot,xc,yc,zc,xcdot,zcdot,...
                                           xgrid,ygrid,zgrid,ni,taui,tauj,...
                                           sch_p_uvlm,fix_p_uvlm,Kinematics,false,false);
    %[coeff,sch_p_uvlm] = UVLM_FlatPlateLEV(x,y,z,xdot,zdot,xc,yc,zc,xcdot,zcdot,...
    %                                       xgrid,ygrid,zgrid,ni,taui,tauj,...
    %                                       sch_p_uvlm,fix_p_uvlm,Kinematics,false,false);
    

    
    %% Calculate Qvec and etadd
    cl = coeff(:,1);
    cm = coeff(:,3);
    a0 = coeff(:,4);
    for jj = 1:ns
        Qvec = Qvec + Phi(:,:,jj)'*pdyn*[cl(jj)*vec_dy(jj)*(2*vec_b(jj));...
                                         cm(jj)*vec_dy(jj)*(2*vec_b(jj))^2;
                                         0];
    end
    
    clvec = cl;
    cmvec = cm;
    a0vec = a0;
end

if t == 0, Qvec = 0*Qvec; end

etadd = (mmod.m_eta_s_etadd)\...
    ((- mmod.m_eta_s_etad)*etad + ...
    (- mmod.m_eta_s_eta)*eta + ...
    Qvec);


dXdt(1:nm,1) = etad;
dXdt(nm+1:2*nm,1) = etadd;

%%
function [x,z,xdot,zdot,xgrid,zgrid,xc,zc,xcdot,zcdot,ni,taui,tauj,Lalpha] = UpdateUVLM_FlatPlate
x = zeros(M+1,N+1);
xdot = zeros(M+1,N+1);
z = zeros(M+1,N+1);
zdot = zeros(M+1,N+1);
xgrid = zeros(M+1,N+1);
zgrid = zeros(M+1,N+1);
xc = zeros(M*N,1);
zc = zeros(M*N,1);
xcdot = zeros(M*N,1);
zcdot = zeros(M*N,1);
Lalpha = zeros(N,1);
xOri = fix_p_uvlm.x;
xgridOri = fix_p_uvlm.xgrid;
xcOri = fix_p_uvlm.xc;

%% Update each strip
% Update grid
for kk = 1:N+1
    if kk == 1 % First points on the first strip
        vecphys = Phi(:,:,kk)*eta;
        vecphys_dot = Phi(:,:,kk)*etad;
    elseif kk == N+1 % Last points on last strip (extrapolated)
        vecphys1 = Phi(:,:,kk-3)*eta;
        vecphys2 = Phi(:,:,kk-1)*eta;
        vecphys = 1.5.*vecphys2 - 0.5*vecphys1;
        vecphys_dot = 1.5*Phi(:,:,kk-1)*etad - 0.5*Phi(:,:,kk-3)*etad;
    else % Middle points with average between two strips
        vecphys = (Phi(:,:,kk)*eta + Phi(:,:,kk-1)*eta)./2;
        vecphys_dot = (Phi(:,:,kk)*etad + Phi(:,:,kk-1)*etad)./2;
    end
    
    hstrip = -vecphys(1);
    alphastrip = vecphys(2);
    hdotstrip = -vecphys_dot(1);
    alphadotstrip = vecphys_dot(2);
                
    %z(:,kk) = hstrip - (pvt*c - x(:,kk))*sin(alphastrip);
    %zgrid(:,kk) = hstrip - (pvt*c - xgrid(:,kk))*sin(alphastrip);
    z(:,kk) =  (hstrip + (xOri(:,kk) - pvt*c)*sin(alphastrip));
    zgrid(:,kk) =  (hstrip + (xgridOri(:,kk) - pvt*c)*sin(alphastrip));
    zdot(:,kk) = (hdotstrip + (xOri(:,kk) - pvt*c).*...
        cos(alphastrip).*alphadotstrip);
    
    %x(:,kk) = -(xOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c);
    %xgrid(:,kk) = -(xgridOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c);
    x(:,kk) = pvt*c + (xOri(:,kk) - pvt*c)*cos(alphastrip);
    xgrid(:,kk) = pvt*c + (xgridOri(:,kk) - pvt*c)*cos(alphastrip);
    xdot(:,kk) = -(xOri(:,kk) - pvt*c)*sin(alphastrip)*alphadotstrip;
end

% Update Control Points
for kk = 1:N
    vecphys = Phi(:,:,kk)*eta;
    vecphys_dot = Phi(:,:,kk)*etad;
    
    hstrip = -vecphys(1);
    alphastrip = vecphys(2);
    hdotstrip = -vecphys_dot(1);
    alphadotstrip = vecphys_dot(2);
    
    cini = M*(kk-1)+1;
    cend = M*kk;
    
%     zc(cini:cend,1) = hstrip - (pvt*c - xc(cini:cend,1)).*sin(alphastrip);
%     zcdot(cini:cend,1) = hdotstrip - (pvt*c - xc(cini:cend,1)).*...
%                                       cos(alphastrip).*alphadotstrip;
    zc(cini:cend,1) = (hstrip + (xcOri(cini:cend,1) - pvt*c).*sin(alphastrip));
    %xc(cini:cend,1) = xcOri(cini:cend,1).*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    xc(cini:cend,1) = pvt*c + (xcOri(cini:cend,1) - pvt*c)*cos(alphastrip);
    
    zcdot(cini:cend,1) = (hdotstrip + (xcOri(cini:cend,1) - pvt*c).*...
        cos(alphastrip).*alphadotstrip);
    %xcdot(cini:cend,1) = -xcOri(cini:cend,1).*sin(alphastrip).*alphadotstrip + pvt*c*sin(alphastrip)*alphadotstrip;                              
    xcdot(cini:cend,1) = -(xcOri(cini:cend,1) - pvt*c)*sin(alphastrip)*alphadotstrip;
    
    
    Lalpha(kk) = alphastrip;
end

% Update Normal Vectors
ni = zeros(3,M*N);
taui = zeros(3,M*N);
tauj = zeros(3,M*N);
for l = 1:M*N
    [i,j] = ind2sub([M N],l);
    
    A = [xgrid(i+1,j+1);ygrid(i+1,j+1);zgrid(i+1,j+1)] -...
        [xgrid(i,j);ygrid(i,j);zgrid(i,j)];
    
    B = [xgrid(i,j+1);ygrid(i,j+1);zgrid(i,j+1)] -...
        [xgrid(i+1,j);ygrid(i+1,j);zgrid(i+1,j)];
    
    ni(:,l) = (cross(A,B)/norm(cross(A,B)))';
    taui(:,l) = ((A-B)./norm(A-B))';
    tauj(:,l) = ((A+B)./norm(A+B))';
end
end

end


