% linear dynamics function
% X = [eta; etad]
function [dXdt] = Dynamics_UVLM_flatplate(t,X,tag,q,nm,ns,mmod,U,rho,Phi,vec_b,vec_dy)
global fix_p_uvlm sch_p_uvlm clvec cmvec a0vec Qvec

M = fix_p_uvlm.M;
N = fix_p_uvlm.N;
c = fix_p_uvlm.c;
%b = fix_p_uvlm.b;
%S = fix_p_uvlm.S;
pvt = fix_p_uvlm.pvt;
xc = fix_p_uvlm.xc;
yc = fix_p_uvlm.yc;
%zc varies every step
x = fix_p_uvlm.x;
y = fix_p_uvlm.y;
%z varies every step
%Sp = fix_p_uvlm.Sp;
%cp = fix_p_uvlm.cp;
%bp = fix_p_uvlm.bp;
xgrid = fix_p_uvlm.xgrid;
ygrid = fix_p_uvlm.ygrid;
%zgrid varies every step
rho = fix_p_uvlm.rho;
%dt = fix_p_uvlm.dt;
%expect_vort = fix_p_uvlm.expect_vort;

eta = X(1:nm,1);
etad = X(nm+1:2*nm,1);

if tag == 1
    Qvec = zeros(nm,1); % vector of generalised forces
    pdyn = 0.5*rho*U^2; % dynamic pressure
        
    [z,zgrid,zc,zcdot,ni,taui,tauj] = UpdateUVLM_FlatPlate;
    
    %% Run UVLM
    [coeff,sch_p_uvlm] = UVLM_FlatPlate(x,y,z,xc,yc,zc,zcdot,...
                           xgrid,ygrid,zgrid,ni,taui,tauj,...
                           t,q,sch_p_uvlm,fix_p_uvlm);
    
    %% Calculate Qvec and etadd
    cl = coeff(:,1);
    cm = coeff(:,3);
    a0 = coeff(:,4);
    for jj = 1:ns
        Qvec = Qvec + Phi(:,:,jj)'*pdyn*[-cl(jj)*vec_dy(jj)*(2*vec_b(jj));...
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
function [z,zgrid,zc,zcdot,ni,taui,tauj] = UpdateUVLM_FlatPlate
z = zeros(M+1,N+1);
zgrid = zeros(M+1,N+1);
zc = zeros(M*N,1);
zcdot = zeros(M*N,1);
%% Update each strip
% Update grid
for kk = 1:N+1
    if kk == 1 % First points on the first strip
        vecphys = Phi(:,:,kk)*eta;
    elseif kk == N+1 % Last points on last strip (extrapolated)
        vecphys1 = Phi(:,:,kk-3)*eta;
        vecphys2 = Phi(:,:,kk-1)*eta;
        vecphys = 1.5.*vecphys2 - 0.5*vecphys1;
    else % Middle points with average between two strips
        vecphys = (Phi(:,:,kk)*eta + Phi(:,:,kk-1)*eta)./2;
    end
    
    hstrip = -vecphys(1);
    alphastrip = vecphys(2);
                
    z(:,kk) = hstrip - (pvt*c - x(:,kk))*sin(alphastrip);
    zgrid(:,kk) = hstrip - (pvt*c - xgrid(:,kk))*sin(alphastrip);
    
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
    
    zc(cini:cend,1) = hstrip - (pvt*c - xc(cini:cend,1)).*sin(alphastrip);
    zcdot(cini:cend,1) = hdotstrip - (pvt*c - xc(cini:cend,1)).*...
                                      cos(alphastrip).*alphadotstrip;
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


