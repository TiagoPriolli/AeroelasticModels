function [x,z,xgrid,zgrid,xc,zc,xcdot,zcdot,ni,taui,tauj,Lalpha] = UpdateUVLM_FlatPlate_Test
x = zeros(M+1,N+1);
z = zeros(M+1,N+1);
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
    elseif kk == N+1 % Last points on last strip (extrapolated)
        vecphys1 = Phi(:,:,kk-3)*eta;
        vecphys2 = Phi(:,:,kk-1)*eta;
        vecphys = 1.5.*vecphys2 - 0.5*vecphys1;
    else % Middle points with average between two strips
        vecphys = (Phi(:,:,kk)*eta + Phi(:,:,kk-1)*eta)./2;
    end
    
    hstrip = -vecphys(1);
    alphastrip = vecphys(2);
                
%     z(:,kk) = hstrip - (pvt*c - x(:,kk))*sin(alphastrip);
%     zgrid(:,kk) = hstrip - (pvt*c - xgrid(:,kk))*sin(alphastrip);
    z(:,kk) =  hstrip + (xOri(:,kk) - pvt*c)*sin(alphastrip);
    zgrid(:,kk) =  hstrip + (xgridOri(:,kk) - pvt*c)*sin(alphastrip);
    
    x(:,kk) = xOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    xgrid(:,kk) = xgridOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
 
    
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
    zc(cini:cend,1) = hstrip + (xcOri(cini:cend,1) - pvt*c).*sin(alphastrip);
    xc(cini:cend,1) = xcOri(cini:cend,1).*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    
    zcdot(cini:cend,1) = hdotstrip + (xcOri(cini:cend,1) - pvt*c).*...
        cos(alphastrip).*alphadotstrip;
    xcdot(cini:cend,1) = -xcOri(cini:cend,1).*sin(alphastrip).*alphadotstrip + pvt*c*sin(alphastrip)*alphadotstrip;                              
    
    Lalpha(kk) = -alphastrip;
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