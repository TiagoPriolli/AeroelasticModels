function [uRing,vRing,wRing] = voringMatrix_FlatPlate(xyz,M,N,circulation,...
                                                    xgrid,ygrid,zgrid)

% C Code Bulshit
x3D = xyz(1).*ones(M+1,N+1);
y3D = xyz(2).*ones(M+1,N+1);
z3D = xyz(3).*ones(M+1,N+1);

%% Function Parameters
Eps = 1e-10;

% r1v e r2v
xDist = x3D - xgrid;
yDist = y3D - ygrid; 
zDist = z3D - zgrid;
% r1d e r2d
DistSize = sqrt(xDist.^2 + yDist.^2 + zDist.^2);

ZerosVline = zeros(M,N+1);
ZerosHline = zeros(M+1,N);

%% Vertical Lines (10.115 from Katz/Plotkin)

% r1r2 (M+1xN+1)
xVline = yDist(1:M,:).*zDist(2:M+1,:) - zDist(1:M,:).*yDist(2:M+1,:);
yVline = zDist(1:M,:).*xDist(2:M+1,:) - xDist(1:M,:).*zDist(2:M+1,:);
zVline = xDist(1:M,:).*yDist(2:M+1,:) - yDist(1:M,:).*xDist(2:M+1,:);
% modr1r2
modVline = (xVline.^2) + (yVline.^2) + (zVline.^2);

% r0
xr0v = xgrid(2:M+1,:) - xgrid(1:M,:);
yr0v = ygrid(2:M+1,:) - ygrid(1:M,:);
zr0v = zgrid(2:M+1,:) - zgrid(1:M,:);

r0r1v = xr0v.*xDist(1:M,:) + yr0v.*yDist(1:M,:) + zr0v.*zDist(1:M,:);
r0r2v = xr0v.*xDist(2:M+1,:) + yr0v.*yDist(2:M+1,:) + zr0v.*zDist(2:M+1,:);

% K
K = 1./(4.*pi.*modVline).*(r0r1v./DistSize(1:M,:) - r0r2v./DistSize(2:M+1,:));
uVline = K.*xVline;
vVline = K.*yVline;
wVline = K.*zVline;

% If statement
VTest = (DistSize(1:M,:) < Eps) | (DistSize(2:M+1,:) < Eps) | (modVline < Eps);
uVline(VTest) = ZerosVline(VTest);
vVline(VTest) = ZerosVline(VTest);
wVline(VTest) = ZerosVline(VTest);

%% Horizontal Lines (10.115 from Katz/Plotkin)

% r1r2 (M+1xN+1)
xHline = yDist(:,1:N).*zDist(:,2:N+1) - zDist(:,1:N).*yDist(:,2:N+1);
yHline = zDist(:,1:N).*xDist(:,2:N+1) - xDist(:,1:N).*zDist(:,2:N+1);
zHline = xDist(:,1:N).*yDist(:,2:N+1) - yDist(:,1:N).*xDist(:,2:N+1);
% modr1r2
modHline = (xHline.^2) + (yHline.^2) + (zHline.^2);

% r0
xr0h = xgrid(:,2:N+1) - xgrid(:,1:N);
yr0h = ygrid(:,2:N+1) - ygrid(:,1:N);
zr0h = zgrid(:,2:N+1) - zgrid(:,1:N);

r0r1h = xr0h.*xDist(:,1:N) + yr0h.*yDist(:,1:N) + zr0h.*zDist(:,1:N);
r0r2h = xr0h.*xDist(:,2:N+1) + yr0h.*yDist(:,2:N+1) + zr0h.*zDist(:,2:N+1);

% K
K = 1./(4.*pi.*modHline).*(r0r1h./DistSize(:,1:N) - r0r2h./DistSize(:,2:N+1));
uHline = K.*xHline;
vHline = K.*yHline;
wHline = K.*zHline;

% If statement
HTest = (DistSize(:,1:N) < Eps) | (DistSize(:,2:N+1) < Eps) | (modHline < Eps);
uHline(HTest) = ZerosHline(HTest);
vHline(HTest) = ZerosHline(HTest);
wHline(HTest) = ZerosHline(HTest);

%% Rings

uRing = (+ uHline(1:end-1,:) - uHline(2:end,:)...
         - uVline(:,1:end-1) + uVline(:,2:end)).*circulation;
vRing = (+ vHline(1:end-1,:) - vHline(2:end,:)...
         - vVline(:,1:end-1) + vVline(:,2:end)).*circulation;
wRing = (+ wHline(1:end-1,:) - wHline(2:end,:)...
         - wVline(:,1:end-1) + wVline(:,2:end)).*circulation;

end