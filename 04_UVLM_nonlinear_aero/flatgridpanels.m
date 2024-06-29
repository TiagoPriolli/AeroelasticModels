function [x0,y0,z0,xc,yc,zc,Sp,cp,bp,ni,x,y,z] = flatgridpanels(c,b,M,N)

%% Panel Edges
Iedges = linspace(0,c,M+1);
Jedges = linspace(0,b,N+1);
[y,x] = meshgrid(Jedges,Iedges);
z = zeros(M+1,N+1);

%% Panel Area
S = b*c;            % Wing area
Sp = S/(N*M);       % Panels area
bp = b/N;           % Panels span
cp = c/M;           % Panels chord

%% Vortex Ring Vertice Position
x0 = cp/4.*ones(M,N+1) + x(1:M,:);
y0 = zeros(M,N+1) + y(1:M,:);
z0 = zeros(M,N+1) + z(1:M,:);

%% Collocation Points
xc = ((3*cp/4) .* ones(M,N)) + x(1:M,1:N);
yc = (1/2*bp .* ones(M,N)) + y(1:M,1:N);
zc = zeros(M,N) + z(1:M,1:N);

%% Normal and Tangent Vectors at Collocation Points

ni = UBuildNormalVec(M,N,x,y,z);
end