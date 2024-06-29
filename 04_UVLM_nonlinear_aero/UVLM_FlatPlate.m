function [coeff,localState] = UVLM_FlatPlate(x,y,z,xc,yc,zc,zcdot,...
                                   xgrid,ygrid,zgrid,ni,taui,tauj,...
                                   t,q,localState,Parameters)

Uref = Parameters.Uref;
M = Parameters.M;
N = Parameters.N;
Sp = Parameters.Sp;
cp = Parameters.cp;
bp = Parameters.bp;

%steps = Parameters.expect_vort;
c = Parameters.c;
pvt = Parameters.pvt;
S = Parameters.S;
rho = Parameters.rho;

NVRING = localState.NVRING;
circwake = localState.circwake;
Xwakegrid = localState.Xwakegrid;
Ywakegrid = localState.Ywakegrid;
Zwakegrid = localState.Zwakegrid;
circd = localState.circd;

dt = Parameters.dt*c/Uref;

%% Kinematics
X0 = -Uref*t;
Y0 = 0;
Z0 = 0;

X0dot = -Uref;
Y0dot = 0;
Z0dot = 0;

P = 0;
Q = 0;
R = 0;

% Inertial to Body Transformation
% Phi = 0; Theta = q; Psi = 0;
tbi = [1 0 0;0 1 0;0 0 1]*...
    ([cos(q) 0 -sin(q);0 1 0;sin(q) 0 cos(q)]*...
    [1 0 0;0 1 0;0 0 1]);

% Body to Inertial Transformation
tib = tbi';


%% Place Trailing Edge Vertice
for j = 1:N+1
    XYZ = tib*[x(end,j); y(end,j); z(end,j)]+[X0;Y0;Z0];
    XYZTEV = XYZ - 0.2.*([X0dot;Y0dot;Z0dot].*dt);
    xyztev = tbi*(XYZTEV - [X0;Y0;Z0]);
    xgrid(M+1,j) = xyztev(1);
    ygrid(M+1,j) = xyztev(2);
    zgrid(M+1,j) = xyztev(3);
    Xwakegrid(1,j) = XYZTEV(1);
    Ywakegrid(1,j) = XYZTEV(2);
    Zwakegrid(1,j) = XYZTEV(3);
end

for j = 1:N
    AI = [xgrid(M+1,j+1);ygrid(M+1,j+1);zgrid(M+1,j+1)] -...
        [xgrid(M,j);ygrid(M,j);zgrid(M,j)];
    BI = [xgrid(M,j+1);ygrid(M,j+1);zgrid(M,j+1)] -...
        [xgrid(M+1,j);ygrid(M+1,j);zgrid(M+1,j)];
    ni(:,j*M) = (cross(AI,BI)/norm(cross(AI,BI)))';
    taui(:,j*M) = ((AI-BI)./norm(AI-BI))';
    tauj(:,j*M) = ((AI+BI)./norm(AI+BI))';
end

%% Update Global Grid Location
XYZc = zeros(M*N,3);
Xgrid = zeros(M+1,N+1);
Ygrid = zeros(M+1,N+1);
Zgrid = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        if j<=N && i<=M
            l = sub2ind([M N],i,j);
            XYZc(l,:) = (tib*[xc(l);yc(l);zc(l)]+...
                [X0;Y0;Z0])';
        end
        Ggrid = tib*[xgrid(i,j);ygrid(i,j);zgrid(i,j)]+[X0;Y0;Z0];
        Xgrid(i,j) = Ggrid(1);
        Ygrid(i,j) = Ggrid(2);
        Zgrid(i,j) = Ggrid(3);
    end
end

%% Influence Coefficients
% Normal Velocity is a combination of self-induced velocity, kinematic
% velocity and wake induced velocity
A = InfluenceCoeff_FlatPlate;

%% Establish Right-Hand Side
[RHS,uvw] = RightHandSide_FlatPlate;

%% Solve Linear Equations
circ = A\RHS;
circ = reshape(circ,M,N);

%% Verify LEV Generation

% Calculate LESP Distribution
A0 = CalcLESP(circ(1,:));

%% Secondary Computations
dP = zeros(M,N);
dF = zeros(3,M,N);

dPhidTaui = [circ(1,:); circ(2:M,:)-circ(1:M-1,:)]./cp;
dPhidTauj = [circ(:,1)-circ(:,1), circ(:,2:N)-circ(:,1:N-1)]./bp;

for j = 1:N
    for i = 1:M
        k = sub2ind([M N],i,j);
        
        dgamadt = (circ(i,j)-circd(i,j))./(dt);
        
        dP(i,j) = rho*(dot(uvw(:,k),taui(:,k).*dPhidTaui(i,j)) + ...
            dot(uvw(:,k),tauj(:,k).*dPhidTauj(i,j)) + ...
            dgamadt);
        
        dF(:,i,j) = -(dP(i,j).*Sp).*ni(:,k);
    end
end

dL = permute(dF(3,:,:),[2,3,1]);
dD = permute(dF(2,:,:),[2,3,1]);
dM = reshape(dL(:).*(pvt*c - (xc(:)-cp/2)),M,N);

L = sum(dL,1)';
D = sum(dD,1)';
Mom = sum(dM,1)';

% FlightCoefficients
cL = L./((S/N)*(1/2*rho*Uref^2));
cD = D./((S/N)*(1/2*rho*Uref^2));
cM = Mom./((S/N)*c*(1/2*rho*Uref^2));

coeff = [cL,cD,cM,A0'];

circd = circ;

%% Vortex Wake Rollup

[newXwakegrid,newYwakegrid,newZwakegrid] = UVLMWakeRollupNoSym;

% Update Global Wake Position and Circulation in the Matrices
Xwakegrid(2:NVRING+2,:) = newXwakegrid(1:NVRING+1,:);
Ywakegrid(2:NVRING+2,:) = newYwakegrid(1:NVRING+1,:);
Zwakegrid(2:NVRING+2,:) = newZwakegrid(1:NVRING+1,:);
circwake(2:NVRING+1,:) = circwake(1:NVRING,:);
circwake(1,:) = circ(M,:);
NVRING = NVRING+1;

% Delete Ring Line that is too far
if all(abs(Xwakegrid(NVRING+1,:) - Xgrid(M+1,:)) > c*10)
    Xwakegrid(NVRING+1,:) = zeros(1,N+1);
    Ywakegrid(NVRING+1,:) = zeros(1,N+1);
    Zwakegrid(NVRING+1,:) = zeros(1,N+1);
    circwake(NVRING,:) = zeros(1,N);
    NVRING = NVRING-1;
end

%% Output States
localState.NVRING = NVRING;
localState.circwake = circwake;
localState.Xwakegrid = Xwakegrid;
localState.Ywakegrid = Ywakegrid;
localState.Zwakegrid = Zwakegrid;
localState.circd = circd;
localState.Xgrid = Xgrid;
localState.Ygrid = Ygrid;
localState.Zgrid = Zgrid;

%
%
%
%% Local Functions
    function A = InfluenceCoeff_FlatPlate
        A = zeros(N*M,N*M);
        circulation = ones(M,N);
        
        for kk = 1:(N*M)
            xyz = [xc(kk);yc(kk);zc(kk)];
            xyzi = [xc(kk);-yc(kk);zc(kk)];
            [uRing,vRing,wRing] = voringMatrix_FPbound...
                 (xyz,M,N,circulation,xgrid,ygrid,zgrid);
            [uRingi,vRingi,wRingi] = voringMatrix_FPbound...
                   (xyzi,M,N,circulation,xgrid,ygrid,zgrid);
%             [uRing,vRing,wRing] = voringMatrix_FlatPlateLOlsen...
%                 (xyz,M,N,circulation,xgrid,ygrid,zgrid);
%             [uRingi,vRingi,wRingi] = voringMatrix_FlatPlateLOlsen...
%                 (xyzi,M,N,circulation,xgrid,ygrid,zgrid);
            
            ucoeff = uRing + uRingi;
            vcoeff = vRing - vRingi;
            wcoeff = wRing + wRingi;
            A(kk,:) = ucoeff(:)'.*(ni(1,kk)*ones(1,M*N)) + vcoeff(:)'.*(ni(2,kk)*ones(1,M*N))...
                + wcoeff(:)'.*(ni(3,kk)*ones(1,M*N));
        end
    end

    function [RHS,uvw] = RightHandSide_FlatPlate
        RHS = zeros(N*M,1);
        uvw = zeros(3,N*M);
        for kk = 1:(N*M)
            
            % Wake induced Velocity
            XYZl = [XYZc(kk,1);XYZc(kk,2);XYZc(kk,3)];
            XYZli = [XYZc(kk,1);-XYZc(kk,2);XYZc(kk,3)];
            if NVRING == 0
                uvwrhs = [0;0;0];
            else
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                uvwrhs = tbi*[sum(URing(:)) + sum(URingi(:));...
                              sum(VRing(:)) - sum(VRingi(:));...
                              sum(WRing(:)) + sum(WRingi(:))];
            end
            
            % Kinematics induced Velocity
            uvwt = tbi*[-X0dot;-Y0dot;-Z0dot]...
                + [-Q*zc(kk) + R*yc(kk);-R*xc(kk) + P*zc(kk);-P*yc(kk) + Q*xc(kk)];
            
            % Total induced Velocity
            uvw(:,kk) = uvwrhs + uvwt + [0;0;zcdot(kk)];
            
            % Right-Hand Side
            RHS(kk) = dot(-uvw(:,kk),ni(:,kk));
        end
    end

    function [newXwakegrid,newYwakegrid,newZwakegrid] = UVLMWakeRollupNoSym
        
        newXwakegrid = zeros(NVRING+1,N+1);
        newYwakegrid = zeros(NVRING+1,N+1);
        newZwakegrid = zeros(NVRING+1,N+1);
        
        for kk = 1:(NVRING+1)*(N+1)
            [ii,jj] = ind2sub([(NVRING+1),(N+1)],kk);
            XYZl = [Xwakegrid(ii,jj);Ywakegrid(ii,jj);Zwakegrid(ii,jj)];
            XYZli = [1;-1;1].*XYZl;
            
            % Bound Rings Induced Velocity

            [UBInduced,VBInduced,WBInduced] = voringMatrix_FPbound(XYZl,M,N,circ,...
                                              Xgrid,Ygrid,Zgrid);
            [UBInducedi,VBInducedi,WBInducedi] = voringMatrix_FPbound(XYZli,M,N,circ,...
                                                 Xgrid,Ygrid,Zgrid);
%             [UBInduced,VBInduced,WBInduced] = voringMatrix_FlatPlateLOlsen(XYZl,M,N,circ,...
%                                               Xgrid,Ygrid,Zgrid);
%             [UBInducedi,VBInducedi,WBInducedi] = voringMatrix_FlatPlateLOlsen(XYZli,M,N,circ,...
%                                                  Xgrid,Ygrid,Zgrid);

            BoundUVWInduced = [sum(UBInduced(:)) + sum(UBInducedi(:));...
                               sum(VBInduced(:)) - sum(VBInducedi(:));...
                               sum(WBInduced(:)) + sum(WBInducedi(:))];
                        
            % Wake Rings Induced Velocity
            if NVRING == 0
                WakeUVWInduced = [0;0;0];
            else
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                WakeUVWInduced = [sum(UWInduced(:)) + sum(UWInducedi(:));...
                                  sum(VWInduced(:)) - sum(VWInducedi(:));...
                                  sum(WWInduced(:)) + sum(WWInducedi(:))];
            end
            
            % Wake Roolup
            newXwakegrid(ii,jj) = XYZl(1) + (BoundUVWInduced(1) + WakeUVWInduced(1))*dt;
            newYwakegrid(ii,jj) = XYZl(2) + (BoundUVWInduced(2) + WakeUVWInduced(2))*dt; 
            newZwakegrid(ii,jj) = XYZl(3) + (BoundUVWInduced(3) + WakeUVWInduced(3))*dt; 
        end
    end

    function LESPValue = CalcLESP(circulation)
        LESPValue = 1.13.*circulation./(-X0dot*c*(acos(1-2*cp/c)+sin(acos(1-2*cp/c))));
    end
end
