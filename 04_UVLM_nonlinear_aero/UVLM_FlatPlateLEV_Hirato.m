function [coeff,localState,forces] = UVLM_FlatPlateLEV_Hirato(x,y,z,xdot,zdot,xc,yc,zc,xcdot,zcdot,...
    xgrid,ygrid,zgrid,ni,taui,tauj,...
    localState,Parameters,Kinematics,viz,showVelocity)
if isempty(viz)
    viz=false;
end
if isempty(showVelocity)
    showVelocity=false;
end

%% Extract Data from Structures
Uref = Parameters.Uref;
M = Parameters.M;
N = Parameters.N;
Sp = Parameters.Sp;
cp = Parameters.cp;
bp = Parameters.bp;

%steps = Parameters.expect_vort;
c = Parameters.c;
b = Parameters.b;
pvt = Parameters.pvt;
S = Parameters.S;
rho = Parameters.rho;
LESPcrit = Parameters.LESPcrit;

NVRING = localState.NVRING;
circwake = localState.circwake;
Xwakegrid = localState.Xwakegrid;
Ywakegrid = localState.Ywakegrid;
Zwakegrid = localState.Zwakegrid;
NLEVRingST = localState.NLEVRing;
circwakeLEV = localState.circwakeLEV;
XwakegridLEV = localState.XwakegridLEV;
YwakegridLEV = localState.YwakegridLEV;
ZwakegridLEV = localState.ZwakegridLEV;
circd = localState.circd;

dt = Parameters.dt*c/Uref;

%% Kinematics
X0 = Kinematics.X0;
Y0 = Kinematics.Y0;
Z0 = Kinematics.Z0;

X0dot = Kinematics.X0dot;
Y0dot = Kinematics.Y0dot;
Z0dot = Kinematics.Z0dot;

P = Kinematics.P;
Q = Kinematics.Q;
R = Kinematics.R;

Phi = Kinematics.Phi;
Theta = Kinematics.Theta;
Psi = Kinematics.Psi;

% Inertial to Body Transformation
tbi = [1 0 0;0 cos(Phi) sin(Phi);0 -sin(Phi) cos(Phi)]*...
    ([cos(Theta) 0 -sin(Theta);0 1 0;sin(Theta) 0 cos(Theta)]*...
    [cos(Psi) sin(Psi) 0;-sin(Psi) cos(Psi) 0;0 0 1]);

% Body to Inertial Transformation
tib = tbi';


%% Place Trailing Edge Vertice
for j = 1:N+1
    XYZ = tib*[x(end,j); y(end,j); z(end,j)]+[X0;Y0;Z0];
    XYZTEV = XYZ - 0.3.*([X0dot;Y0dot;Z0dot].*dt);
    xyztev = tbi*(XYZTEV - [X0;Y0;Z0]);
    xgrid(M+1,j) = xyztev(1);
    ygrid(M+1,j) = xyztev(2);
    zgrid(M+1,j) = xyztev(3);
    Xwakegrid(1,j) = XYZTEV(1);
    Ywakegrid(1,j) = XYZTEV(2);
    Zwakegrid(1,j) = XYZTEV(3);
end
% for l = 1:M*N
%     [i,j] = ind2sub([M N],l);
%
%     A = [xgrid(i+1,j+1);ygrid(i+1,j+1);zgrid(i+1,j+1)] -...
%         [xgrid(i,j);ygrid(i,j);zgrid(i,j)];
%
%     B = [xgrid(i,j+1);ygrid(i,j+1);zgrid(i,j+1)] -...
%         [xgrid(i+1,j);ygrid(i+1,j);zgrid(i+1,j)];
%
%     ni(:,l) = (cross(A,B)/norm(cross(A,B)))';
%     taui(:,l) = ((A-B)./norm(A-B))';
%     tauj(:,l) = ((A+B)./norm(A+B))';
% end

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

% See wake
% if viz
%     visualizeWakeDev(1)
% end

%% Influence Coefficients
% Normal Velocity is a combination of self-induced velocity, kinematic
% velocity and wake induced velocity
% [Eq 4.33]
A = InfluenceCoeff_FlatPlate;

%% Establish Right-Hand Side
pseudoVortexflag = zeros(1,N);
[RHS,uvw] = RightHandSide_FlatPlate;

%% Solve Linear Equations
circ = A\RHS; % [Eq. 4.51]
circ = reshape(circ,M,N);

%% Verify LEV Generation

% Calculate LESP Distribution
A0 = CalcLESP(circ(1,:)); % [Eq. 4.52]
levflag = abs(A0) > LESPcrit;
pseudoVortexflag = levflag;

%LEVPosition

if any(levflag)
    
    % LEV is shed on approppriate strips
    ind = find(levflag);
    
        NLEVRingST(ind) = NLEVRingST(ind) + 1;
        indG = unique([ind ind+1]);
    
        % LEV wakegrid is updated
        isFirst = false(1,N+1);
        isFirst(indG) = (XwakegridLEV(2,indG) == 0);
        XwakegridLEV(3:end,indG) = XwakegridLEV(2:end-1,indG);
        YwakegridLEV(3:end,indG) = YwakegridLEV(2:end-1,indG);
        ZwakegridLEV(3:end,indG) = ZwakegridLEV(2:end-1,indG);
        UVW_LEV = zeros(3,N+1);
        for ll = 1:numel(indG)
            XYZ = tib*[x(1,indG(ll)); y(1,indG(ll)); z(1,indG(ll))]+[X0;Y0;Z0];
            XwakegridLEV(1,indG(ll)) = XYZ(1);
            YwakegridLEV(1,indG(ll)) = XYZ(2);
            ZwakegridLEV(1,indG(ll)) = XYZ(3);
            UVW_LEV(:,indG(ll)) = [-X0dot;-Y0dot;-Z0dot]...
                + tib*([-Q*z(1,indG(ll)) + R*y(1,indG(ll));-R*x(1,indG(ll)) + P*z(1,indG(ll));-P*y(1,indG(ll)) + Q*x(1,indG(ll))]...
                + [xdot(1,indG(ll));0;zdot(1,indG(ll))]);
        end
        % Kinematics Induced Velocity
        XwakegridLEV(2,indG) = 2/3.*XwakegridLEV(1,indG) + 1/3.*XwakegridLEV(3,indG);
        YwakegridLEV(2,indG) = 2/3.*YwakegridLEV(1,indG) + 1/3.*YwakegridLEV(3,indG);
        ZwakegridLEV(2,indG) = 2/3.*ZwakegridLEV(1,indG) + 1/3.*ZwakegridLEV(3,indG);
        if any(isFirst)
            XwakegridLEV(2,isFirst) = XwakegridLEV(1,isFirst) + 0.5.*UVW_LEV(1,isFirst).*dt;
            YwakegridLEV(2,isFirst) = YwakegridLEV(1,isFirst) + 0.5.*UVW_LEV(2,isFirst).*dt;
            ZwakegridLEV(2,isFirst) = ZwakegridLEV(1,isFirst) + 0.5.*UVW_LEV(3,isFirst).*dt;
        end
    
        % LEV circulation matrix is updated
        circwakeLEV(2:end,ind) = circwakeLEV(1:end-1,ind);
    
    % Pre-allocate iteration matrices
    GammaLevIter = zeros(100,N);
    LESPIter = zeros(100,N);
    EPS = 10e-6;
    
    % Initial k-2 iteration -> Original Values of A0 obtained for Gamma=0
    LESPIter(1,:) = A0;
    
    % Initial k-1 iteration
    GammaLevIter(2,ind) = -0.01;
    circwakeLEV(1,ind) = GammaLevIter(2,ind);
    [RHS,uvw] = RightHandSide_FlatPlate;
    circ = A\RHS;
    circ = reshape(circ,M,N);
    LESPIter(2,:) = CalcLESP(circ(1,:));
    
    i = 3;
    
    % LEV strengh iteration
    while any(levflag) %[Eq. 4.47]
        GammaLevIter(i,ind) = (GammaLevIter(i-1,ind) - GammaLevIter(i-2,ind))./...
            (LESPIter(i-1,ind) - LESPIter(i-2,ind)).*(LESPcrit - LESPIter(i-1,ind))...
            + GammaLevIter(i-1,ind);
        circwakeLEV(1,ind) = GammaLevIter(i,ind);
        [RHS,uvw] = RightHandSide_FlatPlate;
        circ = A\RHS;
        circ = reshape(circ,M,N);
        
        LESPIter(i,:) = CalcLESP(circ(1,:));
        A0 = LESPIter(i,:);
        levflag(ind) = ~(abs(LESPIter(i,ind)) < (LESPcrit + EPS) &...
            abs(LESPIter(i,ind)) > (LESPcrit - EPS));
        i = i + 1;
        if i == 101
            save('ErrorWorkspace')
            error('2D Iteration failed...')
        end
    end
    circwake(1,pseudoVortexflag) = circwake(1,pseudoVortexflag) + circwakeLEV(1,pseudoVortexflag);
else
    circwakeLEV(1,:) = zeros(1,N);    
end

%% Secondary Computations
dP = zeros(M,N);
df = zeros(3,M,N);
dF = zeros(3,M,N);
ds = zeros(3,M,N);
dS = zeros(3,M,N);
%alpha = zeros(M,N);

%circP = circ;
% [Eq. 4.76]
%circP(:,pseudoVortexflag) = circ(:,pseudoVortexflag) + ones(M,1).*circwakeLEV(1,pseudoVortexflag);

dPhidTaui = [circ(1,:); circ(2:M,:)-circ(1:M-1,:)]./cp;
dPhidTauj = [circ(:,1), circ(:,2:N)-circ(:,1:N-1)]./bp;

for j = 1:N
    for i = 1:M
        k = sub2ind([M N],i,j);
        
        dcirc = circ(i,j)-circd(i,j);
        if any(pseudoVortexflag)
            dcircLev = circwakeLEV(1,j)-circwakeLEV(2,j);
        else
            dcircLev = 0;
        end
        dgamadt = (dcircLev+dcirc)./(dt);
        
        % [Eq.4.71.5]
        dP(i,j) = rho*(dot(uvw(:,k),taui(:,k).*dPhidTaui(i,j)) + ...
            dot(uvw(:,k),tauj(:,k).*dPhidTauj(i,j)) + ...
            dgamadt);
        
        df(:,i,j) = (dP(i,j).*Sp).*ni(:,k); %Force Vector on Local Reference
        dF(:,i,j) = tib*df(:,i,j); %Force Vector on Global Reference
        ds(:,i,j) = (pi*rho*c/M*(-X0dot)^2*A0(j)^2*bp).*[-cos(localState.Lalpha(j));0;sin(localState.Lalpha(j))]; %Suction Force Vector on Local Reference
        dS(:,i,j) = tib*ds(:,i,j); %Suction Force Vector on Global Reference
        
    end
end

dn = permute(df(3,:,:),[2,3,1]) + permute(ds(3,:,:),[2,3,1]); %Total normal force to panel
dM = reshape(dn(:).*(pvt*c - (xc(:)-cp/2)),M,N); %Panel's moment around strip's elastic axis
Mom = sum(dM,1)'; %Total moment around strip's elastic axis

Lift = sum(sum(permute(dF(3,:,:),[2,3,1]),1));
Suction = sum(sum(permute(dS(3,:,:),[2,3,1]),1));
dL = permute(dF(3,:,:),[2,3,1]) + permute(dS(3,:,:),[2,3,1]); %Total Lift force to panel
L = sum(dL,1)'; %Total Lift force to strip

dD = permute(dF(1,:,:),[2,3,1]) + permute(dS(1,:,:),[2,3,1]); %Total F=Drag force to panel
D = sum(dD,1)'; %Total Drag force to strip

% FlightCoefficients
cL = L./((S/N)*(1/2*rho*Uref^2));
cD = D./((S/N)*(1/2*rho*Uref^2));
cM = Mom./((S/N)*c*(1/2*rho*Uref^2));

coeff = [cL,cD,cM,A0'];
forces = [Lift,Suction];
circd = circ;

%% Vortex Wake Rollup

% See wake
if viz
    visualizeWakeDev
end

% Calculate Wake Rollup
[newXwakegrid,newYwakegrid,newZwakegrid] = UVLMWakeRollup;

if any(NLEVRingST > 0)
    
    % Calculate LEV Wake Rollup
    [newXwakegridLEV,newYwakegridLEV,newZwakegridLEV] = UVLMLEVWakeRollup;
    
    % Update LEV Position
    XwakegridLEV(1:max(NLEVRingST)+1,:) = newXwakegridLEV(1:max(NLEVRingST)+1,:);
    YwakegridLEV(1:max(NLEVRingST)+1,:) = newYwakegridLEV(1:max(NLEVRingST)+1,:);
    ZwakegridLEV(1:max(NLEVRingST)+1,:) = newZwakegridLEV(1:max(NLEVRingST)+1,:);
    
    % Delete LEV too far
    %     Lev2Far = abs(XwakegridLEV(max(NLEVRingST)+1,:) - Xgrid(M+1,:)) > c*10 ...
    %         | XwakegridLEV(max(NLEVRingST)+1,:) == 0;
    %     if all(Lev2Far)
    %         %Lev2Delete = Lev2Far & XwakegridLEV(max(NLEVRingST)+1,:) ~= 0;
    %         XwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
    %         YwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
    %         ZwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
    %         circwakeLEV(max(NLEVRingST),:) = zeros(1,N);
    %         Ring2Delete = NLEVRingST == max(NLEVRingST);
    %         NLEVRingST(Ring2Delete) = NLEVRingST(Ring2Delete) - 1;
    %     end
end

% Update Global Wake Position and Circulation in the Matrices
Xwakegrid(2:NVRING+2,:) = newXwakegrid(1:NVRING+1,:);
Ywakegrid(2:NVRING+2,:) = newYwakegrid(1:NVRING+1,:);
Zwakegrid(2:NVRING+2,:) = newZwakegrid(1:NVRING+1,:);
circwake(2:NVRING+1,:) = circwake(1:NVRING,:);
circwake(1,:) = circ(M,:);
if any(pseudoVortexflag)
    % [Eq. 4.45]
    circwake(1,pseudoVortexflag) = circwake(1,pseudoVortexflag) + circwakeLEV(1,pseudoVortexflag);
end
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
localState.NLEVRing = NLEVRingST;
localState.circwakeLEV = circwakeLEV;
localState.XwakegridLEV = XwakegridLEV;
localState.YwakegridLEV = YwakegridLEV;
localState.ZwakegridLEV = ZwakegridLEV;
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
%            xyzi = [xc(kk);-yc(kk);zc(kk)];
            [uRing,vRing,wRing] = voringMatrix_FPbound...
                (xyz,M,N,circulation,xgrid,ygrid,zgrid);
%             [uRingi,vRingi,wRingi] = voringMatrix_FPbound...
%                 (xyzi,M,N,circulation,xgrid,ygrid,zgrid);
            
            ucoeff = uRing;% + uRingi;
            vcoeff = vRing;% - vRingi;
            wcoeff = wRing;% + wRingi;
            A(kk,:) = ucoeff(:)'.*(ni(1,kk)*ones(1,M*N)) + vcoeff(:)'.*(ni(2,kk)*ones(1,M*N))...
                + wcoeff(:)'.*(ni(3,kk)*ones(1,M*N));
        end
    end

    function [RHS,uvw] = RightHandSide_FlatPlate
        RHS = zeros(N*M,1);
        uvw = zeros(3,N*M);
        for kk = 1:(N*M)
            
            % Point of Interest
            XYZl = [XYZc(kk,1);XYZc(kk,2);XYZc(kk,3)];
%             XYZli = [XYZc(kk,1);-XYZc(kk,2);XYZc(kk,3)];
            
            % Wake induced Velocity
            if NVRING == 0
                uvwTEV = [0;0;0];
            else
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
%                 [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
%                     Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                uvwTEV = tbi*[sum(URing(:));% + sum(URingi(:));...
                    sum(VRing(:));% - sum(VRingi(:));...
                    sum(WRing(:))];% + sum(WRingi(:))];
            end
            
            % LEV induced Velocity
            if all(NLEVRingST == 0)
                uvwLEV = [0;0;0];
            else
                Size = max(NLEVRingST);
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,Size,N,circwakeLEV(1:Size,:),...
                    XwakegridLEV(1:Size+1,:),YwakegridLEV(1:Size+1,:),ZwakegridLEV(1:Size+1,:));
%                 [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,Size,N,circwakeLEV(1:Size,:),...
%                     XwakegridLEV(1:Size+1,:),YwakegridLEV(1:Size+1,:),ZwakegridLEV(1:Size+1,:));
                uvwLEV = tbi*[sum(URing(:));% + sum(URingi(:));...
                    sum(VRing(:));% - sum(VRingi(:));...
                    sum(WRing(:))];% + sum(WRingi(:))];
            end
            
            % Kinematics Induced Velocity
            uvwt = tbi*[-X0dot;-Y0dot;-Z0dot]...
                + [-Q*zc(kk) + R*yc(kk);-R*xc(kk) + P*zc(kk);-P*yc(kk) + Q*xc(kk)]...
                + [xcdot(kk);0;zcdot(kk)];
            
            % Unbounded Induced Velocity
            uvw(:,kk) =  uvwt + uvwTEV + uvwLEV;
            
            % Pseudo Vortex Induced Velocity
            if any(pseudoVortexflag)
                XpseudoGrid = [XwakegridLEV(1,:);Xwakegrid(1,:)];
                YpseudoGrid = [YwakegridLEV(1,:);Ywakegrid(1,:)];
                ZpseudoGrid = [ZwakegridLEV(1,:);Zwakegrid(1,:)];
                pseudocirc = zeros(1,N);
                pseudocirc(1,pseudoVortexflag) = -circwakeLEV(1,pseudoVortexflag);
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,1,N,pseudocirc,...
                    XpseudoGrid,YpseudoGrid,ZpseudoGrid);
%                 [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,1,N,pseudocirc,...
%                     XpseudoGrid,YpseudoGrid,ZpseudoGrid);
                uvwPSEUDO = tbi*[sum(URing(:));% + sum(URingi(:));...
                    sum(VRing(:));% - sum(VRingi(:));...
                    sum(WRing(:))];% + sum(WRingi(:))];
            else
                uvwPSEUDO = 0;
            end
            
            % Right-Hand Side
            uvwrhs = uvw(:,kk) + uvwPSEUDO;
            RHS(kk) = dot(-uvwrhs,ni(:,kk));
        end
    end

    function [newXwakegrid,newYwakegrid,newZwakegrid] = UVLMWakeRollup
        
        newXwakegrid = zeros(NVRING+1,N+1);
        newYwakegrid = zeros(NVRING+1,N+1);
        newZwakegrid = zeros(NVRING+1,N+1);
        
        for kk = 1:(NVRING+1)*(N+1)
            [ii,jj] = ind2sub([(NVRING+1),(N+1)],kk);
            XYZl = [Xwakegrid(ii,jj);Ywakegrid(ii,jj);Zwakegrid(ii,jj)];
%             XYZli = [1;-1;1].*XYZl;
            
            % Bound Rings Induced Velocity
            [UBInduced,VBInduced,WBInduced] = voringMatrix_FPbound(XYZl,M,N,circ,...
                Xgrid,Ygrid,Zgrid);
%             [UBInducedi,VBInducedi,WBInducedi] = voringMatrix_FPbound(XYZli,M,N,circ,...
%                 Xgrid,Ygrid,Zgrid);
            BoundUVWInduced = [sum(UBInduced(:));% + sum(UBInducedi(:));...
                sum(VBInduced(:));% - sum(VBInducedi(:));...
                sum(WBInduced(:))];% + sum(WBInducedi(:))];
            
            % Wake Rings Induced Velocity
            if NVRING == 0
                WakeUVWInduced = [0;0;0];
            else
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
%                 [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
%                     Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                WakeUVWInduced = [sum(UWInduced(:));% + sum(UWInducedi(:));...
                    sum(VWInduced(:));% - sum(VWInducedi(:));...
                    sum(WWInduced(:))];% + sum(WWInducedi(:))];
            end
            
            % LEV Wake Rings Induced Velocity
            if all(NLEVRingST == 0)
                WakeUVWInducedLEV = [0;0;0];
            else
                MeshSize = max(NLEVRingST);
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,MeshSize,N,circwakeLEV(1:MeshSize,:),...
                    XwakegridLEV(1:MeshSize+1,:),YwakegridLEV(1:MeshSize+1,:),ZwakegridLEV(1:MeshSize+1,:));
%                 [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,MeshSize,N,circwakeLEV(1:MeshSize,:),...
%                     XwakegridLEV(1:MeshSize+1,:),YwakegridLEV(1:MeshSize+1,:),ZwakegridLEV(1:MeshSize+1,:));
                WakeUVWInducedLEV = [sum(UWInduced(:));% + sum(UWInducedi(:));...
                    sum(VWInduced(:));% - sum(VWInducedi(:));...
                    sum(WWInduced(:))];% + sum(WWInducedi(:))];
            end
            
            
            CurrentVelocityTerm = [(BoundUVWInduced(1) + WakeUVWInduced(1) + WakeUVWInducedLEV(1));
                (BoundUVWInduced(2) + WakeUVWInduced(2) + WakeUVWInducedLEV(2));
                (BoundUVWInduced(3) + WakeUVWInduced(3) + WakeUVWInducedLEV(3))];
            
            newXwakegrid(ii,jj) = XYZl(1) + (CurrentVelocityTerm(1))*dt;
            newYwakegrid(ii,jj) = XYZl(2) + (CurrentVelocityTerm(2))*dt;
            newZwakegrid(ii,jj) = XYZl(3) + (CurrentVelocityTerm(3))*dt;
            
        end
    end

    function [newXwakegridLEV,newYwakegridLEV,newZwakegridLEV] = UVLMLEVWakeRollup
        
        GridSize = max(NLEVRingST);
        newXwakegridLEV = zeros(GridSize+1,N+1);
        newYwakegridLEV = zeros(GridSize+1,N+1);
        newZwakegridLEV = zeros(GridSize+1,N+1);
        if showVelocity
            figure;
            AxisLEV = axes;
            hold on
            ViewFlow(AxisLEV)
        end
        
        for jj = 1:(N+1)
            iPoints = nnz(XwakegridLEV(:,jj));
            for ii = 1:iPoints
                XYZl = [XwakegridLEV(ii,jj);YwakegridLEV(ii,jj);ZwakegridLEV(ii,jj)];
%                 XYZli = [1;-1;1].*XYZl;
                
                % Bound Rings Induced Velocity
                [UBInduced,VBInduced,WBInduced] = voringMatrix_FPbound(XYZl,M,N,circ,...
                    Xgrid,Ygrid,Zgrid);
%                 [UBInducedi,VBInducedi,WBInducedi] = voringMatrix_FPbound(XYZli,M,N,circ,...
%                     Xgrid,Ygrid,Zgrid);
                BoundUVWInduced = [sum(UBInduced(:));% + sum(UBInducedi(:));...
                    sum(VBInduced(:));% - sum(VBInducedi(:));...
                    sum(WBInduced(:))];% + sum(WBInducedi(:))];
                
                % Wake Rings Induced Velocity
                if NVRING == 0
                    WakeUVWInduced = [0;0;0];
                else
                    [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                        Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
%                     [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
%                         Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                    WakeUVWInduced = [sum(UWInduced(:));% + sum(UWInducedi(:));...
                        sum(VWInduced(:));% - sum(VWInducedi(:));...
                        sum(WWInduced(:))];% + sum(WWInducedi(:))];
                end
                
                % Wake Rings Induced Velocity
                % LEV Wake Rings Induced Velocity
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,GridSize,N,circwakeLEV(1:GridSize,:),...
                    XwakegridLEV(1:GridSize+1,:),YwakegridLEV(1:GridSize+1,:),ZwakegridLEV(1:GridSize+1,:));
%                 [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,GridSize,N,circwakeLEV(1:GridSize,:),...
%                     XwakegridLEV(1:GridSize+1,:),YwakegridLEV(1:GridSize+1,:),ZwakegridLEV(1:GridSize+1,:));
                WakeUVWInducedLEV = [sum(UWInduced(:));% + sum(UWInducedi(:));...
                    sum(VWInduced(:));% - sum(VWInducedi(:));...
                    sum(WWInduced(:))];% + sum(WWInducedi(:))];
                
                DistanceMatrix = sqrt((XYZl(1)-Xgrid(2:end-1,2:end-1)).^2+(XYZl(2)-Ygrid(2:end-1,2:end-1)).^2+(XYZl(3)-Zgrid(2:end-1,2:end-1)).^2);
                MinDistance = min(DistanceMatrix(:));
                if MinDistance < 5e-3
                    BoundUVWInduced(3) = 0;
                end
                CurrentVelocityTerm = [(BoundUVWInduced(1) + WakeUVWInduced(1) + WakeUVWInducedLEV(1));
                    (BoundUVWInduced(2) + WakeUVWInduced(2) + WakeUVWInducedLEV(2));
                    (BoundUVWInduced(3) + WakeUVWInduced(3) + WakeUVWInducedLEV(3))];
                if showVelocity
                    quiver3(XYZl(1),XYZl(2),XYZl(3),BoundUVWInduced(1),BoundUVWInduced(2),BoundUVWInduced(3),'m')
                    quiver3(XYZl(1),XYZl(2),XYZl(3),WakeUVWInduced(1),WakeUVWInduced(2),WakeUVWInduced(3),'c')
                    quiver3(XYZl(1),XYZl(2),XYZl(3),WakeUVWInducedLEV(1),WakeUVWInducedLEV(2),WakeUVWInducedLEV(3),'g')
                end
                newXwakegridLEV(ii,jj) = XYZl(1) + (CurrentVelocityTerm(1))*dt;
                newYwakegridLEV(ii,jj) = XYZl(2) + (CurrentVelocityTerm(2))*dt;
                newZwakegridLEV(ii,jj) = XYZl(3) + (CurrentVelocityTerm(3))*dt;
            end
        end
        if showVelocity
            view(4,10)
            hold off
        end
    end

    function LESPValue = CalcLESP(circulation)
        LESPValue = 1.13.*circulation./(-X0dot*c*(acos(1-2*cp/c)+sin(acos(1-2*cp/c))));
    end

    function LEVPosition
        NLEVRingST = NLEVRingST + ones(N,1);
        if NLEVRingST(1) == 1
            UVW_LEV = zeros(3,N+1);
            for ll = 1:(N+1)
                XYZ = tib*[x(1,ll); y(1,ll); z(1,ll)]+[X0;Y0;Z0];
                XwakegridLEV(1,ll) = XYZ(1);
                YwakegridLEV(1,ll) = XYZ(2);
                ZwakegridLEV(1,ll) = XYZ(3);
                UVW_LEV(:,ll) = [-X0dot;-Y0dot;-Z0dot]...
                    + tib*([-Q*z(1,ll) + R*y(1,ll);-R*x(1,ll) + P*z(1,ll);-P*y(1,ll) + Q*x(1,ll)]...
                    + [xdot(1,ll);0;zdot(1,ll)]);
            end
            XwakegridLEV(2,:) = XwakegridLEV(1,:) + 0.5.*UVW_LEV(1,:).*dt;
            YwakegridLEV(2,:) = YwakegridLEV(1,:) + 0.5.*UVW_LEV(2,:).*dt;
            ZwakegridLEV(2,:) = ZwakegridLEV(1,:) + 0.5.*UVW_LEV(3,:).*dt;
        else
            XwakegridLEV(3:end,:) = XwakegridLEV(2:end-1,:);
            YwakegridLEV(3:end,:) = YwakegridLEV(2:end-1,:);
            ZwakegridLEV(3:end,:) = ZwakegridLEV(2:end-1,:);
            for ll = 1:(N+1)
                XYZ = tib*[x(1,ll); y(1,ll); z(1,ll)]+[X0;Y0;Z0];
                XwakegridLEV(1,ll) = XYZ(1);
                YwakegridLEV(1,ll) = XYZ(2);
                ZwakegridLEV(1,ll) = XYZ(3);
            end
            XwakegridLEV(2,:) = 2/3.*XwakegridLEV(1,:) + 1/3.*XwakegridLEV(3,:);
            YwakegridLEV(2,:) = 2/3.*YwakegridLEV(1,:) + 1/3.*YwakegridLEV(3,:);
            ZwakegridLEV(2,:) = 2/3.*ZwakegridLEV(1,:) + 1/3.*ZwakegridLEV(3,:);
        end
        circwakeLEV(2:end,:) = circwakeLEV(1:end-1,:);
    end

    function ViewFlow(axisHandle)
        mesh(axisHandle,Xgrid,Ygrid,Zgrid,'EdgeColor',[0 0 0]);
        TEVR = NVRING+1;
        if TEVR > 1
            mesh(axisHandle,Xwakegrid(1:TEVR,:),Ywakegrid(1:TEVR,:),Zwakegrid(1:TEVR,:),'EdgeColor',[0 0 1]);ylim([0 0.15]);zlim([-0.1 0.1]);
        end
        LEVR = NLEVRingST+1;
        [XwakeLEV,YwakeLEV,ZwakeLEV] = FixLEVGrid(XwakegridLEV,YwakegridLEV,ZwakegridLEV,NLEVRingST,N);
        if any(LEVR>1)
            mesh(axisHandle,XwakeLEV(1:max(LEVR),:),YwakeLEV(1:max(LEVR),:),ZwakeLEV(1:max(LEVR),:),'EdgeColor',[1 0 0]);ylim([0 0.15]);zlim([-0.1 0.1]);
        end
        if any(pseudoVortexflag)
            for gg=1:numel(pseudoVortexflag)
                if pseudoVortexflag(gg)
                    Xpseudo = [XwakegridLEV(1,gg) XwakegridLEV(1,gg+1);Xwakegrid(1,gg) Xwakegrid(1,gg+1)];
                    Ypseudo = [YwakegridLEV(1,gg) YwakegridLEV(1,gg+1);Ywakegrid(1,gg) Ywakegrid(1,gg+1)];
                    Zpseudo = [ZwakegridLEV(1,gg) ZwakegridLEV(1,gg+1);Zwakegrid(1,gg) Zwakegrid(1,gg+1)];
                    mesh(axisHandle,Xpseudo,Ypseudo,Zpseudo,'EdgeColor',[0 1 0]);ylim([0 0.15]);zlim([-0.1 0.1]);
                end
            end
        end
    end

    function visualizeWakeDev
        subplot(2,1,1)
        TEVR = NVRING+1;
        if TEVR>1
            surf(Xwakegrid(1:TEVR,:),Ywakegrid(1:TEVR,:),Zwakegrid(1:TEVR,:));
            shading interp
            lighting flat
            ylim([0 0.15]);zlim([-0.1 0.1]);
        end
        hold on
        mesh(Xgrid,Ygrid,Zgrid,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
        material shiny
        camlight(0,90);
        LEVR = NLEVRingST+1;
        [XwakeLEV,YwakeLEV,ZwakeLEV] = FixLEVGrid(XwakegridLEV,YwakegridLEV,ZwakegridLEV,NLEVRingST,N);
        if any(LEVR>1)
            surf(XwakeLEV(1:max(LEVR),:),YwakeLEV(1:max(LEVR),:),ZwakeLEV(1:max(LEVR),:),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7,'FaceLighting','gouraud');
            ylim([0 0.15]);zlim([-0.1 0.1]);
        end
        if any(pseudoVortexflag)
            for nn=1:numel(pseudoVortexflag)
                if pseudoVortexflag(nn)
                    Xpseudo = [XwakegridLEV(1,nn) XwakegridLEV(1,nn+1);Xwakegrid(1,nn) Xwakegrid(1,nn+1)];
                    Ypseudo = [YwakegridLEV(1,nn) YwakegridLEV(1,nn+1);Ywakegrid(1,nn) Ywakegrid(1,nn+1)];
                    Zpseudo = [ZwakegridLEV(1,nn) ZwakegridLEV(1,nn+1);Zwakegrid(1,nn) Zwakegrid(1,nn+1)];
                    mesh(Xpseudo,Ypseudo,Zpseudo,'EdgeColor',[0 1 0]);ylim([0 0.15]);zlim([-0.1 0.1]);
                end
            end
        end
        hold off
        view(5,20)
        title('wake')
            
        subplot(2,1,2)
        [XLsp,YLsp] = meshgrid(linspace(0,c,M),linspace(0,b,N));
        surf(XLsp,YLsp,circ','EdgeColor','b','FaceColor','b');ylim([0 0.15]);zlim([-0.005 0.005])
        title('circulation')
        
        getframe;
    end
end
