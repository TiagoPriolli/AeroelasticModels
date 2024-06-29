function [coeff,localState,Force] = UVLM_FlatPlateLEV_Presc(x,y,z,xc,yc,zc,xcdot,zcdot,...
                                      xgrid,ygrid,zgrid,ni,taui,tauj,...
                                      localState,Parameters,Kinematics)

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

% Inertial to Body Transformation
Phi = Kinematics.Phi;
Theta = Kinematics.Theta;
Psi = Kinematics.Psi;

tbi = [1 0 0;0 cos(Phi) sin(Phi);0 -sin(Phi) cos(Phi)]*...
    ([cos(Theta) 0 -sin(Theta);0 1 0;sin(Theta) 0 cos(Theta)]*...
    [cos(Psi) sin(Psi) 0;-sin(Psi) cos(Psi) 0;0 0 1]);

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

% if any(NLEVRingST > 0)
%     ind = find(NLEVRingST > 0);
%     XwakegridLEV(1,[ind ind+1]) = Xgrid(1,[ind ind+1]);
%     YwakegridLEV(1,[ind ind+1]) = Ygrid(1,[ind ind+1]);
%     ZwakegridLEV(1,[ind ind+1]) = Zgrid(1,[ind ind+1]);
% end


%% Influence Coefficients
% Normal Velocity is a combination of self-induced velocity, kinematic
% velocity and wake induced velocity
A = InfluenceCoeff_FlatPlate;

%% Establish Right-Hand Side
pseudoVortexflag = zeros(1,N);
[RHS,uvw] = RightHandSide_FlatPlate;

%% Solve Linear Equations
circ = A\RHS;
circ = reshape(circ,M,N);

%% Verify LEV Generation

% Calculate LESP Distribution
A0 = CalcLESP(circ(1,:));
levflag = abs(A0) > LESPcrit;
pseudoVortexflag = levflag;

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
    XwakegridLEV(1,indG) = Xgrid(1,indG);
    YwakegridLEV(1,indG) = Ygrid(1,indG);
    ZwakegridLEV(1,indG) = Zgrid(1,indG);
    XwakegridLEV(2,indG) = 2/3.*Xgrid(1,indG) + 1/3.*XwakegridLEV(3,indG);
    YwakegridLEV(2,indG) = 2/3.*Ygrid(1,indG) + 1/3.*YwakegridLEV(3,indG);
    ZwakegridLEV(2,indG) = 2/3.*Zgrid(1,indG) + 1/3.*ZwakegridLEV(3,indG);
    if any(isFirst)
        XwakegridLEV(2,isFirst) = Xgrid(1,isFirst) - 0.2.*X0dot.*dt;
        YwakegridLEV(2,isFirst) = Ygrid(1,isFirst) - 0.2.*Y0dot.*dt;
        ZwakegridLEV(2,isFirst) = Zgrid(1,isFirst) - 0.2.*Z0dot.*dt;
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
    while any(levflag)
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
    
end
   
%% Secondary Computations
dP = zeros(M,N);
df = zeros(3,M,N);
dF = zeros(3,M,N);
ds = zeros(3,M,N);
dS = zeros(3,M,N);
%alpha = zeros(M,N);

circ(:,pseudoVortexflag) = circ(:,pseudoVortexflag) + ones(M,1).*circwakeLEV(1,pseudoVortexflag);

dPhidTaui = [circ(1,:); circ(2:M,:)-circ(1:M-1,:)]./cp;
dPhidTauj = [circ(:,1)-circ(:,1), circ(:,2:N)-circ(:,1:N-1)]./bp;

for j = 1:N
    for i = 1:M
        k = sub2ind([M N],i,j);
        
        dgamadt = (circ(i,j)-circd(i,j))./(dt);
        
        dP(i,j) = rho*(dot(uvw(:,k),taui(:,k).*dPhidTaui(i,j)) + ...
            dot(uvw(:,k),tauj(:,k).*dPhidTauj(i,j)) + ...
            dgamadt);
        
        df(:,i,j) = -(dP(i,j).*Sp).*ni(:,k); %Force Vector on Local Reference
        dF(:,i,j) = tib*df(:,i,j); %Force Vector on Global Reference
        ds(:,i,j) = -(pi*rho*c/M*(-X0dot)^2*A0(j)^2*bp).*[cos(localState.Lalpha(j));0;sin(localState.Lalpha(j))]; %Suction Force Vector on Local Reference
        dS(:,i,j) = tib*ds(:,i,j); %Suction Force Vector on Global Reference
        %alpha(i,j) = Theta;% + localState.Lalpha(j);
    end
end

dn = permute(df(3,:,:),[2,3,1]) + permute(ds(3,:,:),[2,3,1]); %Total normal force to panel
dM = reshape(dn(:).*(pvt*c - (xc(:)-cp/2)),M,N); %Panel's moment around strip's elastic axis
Mom = sum(dM,1)'; %Total moment around strip's elastic axis

dL = permute(dF(3,:,:),[2,3,1]) + permute(dS(3,:,:),[2,3,1]); %Total Lift force to panel
L = sum(dL,1)'; %Total Lift force to strip
% Verification Step
Normal = sum(permute(dF(3,:,:),[2,3,1]),1);
Suction = sum(permute(dS(3,:,:),[2,3,1]),1);
Force(1) = sum(Normal);
Force(2) = sum(Suction);
%

dD = permute(dF(1,:,:),[2,3,1]) + permute(dS(1,:,:),[2,3,1]); %Total F=Drag force to panel
D = sum(dD,1)'; %Total Drag force to strip

% FlightCoefficients
cL = L./((S/N)*(1/2*rho*Uref^2));
cD = D./((S/N)*(1/2*rho*Uref^2));
cM = Mom./((S/N)*c*(1/2*rho*Uref^2));

coeff = [cL,cD,cM,A0'];

circd = circ;

%% Vortex Wake Rollup

% Calculate Wake Rollup
[newXwakegrid,newYwakegrid,newZwakegrid] = UVLMWakeRollupNoSym;

if any(NLEVRingST > 0)
    
    % Calculate LEV Wake Rollup
    [newXwakegridLEV,newYwakegridLEV,newZwakegridLEV] = UVLMLEVWakeRollup;
    %
    %
%     subplot(2,1,1),mesh(Xgrid,Ygrid,Zgrid,'FaceColor','none','EdgeColor','b'),hold on
%                    mesh(Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:),'FaceColor','none','EdgeColor','r')
%                    scatter3(XwakegridLEV(XwakegridLEV~=0),YwakegridLEV(XwakegridLEV~=0),ZwakegridLEV(XwakegridLEV~=0),'g'),hold off
%     subplot(2,1,2),mesh(Xgrid,Ygrid,Zgrid,'FaceColor','none','EdgeColor','b'),hold on
%                    mesh(newXwakegrid,newYwakegrid,newZwakegrid,'FaceColor','none','EdgeColor','r')
%                    scatter3(newXwakegridLEV(newXwakegridLEV~=0),newYwakegridLEV(newXwakegridLEV~=0),newZwakegridLEV(newXwakegridLEV~=0),'g'),hold off
    %              
    %
    % Update LEV Position
    XwakegridLEV(1:max(NLEVRingST)+1,:) = newXwakegridLEV(1:max(NLEVRingST)+1,:);
    YwakegridLEV(1:max(NLEVRingST)+1,:) = newYwakegridLEV(1:max(NLEVRingST)+1,:);
    ZwakegridLEV(1:max(NLEVRingST)+1,:) = newZwakegridLEV(1:max(NLEVRingST)+1,:);
    
    % Delete LEV too far
    Lev2Far = abs(XwakegridLEV(max(NLEVRingST)+1,:) - Xgrid(M+1,:)) > c*10 ...
                    | XwakegridLEV(max(NLEVRingST)+1,:) == 0;
    if all(Lev2Far)
        %Lev2Delete = Lev2Far & XwakegridLEV(max(NLEVRingST)+1,:) ~= 0;
        XwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
        YwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
        ZwakegridLEV(max(NLEVRingST)+1,:) = zeros(1,N+1);
        circwakeLEV(max(NLEVRingST),:) = zeros(1,N);
        Ring2Delete = NLEVRingST == max(NLEVRingST);
        NLEVRingST(Ring2Delete) = NLEVRingST(Ring2Delete) - 1;
    end
end

% Update Global Wake Position and Circulation in the Matrices
Xwakegrid(2:NVRING+2,:) = newXwakegrid(1:NVRING+1,:);
Ywakegrid(2:NVRING+2,:) = newYwakegrid(1:NVRING+1,:);
Zwakegrid(2:NVRING+2,:) = newZwakegrid(1:NVRING+1,:);
circwake(2:NVRING+1,:) = circwake(1:NVRING,:);
circwake(1,:) = circ(M,:);
if any(pseudoVortexflag)
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
            xyzi = [xc(kk);-yc(kk);zc(kk)];
            [uRing,vRing,wRing] = voringMatrix_FPbound...
                 (xyz,M,N,circulation,xgrid,ygrid,zgrid);
            [uRingi,vRingi,wRingi] = voringMatrix_FPbound...
                   (xyzi,M,N,circulation,xgrid,ygrid,zgrid);
            
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
            
            % Point of Interest
            XYZl = [XYZc(kk,1);XYZc(kk,2);XYZc(kk,3)];
            XYZli = [XYZc(kk,1);-XYZc(kk,2);XYZc(kk,3)];
            
            % Wake induced Velocity
            if NVRING == 0
                uvwTEV = [0;0;0];
            else
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,NVRING,N,circwake(1:NVRING,:),...
                    Xwakegrid(1:NVRING+1,:),Ywakegrid(1:NVRING+1,:),Zwakegrid(1:NVRING+1,:));
                uvwTEV = tbi*[sum(URing(:)) + sum(URingi(:));...
                              sum(VRing(:)) - sum(VRingi(:));...
                              sum(WRing(:)) + sum(WRingi(:))];
            end
            
            % LEV induced Velocity
            if all(NLEVRingST == 0)
                uvwLEV = [0;0;0];
            else
                Size = max(NLEVRingST);
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,Size,N,circwakeLEV(1:Size,:),...
                    XwakegridLEV(1:Size+1,:),Ywakegrid(1:Size+1,:),Zwakegrid(1:Size+1,:));
                [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,Size,N,circwakeLEV(1:Size,:),...
                    XwakegridLEV(1:Size+1,:),Ywakegrid(1:Size+1,:),Zwakegrid(1:Size+1,:));
                uvwLEV = tbi*[sum(URing(:)) + sum(URingi(:));...
                              sum(VRing(:)) - sum(VRingi(:));...
                              sum(WRing(:)) + sum(WRingi(:))];
            end
            
            % Kinematics Induced Velocity
            uvwt = tbi*[-X0dot;-Y0dot;-Z0dot]...
                + [-Q*zc(kk) + R*yc(kk);-R*xc(kk) + P*zc(kk);-P*yc(kk) + Q*xc(kk)]...
                + [xcdot(kk);0;zcdot(kk)];
            
            % Unbounded Induced Velocity
            uvw(:,kk) =  uvwt + uvwTEV + uvwLEV;
            
            % Pseudo Vortex Induced Velocity
            if any(pseudoVortexflag)
                XpseudoGrid = [Xgrid(1,:);Xgrid(end,:)];
                YpseudoGrid = [Ygrid(1,:);Ygrid(end,:)];
                ZpseudoGrid = [Zgrid(1,:);Zgrid(end,:)];
                pseudocirc = zeros(1,N);
                pseudocirc(1,pseudoVortexflag) = circwakeLEV(1,pseudoVortexflag);
                [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,1,N,pseudocirc,...
                    XpseudoGrid,YpseudoGrid,ZpseudoGrid);
                [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,1,N,pseudocirc,...
                    XpseudoGrid,YpseudoGrid,ZpseudoGrid);
                uvwPSEUDO = tbi*[sum(URing(:)) + sum(URingi(:));...
                                    sum(VRing(:)) - sum(VRingi(:));...
                                    sum(WRing(:)) + sum(WRingi(:))];
            else
                uvwPSEUDO = 0;
            end
            
            % Right-Hand Side
            uvwrhs = uvw(:,kk) + uvwPSEUDO;
            RHS(kk) = dot(-uvwrhs,ni(:,kk));
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
            
            % LEV Wake Rings Induced Velocity
            if all(NLEVRingST == 0)
                WakeUVWInducedLEV = [0;0;0];
            else
                size = max(NLEVRingST);
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,size,N,circwakeLEV(1:size,:),...
                    XwakegridLEV(1:size+1,:),Ywakegrid(1:size+1,:),Zwakegrid(1:size+1,:));
                [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,size,N,circwakeLEV(1:size,:),...
                    XwakegridLEV(1:size+1,:),Ywakegrid(1:size+1,:),Zwakegrid(1:size+1,:));
                WakeUVWInducedLEV = [sum(UWInduced(:)) + sum(UWInducedi(:));...
                                     sum(VWInduced(:)) - sum(VWInducedi(:));...
                                     sum(WWInduced(:)) + sum(WWInducedi(:))];
            end
            
            %Pseudo-Vortex
%             if any(pseudoVortexflag)
%                 XpseudoGrid = [Xgrid(1,:);Xgrid(end,:)];
%                 YpseudoGrid = [Ygrid(1,:);Ygrid(end,:)];
%                 ZpseudoGrid = [Zgrid(1,:);Zgrid(end,:)];
%                 pseudocirc = circwakeLEV(1,:);
%                 [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,1,N,pseudocirc,...
%                     XpseudoGrid,YpseudoGrid,ZpseudoGrid);
%                 [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,1,N,pseudocirc,...
%                     XpseudoGrid,YpseudoGrid,ZpseudoGrid);
%                 uvwrhsPSEUDO = [sum(URing(:)) + sum(URingi(:));...
%                                 sum(VRing(:)) - sum(VRingi(:));...
%                                 sum(WRing(:)) + sum(WRingi(:))];
%             else
%                  uvwrhsPSEUDO = [0;0;0];
%             end
            uvwrhsPSEUDO = [0;0;0];
            
            % Wake Roolup
            newXwakegrid(ii,jj) = XYZl(1) + (BoundUVWInduced(1) + WakeUVWInduced(1) + WakeUVWInducedLEV(1) + uvwrhsPSEUDO(1))*dt;
            newYwakegrid(ii,jj) = XYZl(2) + (BoundUVWInduced(2) + WakeUVWInduced(2) + WakeUVWInducedLEV(2) + uvwrhsPSEUDO(2))*dt;
            newZwakegrid(ii,jj) = XYZl(3) + (BoundUVWInduced(3) + WakeUVWInduced(3) + WakeUVWInducedLEV(3) + uvwrhsPSEUDO(3))*dt;
        end
    end

    function [newXwakegridLEV,newYwakegridLEV,newZwakegridLEV] = UVLMLEVWakeRollup
        
        size = max(NLEVRingST);
        newXwakegridLEV = zeros(size+1,N+1);
        newYwakegridLEV = zeros(size+1,N+1);
        newZwakegridLEV = zeros(size+1,N+1);
        
        for jj = 1:(N+1)
            iPoints = nnz(XwakegridLEV(:,jj));
            for ii = 1:iPoints
                XYZl = [XwakegridLEV(ii,jj);YwakegridLEV(ii,jj);ZwakegridLEV(ii,jj)];
                XYZli = [1;-1;1].*XYZl;

                % Bound Rings Induced Velocity
                [UBInduced,VBInduced,WBInduced] = voringMatrix_FPbound(XYZl,M,N,circ,...
                    Xgrid,Ygrid,Zgrid);
                [UBInducedi,VBInducedi,WBInducedi] = voringMatrix_FPbound(XYZli,M,N,circ,...
                    Xgrid,Ygrid,Zgrid);
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

                % Wake Rings Induced Velocity
                % LEV Wake Rings Induced Velocity
                [UWInduced,VWInduced,WWInduced] = voringMatrix_FlatPlate(XYZl,size,N,circwakeLEV(1:size,:),...
                            XwakegridLEV(1:size+1,:),Ywakegrid(1:size+1,:),Zwakegrid(1:size+1,:));
                [UWInducedi,VWInducedi,WWInducedi] = voringMatrix_FlatPlate(XYZli,size,N,circwakeLEV(1:size,:),...
                            XwakegridLEV(1:size+1,:),Ywakegrid(1:size+1,:),Zwakegrid(1:size+1,:));
                WakeUVWInducedLEV = [sum(UWInduced(:)) + sum(UWInducedi(:));...
                                     sum(VWInduced(:)) - sum(VWInducedi(:));...
                                     sum(WWInduced(:)) + sum(WWInducedi(:))];
                                 
                %Pseudo-Vortex
%                 if any(pseudoVortexflag)
%                     XpseudoGrid = [Xgrid(1,:);Xgrid(end,:)];
%                     YpseudoGrid = [Ygrid(1,:);Ygrid(end,:)];
%                     ZpseudoGrid = [Zgrid(1,:);Zgrid(end,:)];
%                     pseudocirc = circwakeLEV(1,:);
%                     [URing,VRing,WRing] = voringMatrix_FlatPlate(XYZl,1,N,pseudocirc,...
%                         XpseudoGrid,YpseudoGrid,ZpseudoGrid);
%                     [URingi,VRingi,WRingi] = voringMatrix_FlatPlate(XYZli,1,N,pseudocirc,...
%                         XpseudoGrid,YpseudoGrid,ZpseudoGrid);
%                     uvwrhsPSEUDO = [sum(URing(:)) + sum(URingi(:));...
%                                     sum(VRing(:)) - sum(VRingi(:));...
%                                     sum(WRing(:)) + sum(WRingi(:))];
%                 else
%                      uvwrhsPSEUDO = [0;0;0];
%                 end
                uvwrhsPSEUDO = [0;0;0];

                % Wake Roolup
                newXwakegridLEV(ii,jj) = XYZl(1) + (BoundUVWInduced(1) + WakeUVWInduced(1) + WakeUVWInducedLEV(1) + uvwrhsPSEUDO(1))*dt;
                newYwakegridLEV(ii,jj) = XYZl(2) + (BoundUVWInduced(2) + WakeUVWInduced(2) + WakeUVWInducedLEV(2) + uvwrhsPSEUDO(2))*dt;
                newZwakegridLEV(ii,jj) = XYZl(3) + (BoundUVWInduced(3) + WakeUVWInduced(3) + WakeUVWInducedLEV(3) + uvwrhsPSEUDO(3))*dt;
            end
        end
    end

    function LESPValue = CalcLESP(circulation)
        LESPValue = 1.13.*circulation./(-X0dot*c*(acos(1-2*cp/c)+sin(acos(1-2*cp/c))));
        %LESPValue = 1.4.*circulation./(-X0dot*c*(acos(1-2*cp/c)+sin(acos(1-2*cp/c))));
    end
end
