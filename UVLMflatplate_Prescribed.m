global fix_p_uvlm sch_p_uvlm

%% Model Parameters
c = 0.10;                 % Chord
b = 0.15*10;                 % Half-Span: Model uses symmetry to calculate only half the forces.
rho = 1.225;             % Air density
Uref = 0.1;              % Freestream Velocity
lespcrit = 100;        % LESP Critical (use anything over 100 for "pure" UVLM)

Re = Uref*c/18.5e-6;     % Reynolds

M = 20;                  % Number of Chordwise Panels
N = 30*10;                  % Number of Spanwise Panels

% Generate code for faster performance
CodegenVRing(M,N)

% Turn on vortex visualization
viz = false;
showVelocity = false;

%% Simulation Parameters

ExpectVort = 1000;      % Prealocation for number of chorwise vortices on wake
dt_star = 1/M;          % Adimensional time step
dt = dt_star*c/Uref;    % Dimensional time step (s)
t0_star = 0;            % Adimensional start time
t0 = 0;                 % Dimensional start time (s)
tfinal = 7;             % Dimensional final time (s)
tfinal_star = tfinal*Uref/c; % Adimensional final time

alphaTest = [3 25 45];     % Maximum angle of attack (deg)
pvtTest = [0 0.25 0.75]; % Pivot position
alphaMAX = alphaTest(3);
pvt = pvtTest(1);
%% UVLM
% Run selected Test

disp('Running UVLM...')

% bulding structure of step changing parameters
UVLMFlatPlateSetup(c,b,M,N,Uref,ExpectVort,pvt);

fix_p_uvlm.rho = rho;
fix_p_uvlm.dt = dt_star;
fix_p_uvlm.expect_vort = ExpectVort;
fix_p_uvlm.LESPcrit = lespcrit;

niter = round(tfinal_star/dt_star);

%% Grid Data Extraction
xc = fix_p_uvlm.xc;
yc = fix_p_uvlm.yc;
zc = fix_p_uvlm.zc;
xcdot = zeros(size(zc));
zcdot = zeros(size(zc));
x = fix_p_uvlm.x;
y = fix_p_uvlm.y;
z = fix_p_uvlm.z;
xdot = zeros(M+1,N+1);
zdot = zeros(M+1,N+1);
xgrid = fix_p_uvlm.xgrid;
ygrid = fix_p_uvlm.ygrid;
zgrid = fix_p_uvlm.zgrid;

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
Lalpha = zeros(N,1);
OriPos.xc = xc;
OriPos.x = x;
OriPos.xgrid = xgrid;


%% Run UVLM Prescribed Simulation
i_step = 1;
limit = 1;
cl = zeros(N,niter);
cm = zeros(N,niter);
a0 = zeros(N,niter);
Theta = zeros(niter,1);
Normal = zeros(niter,1);
Suction = zeros(niter,1);
Q = zeros(niter,1);

XWakeVideo = zeros(200,N+1,niter);
XLEVWakeVideo = zeros(200,N+1,niter);
XgridVideo = zeros(M+1,N+1,niter);
YWakeVideo = zeros(200,N+1,niter);
YLEVWakeVideo = zeros(200,N+1,niter);
YgridVideo = zeros(M+1,N+1,niter);
ZWakeVideo = zeros(200,N+1,niter);
ZLEVWakeVideo = zeros(200,N+1,niter);
ZgridVideo = zeros(M+1,N+1,niter);

tvec_star = (t0_star:dt_star:tfinal_star)';
tvec = (t0:dt:tfinal)';
G = Eldredge(c,Uref,tvec,1,3,4,6);
maxG = max(G);

fprintf('Sim start: End time = %f, Num of steps = %f\n',tfinal,numel(tvec_star));
for t_star=t0_star:dt_star:tfinal_star
    
    t = t_star*c/Uref;
    %% Kinematics
    Kinematics.X0 = -Uref*t;
    Kinematics.Y0 = 0;
    Kinematics.Z0 = 0;
    
    Kinematics.X0dot = -Uref;
    Kinematics.Y0dot = 0;
    Kinematics.Z0dot = 0;
    
    Kinematics.Phi = 0;
    Kinematics.Psi = 0;
    
    Kinematics.P = 0;
    Kinematics.R = 0;
    
    Kinematics.Theta = alphaMAX*pi/180*G(i_step)/maxG;
    if i_step == numel(G)
        Kinematics.Q = 0;
    else
        Kinematics.Q = (alphaMAX*pi/180*G(i_step+1)/maxG-alphaMAX*pi/180*G(i_step)/maxG)/dt;
    end
%     Kinematics.Theta = 0;
%     Kinematics.Q = 0;
%     [x,xdot,z,zdot,xgrid,zgrid,xc,zc,xcdot,zcdot,ni,taui,tauj,Lalpha] = UpdateUVLM_FlatPlate(M,N,alphaMAX,G,maxG,...
%                      pvt,c,ygrid,i_step,dt,OriPos);
    
    sch_p_uvlm.Lalpha = Lalpha + Kinematics.Theta;
    sch_p_uvlm.PreviousVelocityTerm = [];
    sch_p_uvlm.PreviousVelocityTermLEV = [];
    sch_p_uvlm.StartLevShed = false;
    
    %% Run UVLM
    [coeff,sch_p_uvlm,forces] = UVLM_FlatPlateLEV(x,y,z,xdot,zdot,xc,yc,zc,xcdot,zcdot,...
        xgrid,ygrid,zgrid,ni,taui,tauj,...
        sch_p_uvlm,fix_p_uvlm,Kinematics,viz,showVelocity);
    
    % Calculations in the model use symmetry.
    cl(:,i_step) = coeff(:,1);
    cm(:,i_step) = coeff(:,3);
    a0(:,i_step) = coeff(:,4);
    Normal(i_step) = forces(1);
    Suction(i_step) = forces(2);
    Theta(i_step,1) = sch_p_uvlm.Lalpha(1);
    Q(i_step,1) = Kinematics.Q;
    i_step = i_step + 1;
    
    % Saving grid information for animation
    XWakeVideo(:,:,i_step) = sch_p_uvlm.Xwakegrid(1:200,:);
    XLEVWakeVideo(:,:,i_step) = sch_p_uvlm.XwakegridLEV(1:200,:);
    XgridVideo(:,:,i_step) = sch_p_uvlm.Xgrid;
    YWakeVideo(:,:,i_step) = sch_p_uvlm.Ywakegrid(1:200,:);
    YLEVWakeVideo(:,:,i_step) = sch_p_uvlm.YwakegridLEV(1:200,:);
    YgridVideo(:,:,i_step) = sch_p_uvlm.Ygrid;
    ZWakeVideo(:,:,i_step) = sch_p_uvlm.Zwakegrid(1:200,:);
    ZLEVWakeVideo(:,:,i_step) = sch_p_uvlm.ZwakegridLEV(1:200,:);
    ZgridVideo(:,:,i_step) = sch_p_uvlm.Zgrid;
    
    if i_step > limit
        text = sprintf('Time = %f, Step = %f',t,i_step);
        disp(text)
        limit = limit + 50;
    end
    
end

%%
qinf = rho*Uref^2/2;
l_st = cl.*qinf.*(c*b(1)/N);
L = sum(l_st,1);
CL = L./(qinf*(c*b(1)));
m_st = cm.*qinf.*(c^2*b(1)/N);
M = sum(m_st,1);
CM = M./(qinf*(c^2*b(1)));

timestamp = datestr(now(),'dd_mm_HHMMSS');
file = sprintf(['UVLMPresc_' timestamp]);
save(file)

figure, 
subplot(1,2,1), 
plot(tvec,CL,'g','LineWidth',2),hold on
xlim([tvec(1),tvec(end)])
ylim([-1,5]),yticks(-1:1:5)
ylabel('C_L'),xlabel('time(s)'),title('Leading-edge 45 degree pitch-hold-return')
subplot(1,2,2),
plot(tvec,CM,'g','LineWidth',2),hold on
xlim([tvec(1),tvec(end)])
ylim([-2.5,1]),yticks(-2.5:0.5:1)
ylabel('C_M'),xlabel('time(s)')

%%
function UVLMFlatPlateSetup(c,b,M,N,Uref,ExpectVort,pvt)
global fix_p_uvlm sch_p_uvlm

fix_p_uvlm.Uref = Uref;
fix_p_uvlm.M = M;
fix_p_uvlm.N = N;
fix_p_uvlm.c = c;
fix_p_uvlm.b = b;
fix_p_uvlm.S = b*c;
fix_p_uvlm.pvt = pvt;

[x0,y0,z0,xc,yc,zc,Sp,cp,bp,~,x,y,z] = flatgridpanels(c,b,M,N);

xgrid = zeros(M+1,N+1);
ygrid = zeros(M+1,N+1);
zgrid = zeros(M+1,N+1);
xgrid(1:M,:) = x0;
ygrid(1:M,:) = y0;
zgrid(1:M,:) = z0;

fix_p_uvlm.xc = xc(:);
fix_p_uvlm.yc = yc(:);
fix_p_uvlm.zc = zc(:);
fix_p_uvlm.x = x;
fix_p_uvlm.y = y;
fix_p_uvlm.z = z;
fix_p_uvlm.Sp = Sp;
fix_p_uvlm.cp = cp;
fix_p_uvlm.bp = bp;
fix_p_uvlm.xgrid = xgrid;
fix_p_uvlm.ygrid = ygrid;
fix_p_uvlm.zgrid = zgrid;

sch_p_uvlm.NVRING = 0;
sch_p_uvlm.circwake = zeros(ExpectVort,N);
sch_p_uvlm.Xwakegrid = zeros(ExpectVort+1,N+1);
sch_p_uvlm.Ywakegrid = zeros(ExpectVort+1,N+1);
sch_p_uvlm.Zwakegrid = zeros(ExpectVort+1,N+1);
sch_p_uvlm.NLEVRing = zeros(N,1);
sch_p_uvlm.circwakeLEV = zeros(ExpectVort,N);
sch_p_uvlm.XwakegridLEV = zeros(ExpectVort+1,N+1);
sch_p_uvlm.YwakegridLEV = zeros(ExpectVort+1,N+1);
sch_p_uvlm.ZwakegridLEV = zeros(ExpectVort+1,N+1);
sch_p_uvlm.circd = zeros(M,N);

%CodeGenerationNoSym(M,N,ExpectVort);
end

%%

function [x,xdot,z,zdot,xgrid,zgrid,xc,zc,xcdot,zcdot,ni,taui,tauj,Lalpha] = UpdateUVLM_FlatPlate(M,N,alphaMAX,G,maxG,...
    pvt,c,ygrid,i_step,dt,OriPos)
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
ni = zeros(3,M*N);
taui = zeros(3,M*N);
tauj = zeros(3,M*N);
xOri = OriPos.x;
xgridOri = OriPos.xgrid;
xcOri = OriPos.xc;

alphastrip = -alphaMAX*pi/180*G(i_step)/maxG;
if i_step == numel(G)
    alphadotstrip = (alphaMAX*pi/180*G(i_step)/maxG-alphaMAX*pi/180*G(i_step)/maxG)/dt;
else
    alphadotstrip = (alphaMAX*pi/180*G(i_step+1)/maxG-alphaMAX*pi/180*G(i_step)/maxG)/dt;
end

% Update each strip
% Update grid
for kk = 1:N+1
    
    z(:,kk) = (xOri(:,kk) - pvt*c)*sin(alphastrip);
    zdot(:,kk) = (xOri(:,kk) - pvt*c)*cos(alphastrip)*alphadotstrip;
    zgrid(:,kk) = (xgridOri(:,kk) - pvt*c)*sin(alphastrip);
    
    x(:,kk) = xOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    xdot(:,kk) = -xOri(:,kk).*sin(alphastrip).*alphadotstrip + pvt*c*sin(alphastrip)*alphadotstrip;
    xgrid(:,kk) = xgridOri(:,kk)*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    
end

% Update Control Points
for kk = 1:N
    
    cini = M*(kk-1)+1;
    cend = M*kk;
    
    zc(cini:cend,1) = (xcOri(cini:cend,1) - pvt*c).*sin(alphastrip);
    xc(cini:cend,1) = xcOri(cini:cend,1).*cos(alphastrip) + (1-cos(alphastrip))*pvt*c;
    
    zcdot(cini:cend,1) = (xcOri(cini:cend,1) - pvt*c).*...
        cos(alphastrip).*alphadotstrip;
    xcdot(cini:cend,1) = -xcOri(cini:cend,1).*sin(alphastrip).*alphadotstrip + pvt*c*sin(alphastrip)*alphadotstrip;
    
    Lalpha(kk) = alphaMAX*pi/180*G(i_step)/maxG;
end

% Update Normal Vectors
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