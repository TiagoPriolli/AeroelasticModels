%% Test Setup
TestType = 2;
Uref = 10; % choose flow velocity (must be one of the speeds for which matrices have been determined)

% bulding structure of step changing parameters
global fix_p_uvlm sch_p_uvlm stateout qout

%% Dados Estruturais
switch TestType
    case 1
        % More Flexible Wing - Flutter
        ballast = 5;
        load flatplate_400_27_8_ballast_1_offset_5_nelem_70_aeroelast_ns_10_nU_40
        
    case 2
        % More Flexible Wing - Flutter
        ballast = 5;
        load flatplate_400_27_8_ballast_1_offset_5_nelem_70_aeroelast_ns_20_nU_40
        
    case 3
        % More Flexible Wing - Flutter
        ballast = 5;
        load flatplate_400_27_8_ballast_1_offset_5_nelem_70_aeroelast_ns_30_nU_40
end

%% Model Parameters
c = 2*vec_b(1);
b = vec_y(end) + vec_dy(end)/2;
M = 10;
N = ns;
ExpectVort = 1000;
dt_star = 1/M;%0.015;
pvt = 0.5;

% -> UVLM SETUP
UVLMFlatPlateSetup(c,b,M,N,Uref,ExpectVort,pvt);
CodegenVRing(M,N)

fix_p_uvlm.rho = rho;
fix_p_uvlm.dt = dt_star;
fix_p_uvlm.expect_vort = ExpectVort;
fix_p_uvlm.LESPcrit = 0.03;

t0 = 0;
dt = dt_star*2*vec_b(1)/Uref;
tfinal = 10;
niter = round(tfinal/dt);

%% Kinematics
Kinematics.X0 = -Uref*t0;
Kinematics.Y0 = 0;
Kinematics.Z0 = 0;

Kinematics.X0dot = -Uref;
Kinematics.Y0dot = 0;
Kinematics.Z0dot = 0;

Kinematics.Phi = 0;
Kinematics.Theta = 1*pi/180; % step in the aoa
Kinematics.Psi = 0;

Kinematics.P = 0;
Kinematics.Q = 0;
Kinematics.R = 0;

%%
disp('Running UVLM...')
X0 = zeros(2*nm,1);
X = ode4Full(@(t,X,tag) Dynamics_UVLM_flatplate_LEV(t,X,tag,nm,ns,mmod(find(U==Uref)),...
    Uref,Phi,vec_b,vec_dy,Kinematics),0,dt,tfinal,X0);
t = (t0:dt:tfinal)';
eta = X(:,1:nm);
etad = X(:,nm+1:2*nm);
StatesUVLM = stateout;
disp('Complete!')

%%
% disp('Running Wagner...')
% X0 = zeros(2*nm+2*ns,1);
% Xlinode = ode4Wgner(@(t,X)dynamics_lin(t,X,q,nm,ns,mmod(find(U==Uref))),...
%     0,dt,tfinal,X0);
% tode = (t0:dt:tfinal)';
% Xlinode = (reshape(Xlinode,2*nm+2*ns,size(tode,1)))';
% eta_linode = Xlinode(:,1:nm);
% etad_linode = Xlinode(:,nm+1:2*nm);
% StateWagner = (reshape(qout,nm,size(tode,1)))';
% disp('Complete!')

%%
timestamp = datestr(now(),'dd_mm_HHMMSS');
q = Kinematics.Theta;
file = sprintf(['UVLM_U%d_Q%d_' timestamp],Uref,q*180/pi);
save(file)

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