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

M = 10;
N = ns;
CodegenVRing(M,N)

%% Model Parameters States

% Results Summary
% PartialStopLEV00 -> partial results -> t=0 to t=5.4
% FinalStopLEV00 -> final results -> t=5.4 to t=10
% FinalStopLEV01 -> final results -> t=10 to t=20
% FinalStopLEV02 -> final results -> t=20 to t=40
% PartialStopLEV01 -> partial results -> t=40 to t=42.8 (expanded size of
% LEV matrix)
% PartialStopLEV02 -> partial results -> t=42.8 to t=43

% Start from Partial Results
% load PartialStopLEV02 t fix_p_uvlm sch_p_uvlm yout
% dt = fix_p_uvlm.dt*fix_p_uvlm.c/fix_p_uvlm.Uref;
% t0 = t;
% lastpoint = find(yout(:,2)~=0,1,'last');
% X0 = yout(lastpoint,:)';
% sch_p_uvlm.StartLevShed = false;

% Start from final results
load FinalStopLEV04 fix_p_uvlm sch_p_uvlm eta etad t
sch_p_uvlm.XwakegridLEV = sch_p_uvlm.XwakegridLEV;
sch_p_uvlm.YwakegridLEV = sch_p_uvlm.YwakegridLEV;
sch_p_uvlm.ZwakegridLEV = sch_p_uvlm.ZwakegridLEV;
sch_p_uvlm.circwakeLEV = sch_p_uvlm.circwakeLEV;
dt = fix_p_uvlm.dt*fix_p_uvlm.c/fix_p_uvlm.Uref;
X0 = [eta(end,:) etad(end,:)]';
t0 = t(end,1);
%sch_p_uvlm.StartLevShed = false;
clearvars t eta etad

%% New Simulation Scenario
tfinal = 55;

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
X = ode4Full(@(t,X,tag) Dynamics_UVLM_flatplate_LEV(t,X,tag,nm,ns,mmod(find(U==Uref)),...
    Uref,Phi,vec_b,vec_dy,Kinematics),t0,dt,tfinal,X0);
t = (t0:dt:tfinal)';
eta = X(:,1:nm);
etad = X(:,nm+1:2*nm);
StatesUVLM = stateout;
disp('Complete!')

%%
% disp('Running Wagner...')
% X0 = Xlinode(end,:);
% Xlinode = ode4Wgner(@(t,X)dynamics_lin(t,X,q,nm,ns,mmod(find(U==Uref))),...
%     t0,dt,tfinal,X0);
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


