%% Test Setup
Uref = 10; % choose flow velocity (must be one of the speeds for which matrices have been determined)

%% Dados Estruturais
ballast = 5;
load flatplate_400_27_8_ballast_1_offset_5_nelem_70_aeroelast_ns_20_nU_40


%%
c = 2*vec_b(1);

% bulding structure of step changing parameters
global fix_p sch_p stateout

% -> fixed parameters
fix_p.c = c;
fix_p.Uref = Uref;
fix_p.rho = rho;
fix_p.ns = ns;
fix_p.iter_max=100;  %Max. iterations
fix_p.v_core=0.02;   %Non-dimensional core radius of point vortices
fix_p.n_div=70;      %No. of divisions along chord on airfoil
fix_p.n_aterm=45;  %Number of fourier terms used to compute vorticity at a location on chord
fix_p.del_dist=10;
fix_p.expect_vort = 1000;
fix_p.pvt = 0.5;      %Pivot(0-1)
fix_p.re_ref = 10000;   %Reference Reynolds number
fix_p.lesp_crit = 0.02;  %Critical LESP
fix_p.dt = 0.015; % dt_star = dt*U/c
% Defining Divisions
fix_p.dtheta = pi/(fix_p.n_div-1);
fix_p.theta = ((1:fix_p.n_div)'-1)*fix_p.dtheta;
fix_p.x = (c/2)*(1-cos(fix_p.theta));
% -> variable parameters with strip and time step
for jj = 1:ns
    sch_p(jj).ntev = 0;
    sch_p(jj).nlev = 0;
    sch_p(jj).tev = zeros(fix_p.expect_vort,3);
    sch_p(jj).lev = zeros(fix_p.expect_vort,3);
    sch_p(jj).dist_wind = 0;
    sch_p(jj).kelv_enf = 0;
    sch_p(jj).aterm_prev = zeros(4,1);
    sch_p(jj).levflag = 0;
    sch_p(jj).a0 = 0;
end

dt = fix_p.dt*2*vec_b(1)/Uref;
tfinal = 4;
niter = round(tfinal/dt);

q = 1*pi/180; % step in the aoa
t0 = 0;

%% Run LDVM Testbench Model
X0 = zeros(2*nm,1);
parpool(3);
Xp = ode4state(@(t,X,tag) Dynamics_LDVM_flatplate(t,X,tag,q,nm,ns,mmod(find(U==Uref)),...
                                Uref,rho,Phi,vec_b,vec_dy),0,dt,tfinal,X0);
tp = (t0:dt:tfinal)';
Xp = (reshape(Xp,2*nm,size(tp,1)))';
etap = Xp(:,1:nm);
etadp = Xp(:,nm+1:2*nm);

delete(gcp('nocreate'));

%% Run linear Wagner Model
X0 = zeros(2*nm+2*ns,1);
Xlinode = ode4(@(t,X)dynamics_lin(t,X,q,nm,ns,mmod(find(U==Uref))),...
                                  0,dt,tfinal,X0);
tode = (t0:dt:tfinal)';
Xlinode = (reshape(Xlinode,2*nm+2*ns,size(tode,1)))';
eta_linode = Xlinode(:,1:nm);
etad_linode = Xlinode(:,nm+1:2*nm);

%%
timestamp = datestr(now(),'_ddmm_HHMMSS');
file = sprintf(['LDVMtestbench' timestamp]);
save(file)

