%% Test Setup
TestType = 5;
Uref = 10; % choose flow velocity (must be one of the speeds for which matrices have been determined)

%% Dados Estruturais
switch TestType
    case 1
        % Ballast Centralizado
        ballast = 0;
        load flatplate_350_40_8_ballast_1_offset_0_nelem_70_aeroelast_ns_10_nU_40
        
    case 2
        % Ballast 5mm
        ballast = 5;
        load flatplate_350_40_8_ballast_1_offset_5_nelem_70_aeroelast_ns_10_nU_40
        
    case 3
        % Ballast 10mm
        ballast = 10;
        load flatplate_350_40_8_ballast_1_offset_10_nelem_70_aeroelast_ns_10_nU_40
        
    case 4
        % Ballast 15mm
        ballast = 15;
        load flatplate_350_40_8_ballast_1_offset_15_nelem_70_aeroelast_ns_10_nU_40
        
    case 5
    % Flutter Thin Wing
    % Ballast 5mm
    ballast = 5;
    load flatplate_400_27_8_ballast_1_offset_5_nelem_70_aeroelast_ns_20_nU_40
end

% load partial_third_L002 sch_p fix_p tp etap etadp
% etap_old = [etap etadp];
% tp_old = tp;
% X0 = etap_old(end,:)';
load partial_FourFive sch_p fix_p t yout
X0vec = (reshape(yout,8,29632));
X0 = X0vec(:,end);

for i = 1:ns
    sch_p(i).tev = [sch_p(i).tev;zeros(500,3)];
    sch_p(i).lev = [sch_p(i).lev;zeros(500,3)];
end
    
%%
dt = fix_p.dt*fix_p.c/fix_p.Uref;
%0-4: First
%4-10: Second
%10-15: Fourth
%15-16: Fifth
%16-20: Sixth
t0 = t;
tfinal = 9.2;

q = 1*pi/180; % step in the aoa

clear t X0vec yout

%%
parpool(3);
Xp = ode4state(@(t,X,tag) Dynamics_LDVM_flatplate(t,X,tag,q,nm,ns,mmod(find(U==Uref)),...
                                Uref,rho,Phi,vec_b,vec_dy),t0,dt,tfinal,X0);
tp = (t0:dt:tfinal)';
Xp = (reshape(Xp,2*nm,size(tp,1)))';
etap = Xp(:,1:nm);
etadp = Xp(:,nm+1:2*nm);

delete(gcp('nocreate'));

%%
global stateout
timestamp = datestr(now(),'_ddmm_HHMMSS');
file = sprintf(['LDVMtestbench' timestamp]);
save(file)

