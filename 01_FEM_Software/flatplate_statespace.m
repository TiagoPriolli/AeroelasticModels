%% Routine to calculate all modal matrices, state-space matrices
clear all
% loading modal shapes
%file_modalshapes = 'flatplate_350_40_8_ballast_1_offset_15_nelem_70';
file_modalshapes = 'flatplate_400_27_8_ballast_1_offset_5_nelem_70';
load(file_modalshapes)

% flow properties
rho = 1.1; % air density, kg/m^3

% strips properties
ns = 10; % number of strips
dy = modalstr.Lw/ns;
vec_y = (1:2:2*ns-1)*dy/2;
vec_dy = zeros(1,ns) + dy;
vec_b = zeros(1,ns) + modalstr.cw/2; % strips' semi-chords
vec_a = zeros(1,ns) + 0; % strips' EA positions
vec_cs = zeros(1,ns) + 1; % strips' flap positions (no flap, cs = 1)

% generation of structural, modal matrices
% - 3D matrix, the third dimension is the position of the strip
nm = 4; % number of modes to be used
Phi = zeros(3,nm,ns);
for kk=1:nm
    Fh_kk = griddedInterpolant(modalstr.y,-(modalstr.shape(3:6:end,kk))','linear'); % h(yjj) = -z
    h_interp_kk = Fh_kk(vec_y);
    Phi(1,kk,:) = h_interp_kk;

    Ftheta_kk = griddedInterpolant(modalstr.y,(modalstr.shape(5:6:end,kk))','linear'); % theta(yjj)
    theta_interp_kk = Ftheta_kk(vec_y);
    Phi(2,kk,:) = theta_interp_kk;

    Phi(3,kk,:) = zeros(1,ns); % beta(yjj) = 0 (no flaps)
end

% calculation of numerically exact flutter speed
% options = optimset('Display','iter');
% [flutt.U,flutt.res] = fsolve(@maxeig,15,options,rho,modalstr,nm)

options = optimset('Display','iter');
[flutt.U,flutt.res] = fsolve(@maxeig,15,options,rho,modalstr,nm,ns,Phi,vec_dy,vec_b,vec_a,vec_cs)


% stability check for a range of velocities
U = sort([1:0.5:20,flutt.U]); % a vector of flow velocities, m/s
nU = length(U);

% allocating space for aerodynamic and modal matrices
maed.A1 = zeros(3,3,ns);
maed.A2 = zeros(3,3,ns);
maed.A3 = zeros(3,3,ns);
maed.A4 = zeros(3,2,ns);
maed.A1_nc = maed.A1;
maed.A2_nc = maed.A2;
maed.A3_nc = maed.A3;
maed.B1 = zeros(2,3,ns);
maed.B2 = zeros(2,3,ns);
maed.B3 = zeros(2,3,ns);
maed.B4 = zeros(2,2,ns);

stab.A_anl = zeros(2*nm+2*ns);
stab.A_num = zeros(2*nm+2*ns);


for ii = 1:nU
    for jj = 1:ns
        % generation of aerodynamic matrices, fields of structure maed
        % - 3D matrices, the third dimension is the position of the strip
        amat = potential_2D_coefs(rho,U(ii),vec_b(jj),vec_a(jj),vec_cs(jj));
        maed(ii).A1(:,:,jj) = amat.A1;
        maed(ii).A2(:,:,jj) = amat.A2;
        maed(ii).A3(:,:,jj) = amat.A3;
        maed(ii).A4(:,:,jj) = amat.A4;
        maed(ii).A1_nc(:,:,jj) = amat.A1_nc;
        maed(ii).A2_nc(:,:,jj) = amat.A2_nc;
        maed(ii).A3_nc(:,:,jj) = amat.A3_nc;
        maed(ii).B1(:,:,jj) = amat.B1;
        maed(ii).B2(:,:,jj) = amat.B2;
        maed(ii).B3(:,:,jj) = amat.B3;
        maed(ii).B4(:,:,jj) = amat.B4;
    end



    % generation of all matrices involved on the final modal problem
    % - aerodynamics, structural dynamics, and step alpha disturbance q
    % - final system of equation:
    %
    %     - modal amplitudes:
    %  m_eta_s_etadd*(etadd) + 
    %  m_eta_s_etad*(etad) + 
    %  m_eta_s_eta*(eta) = m_eta_a_etadd*(etadd) +
    %                      m_eta_a_etad*(etad) +
    %                      m_eta_a_eta*(eta) +
    %                      m_eta_a_lbda*(lbda) + 
    %                      m_eta_q*(q)
    %  
    %     - lag states:
    %  lbdad = m_lbda_etadd*(etadd) + 
    %          m_lbda_etad*(etad) +
    %          m_lbda_eta*(eta) +
    %          m_lbda_lbda*(lbda) +
    %          m_lbda_q*(q)
    %
    % >> matrices involved in the modal eom are fields of the structure mmod

    mmod(ii).m_eta_s_etadd = diag(modalstr.mu(1:nm));                                   % structural modal mass matrix
    mmod(ii).m_eta_s_etad = 2*zeros(nm,nm)*diag(modalstr.omega(1:nm));                  % structural damping matrix
    mmod(ii).m_eta_s_eta = (diag(modalstr.omega(1:nm)))^2;                              % structural stiffness matrix
    mmod(ii).m_eta_a_etadd = zeros(nm);
    mmod(ii).m_eta_a_etad = zeros(nm);
    mmod(ii).m_eta_a_eta = zeros(nm);
    mmod(ii).m_eta_a_lbda = zeros(nm,2*ns);
    mmod(ii).m_lbda_etadd = zeros(2*ns,nm);
    mmod(ii).m_lbda_etad = zeros(2*ns,nm);
    mmod(ii).m_lbda_eta = zeros(2*ns,nm);
    mmod(ii).m_lbda_lbda = zeros(2*ns,2*ns);
    mmod(ii).m_eta_q = zeros(nm,1);
    mmod(ii).m_lbda_q = zeros(2*ns,1);
    
    
    for jj = 1:ns
        mmod(ii).m_eta_a_etadd = mmod(ii).m_eta_a_etadd + ...
                        Phi(:,:,jj).'*maed(ii).A1(:,:,jj)*Phi(:,:,jj)*vec_dy(jj);
        mmod(ii).m_eta_a_etad = mmod(ii).m_eta_a_etad + ...
                       Phi(:,:,jj).'*maed(ii).A2(:,:,jj)*Phi(:,:,jj)*vec_dy(jj);                  
        mmod(ii).m_eta_a_eta = mmod(ii).m_eta_a_eta + ...
                      Phi(:,:,jj).'*maed(ii).A3(:,:,jj)*Phi(:,:,jj)*vec_dy(jj);
        mmod(ii).m_eta_a_lbda(:,2*(jj-1)+1:2*(jj-1)+2) = Phi(:,:,jj).'*maed(ii).A4(:,:,jj)*vec_dy(jj);

        mmod(ii).m_lbda_etadd(2*(jj-1)+1:2*(jj-1)+2,1:nm) = maed(ii).B1(:,:,jj)*Phi(:,:,jj);
        mmod(ii).m_lbda_etad(2*(jj-1)+1:2*(jj-1)+2,1:nm) = maed(ii).B2(:,:,jj)*Phi(:,:,jj);
        mmod(ii).m_lbda_eta(2*(jj-1)+1:2*(jj-1)+2,1:nm) = maed(ii).B3(:,:,jj)*Phi(:,:,jj);
        mmod(ii).m_lbda_lbda(2*(jj-1)+1:2*(jj-1)+2,2*(jj-1)+1:2*(jj-1)+2) = maed(ii).B4(:,:,jj);

        mmod(ii).m_eta_q = mmod(ii).m_eta_q + ...
                  Phi(:,:,jj).'*maed(ii).A3(:,:,jj)*[0;1;0]*vec_dy(jj);
        mmod(ii).m_lbda_q(2*(jj-1)+1:2*(jj-1)+2,1) = maed(ii).B3(:,:,jj)*[0;1;0];
    end



    % building linear state-space matrix A (analitically)
    stab(ii).A_anl = [zeros(nm),eye(nm),zeros(nm,2*ns);
             (mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\(mmod(ii).m_eta_a_eta - mmod(ii).m_eta_s_eta),...
             (mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\(mmod(ii).m_eta_a_etad - mmod(ii).m_eta_s_etad),...
             (mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\mmod(ii).m_eta_a_lbda;
             mmod(ii).m_lbda_eta + mmod(ii).m_lbda_etadd*((mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\(mmod(ii).m_eta_a_eta - mmod(ii).m_eta_s_eta)), ...
             mmod(ii).m_lbda_etad + mmod(ii).m_lbda_etadd*((mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\(mmod(ii).m_eta_a_etad - mmod(ii).m_eta_s_etad)),...
             mmod(ii).m_lbda_lbda + mmod(ii).m_lbda_etadd*((mmod(ii).m_eta_s_etadd - mmod(ii).m_eta_a_etadd)\mmod(ii).m_eta_a_lbda)];


    % linearising numerically
    stab(ii).A_num = zeros(2*nm+2*ns);
    dX = 1e-8;
    for kk = 1:(2*nm+2*ns)
        X0 = zeros(2*nm+2*ns,1);
        X0(kk) = X0(kk)  + dX;
        dXdt_plus = dynamics_lin(0,X0,0,nm,ns,mmod(ii));

        X0 = zeros(2*nm+2*ns,1);
        X0(kk) = X0(kk) - dX;
        dXdt_minus = dynamics_lin(0,X0,0,nm,ns,mmod(ii));

        stab(ii).A_num(:,kk) = (dXdt_plus - dXdt_minus)/(2*dX);
    end
end




%% saving output file


file_out = [file_modalshapes,'_aeroelast_ns_',num2str(ns),'_nU_',num2str(nU)];
save(file_out,'file_modalshapes',...
              'maed',...
              'mmod',...
              'nm',...
              'ns',...
              'nU',...
              'Phi',...
              'rho',...
              'stab',...
              'U',...
              'vec_a',...
              'vec_b',...
              'vec_cs',...
              'vec_dy',...
              'vec_y',...
              'flutt');
          
          
%% verification of modal shapes interpolation

Phi_full = zeros(3,nm,length(modalstr.y)); 
for jj = 1:length(modalstr.y)
    Phi_full(:,:,jj) = [-modalstr.shape(3+6*(jj-1),1:nm);    % h(yjj) = -z
                    modalstr.shape(5+6*(jj-1),1:nm);    % theta(yjj)
                    zeros(1,nm)];                       % beta(yjj) = 0 (no flaps)
end

figure
subplot(2,4,1)
plot(modalstr.y,squeeze(Phi_full(1,1,:)),'b',vec_y,squeeze(Phi(1,1,:)),'rs')
subplot(2,4,2)
plot(modalstr.y,squeeze(Phi_full(2,1,:)),'b',vec_y,squeeze(Phi(2,1,:)),'rs')
subplot(2,4,3)
plot(modalstr.y,squeeze(Phi_full(1,2,:)),'b',vec_y,squeeze(Phi(1,2,:)),'rs')
subplot(2,4,4)
plot(modalstr.y,squeeze(Phi_full(2,2,:)),'b',vec_y,squeeze(Phi(2,2,:)),'rs')
subplot(2,4,5)
plot(modalstr.y,squeeze(Phi_full(1,3,:)),'b',vec_y,squeeze(Phi(1,3,:)),'rs')
subplot(2,4,6)
plot(modalstr.y,squeeze(Phi_full(2,3,:)),'b',vec_y,squeeze(Phi(2,3,:)),'rs')
subplot(2,4,7)
plot(modalstr.y,squeeze(Phi_full(1,4,:)),'b',vec_y,squeeze(Phi(1,4,:)),'rs')
subplot(2,4,8)
plot(modalstr.y,squeeze(Phi_full(2,4,:)),'b',vec_y,squeeze(Phi(2,4,:)),'rs')

