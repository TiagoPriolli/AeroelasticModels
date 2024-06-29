function f = maxeig(U,rho,modalstr,nm)

nU = length(U);

% wing strips' properties
ns = length(modalstr.y)-1; % number of strips
vec_y = (modalstr.y(1:end-1) + modalstr.y(2:end))/2; % strips' y position yi
vec_dy = modalstr.y(2:end) - modalstr.y(1:end-1); % strips' widths
vec_b = zeros(1,ns) + modalstr.cw/2; % strips' semi-chords
vec_a = zeros(1,ns) + 0; % strips' EA positions
vec_cs = 0*vec_y + 1; % strips' flap positions (no flap, cs = 1)

% generation of structural, modal matrices
% - 3D matrix, the third dimension is the position of the strip
Phi = zeros(3,nm,ns); 
for jj = 1:ns
    Phi(:,:,jj) = [-modalstr.shape(3+6*(jj-1),1:nm);    % h(yjj) = -z
                    modalstr.shape(5+6*(jj-1),1:nm);    % theta(yjj)
                    zeros(1,nm)];                       % beta(yjj) = 0 (no flaps)
end

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

end

f = max(real(eig(stab.A_anl)));
if U <=2 % knowing that the flutter velocity for the plate is certainly greater than 2
    f = -100;
end