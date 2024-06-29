clear all
% loading modal shapes
%file_modalshapes = 'flatplate_350_40_8_ballast_1_offset_15_nelem_70';
file_modalshapes = 'Wing';
load(file_modalshapes)

% strips properties
ns = 20; % number of strips
dy = modalstr.Lw/ns;
vec_y = (1:2:2*ns-1)*dy/2;
vec_dy = zeros(1,ns) + dy;
vec_b = zeros(1,ns) + modalstr.cw/2; % strips' semi-chords
vec_a = zeros(1,ns) + 0; % strips' EA positions
vec_cs = zeros(1,ns) + 1; % strips' flap positions (no flap, cs = 1)

% generation of structural, modal matrices
% - 3D matrix, the third dimension is the position of the strip
nm = 5; % number of modes to be used
Phi_Wg = zeros(3,nm,ns);

for kk=1:nm
    Fh_kk = griddedInterpolant(modalstr.y,-(modalstr.shape(3:6:end,kk))','linear'); % h(yjj) = -z
    h_interp_kk = Fh_kk(vec_y);
    Phi_Wg(1,kk,:) = h_interp_kk;
    
    Ftheta_kk = griddedInterpolant(modalstr.y,(modalstr.shape(5:6:end,kk))','linear'); % theta(yjj)
    theta_interp_kk = Ftheta_kk(vec_y);
    Phi_Wg(2,kk,:) = theta_interp_kk;

    Phi_Wg(3,kk,:) = zeros(1,ns); % beta(yjj) = 0 (no flaps)
end

save WingStrips Phi_Wg