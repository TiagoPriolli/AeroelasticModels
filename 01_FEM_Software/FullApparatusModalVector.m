clear all
% loading modal shapes
%file_modalshapes = 'flatplate_350_40_8_ballast_1_offset_15_nelem_70';
file_modalshapes = 'FullApparatus';
load(file_modalshapes)

% strips properties
ns = 20; % number of strips
dy = modalstr.wing.Lw/ns;
vec_y = (1:2:2*ns-1)*dy/2;
vec_dy = zeros(1,ns) + dy;
vec_b = zeros(1,ns) + modalstr.wing.cw/2; % strips' semi-chords
vec_a = zeros(1,ns) + 0; % strips' EA positions
vec_cs = zeros(1,ns) + 1; % strips' flap positions (no flap, cs = 1)

% generation of structural, modal matrices
% - 3D matrix, the third dimension is the position of the strip
nm = 5; % number of modes to be used
Phi_App = zeros(3,nm,ns);

for kk=1:nm
    zshape = modalstr.shape(3:6:end,kk);
    Fh_kk = griddedInterpolant(modalstr.y(1,end-70:end)-modalstr.y(1,end-70),-(zshape(end-70:end,1))','linear'); % h(yjj) = -z
    h_interp_kk = Fh_kk(vec_y);
    Phi_App(1,kk,:) = h_interp_kk;
    
    yshape = modalstr.shape(5:6:end,kk);
    Ftheta_kk = griddedInterpolant(modalstr.y(1,end-70:end)-modalstr.y(1,end-70),(yshape(end-70:end))','linear'); % theta(yjj)
    theta_interp_kk = Ftheta_kk(vec_y);
    Phi_App(2,kk,:) = theta_interp_kk;

    Phi_App(3,kk,:) = zeros(1,ns); % beta(yjj) = 0 (no flaps)
end

save FullAppStrips Phi_App