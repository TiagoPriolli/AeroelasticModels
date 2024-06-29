clear

load FullAppStrips Phi_App
load WingStrips Phi_Wg

n_mode_wg = 1; %[1 1],[2 4],[3 5]
n_mode = 5;

hvec = permute(Phi_App(1,n_mode,:),[3,1,2]);
alphavec = permute(Phi_App(2,n_mode,:),[3,1,2]);

hvec_wg = permute(Phi_Wg(1,n_mode_wg,:),[3,1,2]);
alphavec_wg = permute(Phi_Wg(2,n_mode_wg,:),[3,1,2]);

test01 = dot(hvec,hvec_wg)
test02 = dot(alphavec,alphavec_wg)

