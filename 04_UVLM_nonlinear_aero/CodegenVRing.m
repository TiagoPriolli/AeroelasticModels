function CodegenVRing(M,N)

%[uRing,vRing,wRing] = voringMatrix_FlatPlate(xyz,M,N,circulation,...
%                                                    xgrid,ygrid,zgrid)
cfg = coder.config('mex');
cfg.GenerateReport = true;

ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0,[3 1]);               % xyz
ARGS{1}{2} = coder.Constant(M);                   % M
ARGS{1}{3} = coder.Constant(N);                   % N
ARGS{1}{4} = coder.typeof(0,[M N]);               % circulation (Bound)
ARGS{1}{5} = coder.typeof(0,[(M+1) (N+1)]);       % xgrid (Bound)
ARGS{1}{6} = coder.typeof(0,[(M+1) (N+1)]);       % ygrid (Bound)
ARGS{1}{7} = coder.typeof(0,[(M+1) (N+1)]);       % zgrid (Bound)

% Change to CodeGen Folder
cd('04_UVLM_nonlinear_aero');

% Generate Code
disp('Generating Voring Code...')
codegen -config cfg voringMatrix_FlatPlate -args ARGS{1} -o voringMatrix_FPbound

% Return to Project Root
cd('..');
