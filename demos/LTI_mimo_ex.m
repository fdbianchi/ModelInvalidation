
% Model invalidation:
%
% Simple LTI example to check tools capabilities

% fbianchi - 15/07/2018

clearvars, close all, % clc

N  = 80;    % number of samples
Ts = 0.01;

Gc = ss(diag([-1,-5,-10]),[1 1;0 1;1 0],[1 0 2;0 1 1],0);
Gd = c2d(Gc,Ts);
ns = order(Gd);      % number of states
[no,ni] = size(Gd);  % number of outputs, inputs

% uncertainty weight
% Wd  = 0.05*eye(no);
% Wdc = tf(0.1*[1 1],[0.01 1])*eye(no);
Wdc = makeweight(0.1,1,10)*eye(no);
Wd  = c2d(Wdc,Ts);

% real system
Delta = usample(ultidyn('Delta',[no ni]));
nDelta = norm(Delta,inf);
Grc = Gc*(1+Wdc*Delta);
Gr  = c2d(Grc,Ts);
sigma(Gc,Grc)

% data
tk = (0:N-1)*Ts;
pw = 10;
uk = zeros(ni,N);
for ii = 1:ni
    uk(ii,2*ii*pw+(1:pw)) = ones(1,pw);
end

% initial conditions
x0 = [rand(order(Gd),1); 0; 0; 0];
x0 = rand(order(Gr),1);
x0 = zeros(order(Gr),1);

% output
yk = lsim(Gr,uk,tk,x0)';

% with norm-2 noise
opt.norm  = inf;
norm_nk_d = 0.01;
rng(56398);
nk = rand(size(yk));
if opt.norm == inf
    nk = nk*norm_nk_d;
else
    norm_nk = norm(nk(:))/length(nk(:));
    nk = nk*norm_nk_d/norm_nk;
end
yk = yk + nk;

opt.verb   = 1;
opt.debug  = 1;
opt.solver = 'sedumi';
opt.is_xO  = any(x0);
opt.obj    = [0  norm_nk_d];
% opt.obj    = [1 0];
% opt.obj    = [0.5 0.5];
[dmx,wmx,w,d] = modinval(Gd,Wd,uk,yk,[],opt,x0(1:3,1));

fprintf('\n ||Delta||_inf = %5.4f\n',nDelta)

