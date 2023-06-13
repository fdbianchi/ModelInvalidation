
% Model invalidation:
%
% Simple LPV example to check tools capabilities

% fbianchi - 15/07/2018

clearvars, close all, % clc

N  = 100;    % number of samples
Ts = 1;

% -----------------------------------------------------------------------
% lmitool functions

% LPV model
A0 = [0 1; -0.3  -0.6];
A1 = [0 0;  0    -0.6];
A2 = [0 0; -0.01  0];
B  = eye(2); 
C  = eye(2);
D  = zeros(2);
s0 = ltisys(eye(2)+Ts*A0,Ts*B,C,D);
s1 = ltisys(Ts*A1,zeros(2),zeros(2),zeros(2),0);
s2 = ltisys(Ts*A2,zeros(2),zeros(2),zeros(2),0);

p1mn = -1; p1mx = 1;
p2mn = -1; p2mx = 1;
pv = pvec('box',[p1mn p1mx;p2mn p2mx]);
saff = psys(pv,[s0 s1 s2]);

[~,nv,ns,ni,no] = psinfo(saff); 

% uncertainty weight
% Wd  = 1*eye(no);
% Wdc = tf(0.1*[1 1],[0.01 1])*eye(no);
Wdc = makeweight(0.1,1,10)*eye(no);
Wd  = c2d(Wdc,Ts);

% real system
Delta  = usample(ultidyn('Delta',[no ni]));
[ad,bd,cd,dc] = ssdata(Delta);
nDelta = norm(Delta,inf);

sim('LPV_mimo_ex_mdl.slx');

% invalidation data
yk = data.y';
uk = data.u(:,1:2)';
pk = data.u(:,3:4)';

% invalidation
opt.norm   = inf;
opt.N      = N;
opt.verb   = 1;
opt.debug  = 1;
opt.solver = 'sedumi';
% opt.solver = 'sdpt3';
opt.is_xO  = 0;%any(x0);
opt.obj    = [0  0.01];
opt.obj    = [1 0];
% opt.obj    = [0.5 0.5];
[dmx,wmx,w,d] = modinval(saff,Wd,uk,yk,pk,opt);

fprintf('\n ||Delta||_inf = %5.4f\n',nDelta)


% -----------------------------------------------------------------------
% using LPVtools
if (exist('p_ss','class') == 8)

    % LPV model
    % - matrices
    Ad(:,:,1) = eye(2)+Ts*A0;
    Ad(:,:,2) = Ts*A1;
    Ad(:,:,3) = Ts*A2;
    Bd  = Ts*B; 
    % - parameter set
    p1mn = -1; p1mx = 1;
    p2mn = -1; p2mx = 1;
    pv = pset.Box([p1mn p1mx;p2mn p2mx]);
    pdG = pass(Ad,Bd,C,D,pv);

    [no,ni,ns,~,nv] = size(pdG); 

    % uncertainty weight
    Wdc = makeweight(0.1,1,10)*eye(no);
    Wd  = c2d(Wdc,Ts);

    % invalidation
    opt.norm   = inf;
    opt.N      = N;
    opt.verb   = 1;
    opt.debug  = 1;
    opt.solver = 'sedumi';
    % opt.solver = 'sdpt3';
    opt.is_xO  = 0;%any(x0);
    opt.obj    = [0  0.01];
    opt.obj    = [1 0];
    % opt.obj    = [0.5 0.5];
    figure
    [dmx,wmx,w,d] = modinval(pdG,Wd,uk,yk,pk,opt);

    fprintf('\n ||Delta||_inf = %5.4f\n',nDelta)

end
