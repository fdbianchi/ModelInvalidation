
% Model invalidation:
%
% Simple LTI example to check tool capabilities

% fbianchi - 13/07/2018

clearvars, close all, % clc

N  = 50;     % number of samples
Ts = 0.01;   % sampling time

% -----------------------------------------------------------------------
% Plant + model

Gc = ss(diag([-1,-5,-10]),[1;1;1],[1 2 1],0);
Gd = c2d(Gc,Ts);
ns = order(Gd);      % number of states
[no,ni] = size(Gd);  % number of outputs, inputs

% uncertainty weight
% Wd  = 0.05*eye(no);
% Wdc = tf(0.1*[1 1],[0.01 1])*eye(no);
Wdc = makeweight(0.1,1,10);
Wd = c2d(Wdc,Ts);

% real system
Delta = usample(ultidyn('Delta',[no ni]));
Grc = Gc*(1+Wdc*Delta);
Gr  = c2d(Grc,Ts);
% bodemag(Gc,Grc)

% -----------------------------------------------------------------------
% input-output data

tk = (0:N-1)*Ts;
% input -> pulse
pw = 10;
uk = zeros(ni,N);
for ii = 1:ni
    uk(ii,2*ii*pw+(1:pw)) = ones(1,pw);
end

% initial conditions
x01 = zeros(order(Gr),1);   % zero initial conditions
x02 = rand(order(Gr),1);    % random initial conditions

% output
yk1 = lsim(Gr,uk,tk,x01)';  % w/ zero initial conditions
yk2 = lsim(Gr,uk,tk,x02)';  % w/ random initial conditions

% with norm-2 noise
opt.norm   = 2;
nk = rand(size(yk1)); 
if opt.norm == inf
    yk1 = yk1 + nk*0.01;
    yk2 = yk2 + nk*0.01;
else
    norm_nk = norm(nk);%/length(nk);
    yk1 = yk1 + nk*0.01/norm_nk;
    yk2 = yk2 + nk*0.01/norm_nk;
end

opt.verb   = 1;             % result details ON
opt.debug  = 1;             % check signal reconstruction ON
opt.solver = 'sedumi';      % solver for LMIs
opt.obj    = [0  0.01];     % obj: minimize uncertainty bound

% model invalidation w/zero initial conditions
figure
[dmx,wmx,w,d] = modinval(Gd,Wd,uk,yk1,[],opt);
title('Without initial conditions')

% model invalidation w/random initial conditions ==> infeasible
opt.is_xO = false;     % ignoring it when solving
[dmx,wmx,w,d] = modinval(Gd,Wd,uk,yk2,[],opt);


opt.is_xO = true;      % considering it when solving
figure
[dmx,wmx,w,d] = modinval(Gd,Wd,uk,yk2,[],opt);
title('With initial conditions')
