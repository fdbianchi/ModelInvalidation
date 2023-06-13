
function [dmx,wmx,w,diagnostic] = modinval(sys,Wd,uk,yk,pk,opt,x0)

% 
% MODINVAL implements the (In)validation of LPV or LTI model with LTI 
% uncertainty
%
% Use:
%   [dmx,wmx,w] = MODINVAL(sys,Wd,uk,yk,pk,opt)
%
% Input arguments:
% - sys: LPV model (in affine form) or LTI (lmitool objects or constrol
%                   systems toolbox objectes)
% - Wd: Diagonal matrix with the uncertainty weight (additive or
%       multiplicative uncertainty)
% - uk, yk: input-output datum
% - pk: varying parameter datum (use [] in LTI case)
% - opt: struct with settings
%     opt.obj    : a 2x1 vector for setting the optimisation objetive
%               [wmx^2 0]     = minimisation of the noise norm (default)
%               [0     dmx^2] = minimisation of the perturbation norm
%               [a     b]     = minimisation of the sum of previous bounds
%                               with a and b weights
%     opt.norm   : norm for noise set (2=2-norm, inf=inf-norm)
%     opt.verb   : optimization information (1=verbose)
%     opt.solver : solver ('sedumi', 'sdpt3')
%     opt.N      : last sample index
%     opt.n      : first sample index
%     opt.is_xO  : 1 = compute initial conditions
%     opt.debug  : 1 = debug mode 


% fbianchi - 2008-03-13
% fbianchi - 2018-03-29 - rev 1.0
% fbianchi - 2023-06-12 - rev 1.1


% Default output argument
dmx = []; wmx = []; w = [];

% -------------------------------------------------------------------------
% argument tests
sysInfo = checkSys(sys);
if strcmp(sysInfo.typ,'lpv') && (nargin < 5)
    disp('use: [dmx,wmx,w] = modinval(Sys,Wd,uk,yk,pk)');
    return,
elseif (nargin < 4)
    disp('use: [dmx,wmx,w] = modinval(Sys,Wd,uk,yk)');
    return,
end

if (nargin < 6)
    opt.obj    = [1 0];         % minimisation of the noise bound (default)
    opt.norm   = 2;             % norm for noise set
    opt.verb   = 0;             % verbose
    opt.solver = 'sedumi';      % solver
    opt.N      = length(uk);    % last sample index
    opt.n      = 1;             % first sample index
    opt.is_xO  = 1;             % compute initial conditions
    opt.debug  = 0;             % 1 = debug mode
else    
    if ~isstruct(opt)
        error('OPT must be a struct as explained in help')
    else
        if ~isfield(opt,'obj')
            opt.obj    = [1 0];         % minimisation of the noise bound (default)
        end
        if ~isfield(opt,'norm')
            opt.norm   = 2;             % norm for noise set
        end
        if ~isfield(opt,'verb')
            opt.verb   = 0;             % verbose
        end
        if ~isfield(opt,'solver')
            opt.solver = 'sedumi';      % solver
        end
        if ~isfield(opt,'N')
            opt.N      = length(uk);    % last sample index
        end
        if ~isfield(opt,'n')
            opt.n      = 1;             % first sample index
        end
        if ~isfield(opt,'is_xO')
            opt.is_xO  = 1;             % compute initial conditions
        end
        if ~isfield(opt,'debug')
            opt.debug  = 0;             % 1 = debug mode
        end
    end
end


% -------------------------------------------------------------------------
% LMI solver settings
switch opt.solver
    case 'sedumi'
        lmiopt = sdpsettings('solver',opt.solver,...'sedumi.eps',1e-5,...
            'verbose', opt.verb);
    case 'sdpt3'
        lmiopt = sdpsettings('solver',opt.solver,'verbose', opt.verb);
    otherwise
        error('Solver %s not soported',opt.solver);
end

% =========================================================================
% parameter set information

% number of samples
N = opt.N;
% First sample considered in the invalidation
n = opt.n;
% Dimensions of the model
ns = sysInfo.ns;
ny = sysInfo.ny;

% =========================================================================
% Toeplitz matrices

% Toeplitz matrix for the uncertainty weights 
T_Wd = stoeplitz(Wd,N-n+1);%N+n-1);
% Toeplitz matrix for the input
Tu   = mtoeplitz(uk(:,n:N+n-1,:));

% Construction of the model Toeplitz matrices
if strcmp(sysInfo.typ,'lpv')
    % LPV model case
    [T_G, T_G0] = ltoeplitz(sys,pk(:,n:N));
else
    % LTI model case
    [T_G, T_G0] = stoeplitz(sys,N-n+1);
end

% -------------------------------------------------------------------------
% Optimisation problem

% yalmip initialization
yalmip('clear')
lmis   = [];
eigtol = 1e-8;

% Objective
objlmi = 0;
if (opt.obj(1) > 0) && (opt.obj(2) == 0) 
    wmax = opt.obj(1)^2;     % Uncertainty bound
    dmax = sdpvar(1);        % Noise energy bound
    objlmi = objlmi + dmax;
    if (opt.norm == inf)
        objStr = sprintf('\tmin  ||dmax||_inf\n\t\ts.t. wmax <= %5.4f\n',opt.obj(1));
    else
        objStr = sprintf('\tmin  ||dmax||_2\n\t\ts.t. wmax <= %5.4f\n',opt.obj(1));
    end
    
elseif (opt.obj(1) == 0) && (opt.obj(2) > 0)
    wmax = sdpvar(1);        % Uncertainty bound
    objlmi = objlmi + wmax; 
    if (opt.norm == inf)
        dmax   = opt.obj(2); % Noise amplitud bound    
        objStr = sprintf('\tmin   wmax\n\t\ts.t. ||dmax||_inf <= %5.4f\n',opt.obj(2));
    else
%         dmax = ((N-n+1)*opt.obj(2))^2; % Noise energy bound
        dmax = (opt.obj(2))^2;
        objStr = sprintf('\tmin   wmax\n\t\ts.t. ||dmax||_2 <= %5.4f\n',opt.obj(2));
    end        
    
elseif (opt.obj(1) > 0) && (opt.obj(2) > 0)
    dmax = sdpvar(1);        % Noise energy bound
    wmax = sdpvar(1);        % Uncertainty bound
    objlmi = objlmi + opt.obj(1)*wmax + opt.obj(2)*dmax; 
    if (opt.norm == inf)
        objStr = sprintf('\tmin %3.2f||dmax||_inf + %3.2f wmax\n',opt.obj(1),opt.obj(2));
    else
        objStr = sprintf('\tmin %3.2f||dmax||_2 + %3.2f wmax\n',opt.obj(1),opt.obj(2));
    end

else
    error('Objective not valid')
end

% Disturbance sequence
wk = sdpvar(ny,N-n+1,'full'); 
if (opt.is_xO == 1)
    x0 = sdpvar(ns,1);   % initial conditions
    d0 = sdpvar(1);      % initial conditions bound
    ek = reshape(yk(:,n:N),(N-n+1)*ny,1) - T_G0*x0 - T_G*Tu(:,1);
    dk = ek - reshape(wk,(N-n+1)*ny,1);
else
    dk = reshape(yk(:,n:N),(N-n+1)*ny,1) - T_G*Tu(:,1) - ...
         reshape(wk,(N-n+1)*ny,1);
end

% Constraint: initial conditions
if (opt.is_xO == 1)
    e0 = ek(n:n+ny-1);
%     lmis = [lmis, [d0 e0';e0 eye(ny)] >= eigtol*eye(ny+1)];
    lmis = [lmis, -d0 <= e0 <= d0];
    objlmi = objlmi + d0; 
end

% Constraint: noise norm
if (opt.norm == inf)
    % - infty-norm
    lmis = [lmis, -dmax <= dk <= dmax ];
else
    % - 2-norm:
    lmis = [lmis, ...
        [dmax dk(:)';dk(:) eye((N-n+1)*ny)] >= eigtol*eye((N-n+1)*ny+1) ];
end

% Constraint: Uncertainty
T_Gd = T_Wd*T_G;
Tw   = mtoeplitz(wk);    % Noise Toeplitz
M11  = (T_Gd*Tu)'*T_Gd*Tu;
Mat  = [M11, Tw'; Tw, eye(ny*(N-n+1))*wmax];
lmis = [lmis, 0.5*(Mat+Mat') >= 0*eigtol*eye(size(Mat))]; % <= if D=0, 
                                                          %    T_G is singular

if ~isa(objlmi,'sdpvar'),
    objlmi = []; 
end

% Solving LMIs
diagnostic = optimize(lmis,objlmi,lmiopt);

% -------------------------------------------------------------------------
% Results

if ~((diagnostic.problem == 0) || (diagnostic.problem == 4))

    if opt.verb
        fprintf('\n-------------------------')
        fprintf('\nThe problem could not be solved:')
        fprintf('\n')
        fprintf('\n\tSolver says: %s\n', diagnostic.info)
        fprintf('\n-------------------------')
        fprintf('\n')
    end
    
else
    w   = value(wk);
    d   = value(dk);
    wmx = sqrt(value(wmax));
    if (opt.norm == inf)
        % - infty-norm
        dmx = value(dmax);
        str_dmax = sprintf('(real ||dk||_inf = %5.4f)',max(max(abs(d))));
    else
        % - 2-norm:
        dmx = sqrt(value(dmax));%/(N-n+1);
        str_dmax = sprintf('(real ||dk||_2 = %5.4f)',sqrt(d(:)'*d(:))/(N-n+1));
    end
    
    if opt.verb
        fprintf('\n-------------------------')
        fprintf('\nInvalidation results:')
        fprintf('\n\tOpmitization: \n\t%s', objStr)
        fprintf('\n\tSolver says: %s\n', diagnostic.info)
        fprintf('\n\twmx = %5.4f',wmx)
        fprintf('\n\tdmx = %5.4f %s',dmx,str_dmax)
        fprintf('\n-------------------------')
        fprintf('\n')
    end
    
    % ====================================================================
    if opt.debug
        
        tk = 1:N-n+1;
        
        yt = T_G*Tu(:,1) + value(dk);
        if (opt.is_xO == 1)
            yt = yt + T_G0*value(x0);
        end
        yt = reshape(yt,ny,N-n+1);
        
        set(gcf,'Position',[680 125 764 853]);
        
        if ispsys(sys)
            nplots = 5;
        else
            nplots = 4;
        end
        
        m = 1; subplot(nplots,1,m)
        stairs(tk,uk(:,n:N)');
        xlabel('samples'); ylabel('uk');
        
        m = m + 1; subplot(nplots,1,m)
        [tk1,yk1] = stairs(tk,yk(:,n:N)');
        [tk2,yk2] = stairs(tk,yt');
        hl = plot(tk1,yk1,'r',tk2,yk2,'b');
        xlabel('samples'); ylabel('yk');
        legend(hl([1 ny+1]),'data','toeplitz')
        
        if ispsys(sys)
            m = m + 1; subplot(nplots,1,m)
            stairs(tk,pk(:,n:N)');
            xlabel('samples'); ylabel('pk');
        end
        
        m = m + 1; subplot(nplots,1,m)
        stairs(tk,reshape(value(dk),ny,N-n+1)');
        xlabel('samples'); ylabel('dk');
        
        m = m + 1; subplot(nplots,1,m)
        stairs(tk,value(wk)');
        xlabel('samples'); ylabel('wk');
        
    end
end




