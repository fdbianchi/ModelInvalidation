function [Tg, Tg0] = stoeplitz(G,N,debug)

% STOEPLITZ produces the Toeplitz matrix corresponding to the convolution
% kernel of an LTI system
%
% Use:  
%   T_G = STOEPLITZ(G,N)
%
%   G: LTI system in lti or matrix form
%   N: final discrete-time
%

% fbianchi - 10/04/2018  


if (nargin < 3)
    debug = 0; % 0 to avoid debug
end

% system matrices
if isa(G,'lti')
    [A,B,C,D] = ssdata(G);
    [ny,nu]   = size(G);
    ns = order(G);
elseif isnumeric(G) && G(end,end) == -Inf
    [A,B,C,D]  = ltiss(G);
    [ns,nu,ny] = sinfo(G);
elseif isnumeric(G)
    ns = 0;
    [ny,nu] = size(G);
    A = zeros(ns);    B = zeros(ns,nu);
    C = zeros(ny,ns); D = G;
else
    error('G is not a system')
end

% memory pre-allocation 
Tg  = zeros(N*ny,N*nu);
Tg0 = zeros(N*ny,ns);

% 1st row
Tg(1:ny,1:nu)  = D;   auxTg  = B;    % T_G
Tg0(1:ny,1:ns) = C;   auxTg0 = A; 

% 2nd row
if (N > 1)
    Tg((1:ny)+ny,1:2*nu) = [C*auxTg D];
    auxTg = [A*auxTg B];

    Tg0((1:ny)+ny,1:ns)  = C*auxTg0; 
    auxTg0 = A*auxTg0; 
end

% i-th rows
if (N > 2)
    for ii=3:N
        Tg((1:ny)+ny*(ii-1),1:ii*nu) = [C*auxTg D];
        auxTg = [A*auxTg B];

        Tg0((1:ny)+ny*(ii-1),1:ns)   = C*auxTg0; 
        auxTg0 = A*auxTg0; 
    end
end


% =======================================================================
% for debug
%
if debug

    G = ss(A,B,C,D,1);
    
    x0 = rand(ns,1);
    
    % time
    t = (0:N-1);

    % pulse signal
    u = zeros(nu,N);
    u(nu,1:round(N/10)) = 1;

    % outputs
    y1  = lsim(G,u,t,x0);
    Tu  = mtoeplitz(u);
    Tx0 = repmat(x0,N,1);
    y2u = Tg*Tu;
    y20 = Tg0*x0;
    y2  = y20 + y2u(:,1); y2 = reshape(y2,ny,N)';
    
    % plots
    subplot(211)
    plot(t,u);
    xlabel('samples'); ylabel('u');
    
    subplot(212)
    hl = plot(t,y1,'r',t,y2,'b.');
    xlabel('samples'); ylabel('y');
    legend(hl([1 ny+1]),'lsim','toeplitz')

end
