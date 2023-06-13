function [Tg,Tg0] = ltoeplitz(G,pk)

% LTOEPLITZ produces the Toeplitz matrix corresponding to the convolution
% kernel of an LPV affine system
%
% Use:  
%   T_G = STOEPLITZ(G,pk)
%
%   G: affine LPV system (LMI)
%   pk: varying parameter vector
%

% fbianchi - 10/04/2018  


debug = 0; % 0 to avoid debug

% number of samples
N = size(pk,2);

% system info
if isa(G,'p_ss')
    if ~isa(G,'pass')
        error('G must be an affine LPV model')
    end
    [ny,nu,ns] = size(G);
    newLPV = true;
elseif ispsys(G)
    [typ,~,ns,nu,ny] = psinfo(G);
    if ~strcmp(typ,'aff')
        error('G must be an affine LPV model')
    end
    newLPV = false;
else
    error('SYS is not a valid LPV system description')
end

% memory pre-allocation 
Tg(N*ny,N*nu) = 0;
Tg0(N*ny,ns)  = 0;

% 1st row
if newLPV
    [A,B,C,D] = ssdata(ss(G,pk(:,1)));
else
    [A,B,C,D] = ltiss(psinfo(G,'eval',pk(:,1)));
end
Tg(1:nu,1:ny)  = D;   auxTG  = B;   % T_G
Tg0(1:ny,1:ns) = C;   auxTg0 = A; 

% 2nd row
if (N > 1)
    if newLPV
        [A,B,C,D] = ssdata(ss(G,pk(:,2)));
    else
        [A,B,C,D] = ltiss(psinfo(G,'eval',pk(:,2)));
    end
    Tg((1:ny)+ny,1:2*nu) = [C*auxTG D]; auxTG = [A*auxTG B];

    Tg0((1:ny)+ny,1:ns)  = C*auxTg0; 
    auxTg0 = A*auxTg0; 
end

% i-th rows
if (N > 2)
    for ii=3:N
        if newLPV
            [A,B,C,D] = ssdata(ss(G,pk(:,ii)));
        else
            [A,B,C,D] = ltiss(psinfo(G,'eval',pk(:,ii)));
        end
        Tg((1:ny)+ny*(ii-1),1:ii*nu) = [C*auxTG D];
        auxTG = [A*auxTG B];

        Tg0((1:ny)+ny*(ii-1),1:ns)   = C*auxTg0; 
        auxTg0 = A*auxTg0; 
    end
end


% =======================================================================
% for debug
%
if debug

    x0 = rand(ns,1);

    % time
    t = 0:N-1;

    % step signal
    u = ones(nu,1)*(t <= round(N/10));

    % outputs
    x  = zeros(ns,1);
    y(ny,N) = 0;
    for ii=1:length(pk)
        if newLPV
            [A,B,C,D] = ssdata(ss(G,pk(:,ii)));
        else
            [A,B,C,D] = ltiss(psinfo(G,'eval',pk(:,ii)));
        end
        y(:,ii) = C*x + D*u(:,ii);
        x = A*x + B*u(:,ii);
    end
%     y   = pdsimuld(G,u,pk,x0);
    
    Tu  = mtoeplitz(u);
    ytu = Tg*Tu;
    yt0 = Tg0*x0;
    yt  = yt0 + ytu(:,1); 
    yt  = reshape(yt,ny,N);
    
    % plots
    subplot(211)
    plot(t,u);
    xlabel('samples'); ylabel('u');
    
    subplot(212)
    plot(t,y,'r',t,yt,'b.')
    xlabel('samples'); ylabel('y');
    legend('lsim','toeplitz')

end




