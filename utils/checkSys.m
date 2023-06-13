function sysInfo = checkSys(sys)

% CHECKSYS serves to check the system description

% fbianchi - 2022-06-12

if isa(sys,'p_ss')
    % lpv objects
    sysInfo.sys = sys;
    sysInfo.typ = 'lpv';
    if isa(sys,'pass')
        sysInfo.styp = 'aff';
    elseif isa(sys,'ppss')
        sysInfo.styp = 'pol';
    else
        sysInfo.styp = 'gral';
    end
    [ny,nu,ns] = size(sys);
    sysInfo.ns = ns;
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sysInfo.nv = size(sys,3);
    
elseif (isa(sys,'ss') || isa(sys,'tf') || isa(sys,'zpk'))
    % ss/tf/zpk systems
    sysInfo.typ = 'lti';
    sysInfo.ns = order(sys);
    [ny,nu] = iosize(sys);
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    
elseif ispsys(sys)
    % psys systems
    sysInfo.typ = 'lpv';
    [type,nv,ns,nu,ny] = psinfo(sys);
    if strcmp(type,'aff')
        sysInfo.styp = 'aff';
    else
        sysInfo.styp = 'pol';
    end
    sysInfo.ns = ns;
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sysInfo.nv = nv;
    
elseif isnumeric(sys)
    % lmitool/robust (old) objects
    sysInfo.typ = 'mat';
    sys = mat2lti(sys);
    sysInfo.ns = order(sys);
    [nu,ny] = iosize(sys);
    sysInfo.ny = ny;
    sysInfo.nu = nu;
    sysInfo.nv = 1;
else
    error('STANDARDIZESYS:InputError',...
        'SYS is not a valid system description')
end
