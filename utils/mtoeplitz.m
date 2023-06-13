function Tu = mtoeplitz(U)

% MTOEPLITZ produces the Toeplitz matrix corresponding to signal vector uk
%
% Use:  
%   Tu = MTOEPLITZ(u)

% fbianchi - 10/04/2018  


[r,c] = size(U);

data = [zeros(r*(c-1),1);U(:)];                  % build vector of user data

cid = (0:r*c-1)';
rid = r*(c-1)+1:-r:1;

idTu = cid(:,ones(c,1)) + rid(ones(r*c,1),:);    % Toeplitz subscripts

Tu = data(idTu);                                 % actual data
