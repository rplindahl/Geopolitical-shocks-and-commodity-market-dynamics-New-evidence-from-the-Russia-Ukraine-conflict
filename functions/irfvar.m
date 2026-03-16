% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

% B: u = B*w

function [IRF]=irfvar(A,B,p,h)

q=size(B,1);
J=[eye(q,q) zeros(q,q*(p-1))];
IRF=reshape(J*A^0*J'*B,q^2,1);

for i=1:h
	IRF=([IRF reshape(J*A^i*J'*B,q^2,1)]);
end

