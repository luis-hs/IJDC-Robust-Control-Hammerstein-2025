function output = dsolvelmi_rayouf(A,B, minf, msup)
nx = size(A, 1);
nu = size(B, 2);
X = sdpvar(nx, nx, 'symmetric'); 
R = sdpvar(nu, nx, 'full');

LMIs = [];
% The conditions are in equation (26) of the article by (Rayouf et al.,2018)

LMIs = [LMIs; [X, (A*X - minf*B*R)'; (A*X - minf*B*R), X] >= 0;
              [X, (A*X - msup*B*R)'; (A*X - msup*B*R), X] >= 0 ];
solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
p = min(checkset(LMIs));
output.feas = 0;
if p > 0  
    output.X = double(X);
    output.R = double(R);
    output.feas = 1;
    output.minf = double(minf);
    output.msup = double(msup);
end
end