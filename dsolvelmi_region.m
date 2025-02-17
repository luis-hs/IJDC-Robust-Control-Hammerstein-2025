% Luís Henrique dos Santos 16/02/2025
% Implementation of the numerical example from the article A new approach
% for robust control based on parametric Hammerstein models

function output = dsolvelmi_region(A, B, minf, msup, ts, eta, v0)
nx = size(A, 1);
nu = size(B, 2);
P = sdpvar(nx, nx, 'symmetric');
H = P;
S = sdpvar(nu, nx, 'full');
B1 = minf*B;
B2 = msup*B;

LMIs = [];
%--------------------------------------------------------------------------
% Equations for the region of a circle centered at the origin | Limits ts
%--------------------------------------------------------------------------
sigma = 4/ts; % 2% criteria ts = 4/sigma
ro = exp(-sigma);
R11 = -(ro^2);
R12 = 0;
R22 = 1;

% I want MUP = 4.5%, so MU = 0.045 = e^(-pi*zeta/sqrt(1 - zeta^2))
% zeta = -[ln(MU)]/sqrt[pi^2 + ln^2(MU)] I calculate and go to zgrid
% Equations for the region of a cone
xv = 1;         % cone vertex on the horizontal axis
xe = 0.7;       % chosen on zgrid after finding zeta = 0.7
ye = 0.225;     % evaluation of y at xe = 0.7
gama = atan(ye/(1 -xe)); % opening angle of the cone/2 in radians
R11c = [-xv*sin(gama)*2   0  ; 0  -xv*sin(gama)*2];
R12c = [sin(gama) cos(gama)  ; -cos(gama)  sin(gama)];
R22c = 0.001*eye(2);

% Equations for the region of an ellipse
zeta = 0.7;
phi = acos(zeta);  % damping factor angle in radians
x0 = -exp(-pi/(tan(phi)));
xse = (1 + x0)/2;
ak = (1 - x0)/2;
bk = ye*ak/sqrt(ak^2 -(xe -xse)^2);
R11e = [-1  -xse/ak; -xse/ak -1];
R12e = [0 (1/ak -1/bk)/2; (1/ak + 1/bk)/2 0];
R22e = 0.001*eye(2);

% Assembly of the intersection of two regions: cone and ellipse
Z = zeros(2);
R11ce = [R11c Z; Z' R11e];
R12ce = [R12c Z; Z' R12e];
R22ce = [R22c Z; Z' R22e];

% % Joining cone-ellipse with the circle centered at the origin
% Z = zeros(4,1);
% R11 = [R11ce Z; Z' R11cc];
% R12 = [R12ce Z; Z' R12cc];
% R22 = [R22ce Z; Z' R22cc];

M11 = kron(R11, P) + kron(R12, (A*H + B1*S)) + kron(R12', (A*H + B1*S)');
M22 = kron(R22, (P - H - H'));
M12 = kron(R12', (P - H')) + kron(R22, (A*H + B1*S));
LMIs = [LMIs; [M11,M12; M12',M22]<=0];

M11 = kron(R11, P) + kron(R12, (A*H + B2*S)) + kron(R12', (A*H + B2*S)');
M22 = kron(R22, (P - H - H'));
M12 = kron(R12', (P - H')) + kron(R22, (A*H + B2*S));
LMIs = [LMIs; [M11,M12; M12',M22]<=0];

%% Inserts the controller gain limitations
LMIs = [LMIs; H >= eye(nx,nx)];
LMIs = [LMIs; [H, S'; S, eta*(v0^2)] >= 0];
%% Solves the LMIs
solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
% checkset(LMIs) % Uncomment to see the checkset
p = min(checkset(LMIs));
output.feas = 0;
output.p = p;
% Verifies if a P has been found
if p > 0
    
    output.H = double(H);
    output.S = double(S);
    output.feas = 1;
    %% Plots the regions and the poles
    % plot the ts limit
    xCentro = 0;
    yCentro = 0;
    theta = 0 : 0.01 : 2*pi;
    r = exp(-4/ts);
    x = r * cos(theta) + xCentro;
    y = r * sin(theta) + yCentro;
    p1 = plot(x, y, 'r--', 'linewidth', 2);
    axis square;
    xlim([-1 1]);
    ylim([-1 1]);
    grid on;
    axis equal;
    hold on
    % plot the ellipse until it touches the cone
    xc = xse;
    yc = 0;
    xr = ak;
    yr = bk;
    theta = atan(ye/(xe-xc))+0.33 : 0.001 : 2*pi - atan(ye/(xe-xc))-0.33 ;
    x = xr * cos(theta) + xc;
    y = yr * sin(theta) + yc;
    hold on;
    p2 = plot(x, y, 'LineWidth', 2, 'Color', [0,0,1]);
    axis square;
    xlim([-1 1]);
    ylim([-1 1]);
    grid on;
    hold on
    
    % plot the cone until it touches the ellipse
    line([xv, xe], [0, ye], 'LineWidth', 2, 'Color', [0,0,1]);
    line([xv, xe], [0, -ye], 'LineWidth', 2, 'Color', [0,0,1]);
    
    % plot the poles of CL
    p3 = plot(real(eig(A + minf*B*output.S*pinv(output.H))), imag(eig(A + minf*B*output.S*pinv(output.H))), 'mx', 'MarkerSize',20,'linewidth',3);
    p4 = plot(real(eig(A + msup*B*output.S*pinv(output.H))), imag(eig(A + msup*B*output.S*pinv(output.H))), 'gs', 'MarkerSize',20 ,'linewidth',3);
    hold on
    
    % plot the cardioid for the defined zeta
    sigma = 0:0.01:pi-0.2;
    zx = exp(-sigma).*cos(sigma.*tan(pi-phi));
    zy = exp(-sigma).*sin(sigma.*tan(pi-phi));
    plot(zx,zy, 'k:', 'linewidth', 2);
    a = get (gca,'XTickLabel');
    set(gca,'XTickLabel', a,'fontsize',22,'FontName','Times')
    zx = exp(-sigma).*cos(sigma.*tan(-pi+phi));
    zy = exp(-sigma).*sin(sigma.*tan(-pi+phi));
    p5 = plot(zx,zy, 'k:', 'linewidth', 2);
    legend([p1 p2 p3 p4 p5], 'limit $T_s$','approx. limit MUP','poles $G_1$ CL', 'poles $G_2$ CL','MUP limit','FontName','Times')
    legend({},'FontSize',22,'FontName','Times', 'Interpreter', 'latex')
    set(gcf,'Color', 'white')
    xticklabels('auto')
    set(gcf,'Color', 'white', 'Position',  [100, 100, 700, 350])
    title('CL allocation region')
    hold off
    
end
end
