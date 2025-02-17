% Luís Henrique dos Santos 16/02/2025
% Implementation of the numerical example from the article A new approach
% for robust control based on parametric Hammerstein models

clearvars -except up_zic_al
load('IdentTank.mat') % Loads the identified model of the two-tank cascade system
% into the variables Ao, Bo, Co, and alpha_hat

close all
clc

try
    usmfr = up_zic_al;
catch ME
    disp('First execute proposed_control to provide up_zic_al')
end

minf = 0.85;   
msup = 1.15;   
timef = 1001;
ycr = 0.2;

yr = linspace(ycr,ycr,timef);
yr(401:600) = 0.1*ones(size(yr(401:600)));
yr(601:800) = 0.15*ones(size(yr(601:800)));
yr(801:end) = 0.05*ones(size(yr(801:end)));
Co1 = Co(2,:); 
n = 1;         

%--------------------------------------------------------------------------
%           Aumentando o sistema com a nova variável de estado
%--------------------------------------------------------------------------
A = Ao;
B = Bo;
C = Co1;


beta = [1/alpha_hat(1) , -alpha_hat(2)/(alpha_hat(1)^3)];

%--------------------------------------------------------------------------
%        Calcula o controla proporcional K que estabilize o sistema
%--------------------------------------------------------------------------
controlador = dsolvelmi_rayouf(A, B, minf, msup);
if controlador.feas ~= 0
    P = pinv(controlador.X);
    K = controlador.R*P;
end
x0nn = [0.1;0.05];
x0 = [0;0];
%--------------------------------------------------------------------------
%       Inicializa todas variáveis usadas na simulação do sistema
%--------------------------------------------------------------------------
G1 = A - minf*B*K;
G2 = A - msup*B*K;

timefnn = timef;
ns = 100;
ymn = zeros(ns,timef-1);
ysn = zeros(ns,timef-1);
yman = zeros(ns,timef-1);
ysan = zeros(ns,timef-1);
ymnnn = zeros(ns,timefnn-1);
ysnnn = zeros(ns,timefnn-1);
ymannn = zeros(ns,timefnn-1);
ysannn = zeros(ns,timefnn-1);
for mc = 1:ns

N = zeros(1, timef);
Na = zeros(1, timef);
u = zeros(1, timef);

x1 = zeros(size(A, 2),timef);
yMmin = zeros(1, timef);
d1 = zeros(1, timef);
d1a = zeros(1, timef);

x2 = zeros(size(A, 2),timef);
yMmax = zeros(1, timef);
d2 = zeros(1, timef);
d2a = zeros(1, timef);

mu1 = zeros(1, timef);
mu2 = mu1;
sigma = zeros(1, timef);
mu1a = zeros(1, timef);
mu2a = mu1a;
sigmaa = zeros(1, timef);

v = zeros(1, timef);
d = zeros(1, timef);

xm = zeros(size(A, 2),timef);
xm(:,1) = x0;
ym = zeros(1, timef);
vm_hat = zeros(1, timef);
vm = zeros(1, timef);
um = zeros(1, timef);
xma = zeros(size(A, 2),timef);
xma(:,1) = x0;
yma = zeros(1, timef);
vm_hata = zeros(1, timef);
vma = zeros(1, timef);
uma = zeros(1, timef);


xs = zeros(size(A, 2), timef);
xs(:,1) = x0;
us = zeros(1, timef);
ysys = zeros(1, timef);
vs_hat = zeros(1, timef);


xsa = zeros(size(A, 2), timef);
xsa(:,1) = x0;
usa = zeros(1, timef);
satua = zeros(1, timef);
ysysa = zeros(1, timef);
vs_hata = zeros(1, timef);
vsa = zeros(1, timef);


% Simula condições iniciais não nulas
xmnn = zeros(size(A, 2),timef);
xmnn(:,1) = x0nn;
ymnn = zeros(1, timefnn);
vm_hatnn = zeros(1, timefnn);
vmnn = zeros(1, timefnn);
umnn = zeros(1, timefnn);
xmann = zeros(size(A, 2),timefnn);
xmann(:,1) = x0nn;
ymann = zeros(1, timefnn);
vm_hatann = zeros(1, timefnn);
vmann = zeros(1, timefnn);
umann = zeros(1, timefnn);
xsnn = zeros(size(A, 2), timefnn);
xsnn(:,1) = x0nn;
usnn = zeros(1, timefnn);
ysysnn = zeros(1, timefnn);
vs_hatnn = zeros(1, timefnn);

xsann = zeros(size(A, 2), timefnn);
xsann(:,1) = x0nn;
usann = zeros(1, timefnn);
ysysann = zeros(1, timefnn);
vs_hatann = zeros(1, timefnn);



Ninf = ( 1/minf )/( C* ( (eye(size(G1)) - G1) \ B ) );
Nsup = ( 1/msup )/( C* ( (eye(size(G2)) - G2) \ B ) );


ysysm = ycr;
ymm = ycr;
ysysma = ycr;
ymma = ycr;

w = randn(timef,1)*0.01;
erroa = zeros(1,timef);
erro = zeros(1,timef);
itsea = zeros(1,timef);
itse = zeros(1,timef);
pert = zeros(1,timef);
pert(250:end) = ones(size(pert(250:end)));
%--------------------------------------------------------------------------
%  Fecha a malha para o mesmo controlador K nos dois sistemas
%--------------------------------------------------------------------------
% A simulação foi feita passo a passo, ou seja, calcula tudo a cada
% iteração e parte para a próxima iteração com valores atualizados
for k = 1:timef
    
    x1(:,k+1) = G1*x1(:,k) + minf*B*Ninf*ycr;  
    yMmin(k) = C*x1(:,k);
    
    x2(:,k+1) = G2*x2(:,k) + msup*B*Nsup*ycr;  
    yMmax(k) = C*x2(:,k);
    
    d1(k) = abs(ysys(k) - yMmin(k));          
    d2(k) = abs(ysys(k) - yMmax(k));          
    
    if (k >= 3)
        mu1(k) = d2(k)/(d1(k)+d2(k));         
        mu2(k) = 1 - mu1(k);                  
    else
        mu1(k) = 0.5;                         
        mu2(k) = 1 - mu1(k);
    end
    sigma(k) = minf*mu1(k) + msup*mu2(k);
    

    N(k) = ( 1/sigma(k) )/( C* ( (eye(size(G1)) - (mu1(k)*G1 + mu2(k)*G2)) \ B ) );
    
    %--------------------------------------------------------------------------
    %          Fecha a malha manualmente no modelo de Hammerstein
    %--------------------------------------------------------------------------
    vm_hat(k) = -K*xm(:,k) + N(k)*yr(k) +pert(k);        
    vm_hat(k) = sigma(k)*vm_hat(k);   
    %um(k) = -0.182796 + 0.0107527*sqrt(465000*vm_hat(k) + 289);
    um(k) = -0.18 + 0.01*sqrt(465000*vm_hat(k) + 289);
    if um(k) > 5                   
        um(k) = 5;                 
    elseif um(k) < -5               
        um(k) = -5;                 
    end
    vm(k) = alpha_hat(1)*um(k) + alpha_hat(2)*(um(k)^2);
    xm(:,k+1) = A*xm(:,k) + B*vm(k) ; 
    ym(k) = C*xm(:,k) + w(k);
    
    
    if ym(k) < ycr*0.98 || ym(k) > ycr*1.02 
        tsm = k+1;
    end
    if ym(k) > ycr && ym(k) > ymm
        ymm = ym(k);
    end
    
    %--------------------------------------------------------------------------
    %             Fecha a malha manualmente para o sistema real
    %--------------------------------------------------------------------------
    vs_hat(k) = -K*xs(:,k)*H  + N(k)*yr(k)+pert(k);    
    vs_hat(k) = sigma(k)*vs_hat(k);  
    %us(k) = -0.182796 + 0.0107527*sqrt(465000*vs_hat(k) + 289);
    us(k) = -0.18 + 0.01*sqrt(465000*vs_hat(k) + 289);
    us(k) = sign(real(us(k)))*abs(us(k));
    xs(:,k+1) = dv_Tank(us(k), xs(:, k), Ts);
    ysys(k) = xs(2,k) + w(k);
    
    if ysys(k) < ycr*0.98 || ysys(k) > ycr*1.02 
        tssys = k+1;
    end
    if ysys(k) > ycr && ysys(k)>ysysm
        ysysm = ysys(k);
    end
    
    
    
    d1a(k) = abs(ysysa(k) - yMmin(k));          
    d2a(k) = abs(ysysa(k) - yMmax(k));          
    
    if (k >= 3)
        mu1a(k) = d2a(k)/(d1a(k)+d2a(k));       
        mu2a(k) = 1 - mu1a(k);                  
    else
        mu1a(k) = 0.5;                         
        mu2a(k) = 1 - mu1a(k);
    end
    sigmaa(k) = minf*mu1a(k) + msup*mu2(k);
  
    Na(k) = ( 1/sigmaa(k) )/( C* ( (eye(size(G1)) - (mu1a(k)*G1 + mu2a(k)*G2)) \ B ) );
    
    %--------------------------------------------------------------------------
    %          Fecha a malha manualmente no modelo de Hammerstein
    %--------------------------------------------------------------------------
    vm_hata(k) = -K*xma(:,k) + Na(k)*yr(k)+pert(k);       
    vm_hata(k) = sigmaa(k)*vm_hata(k);   
    uma(k) =  beta(1)*vm_hata(k) + beta(2)*vm_hata(k)^2;
    uma(k) = sign(real(uma(k)))*abs(uma(k));
    if uma(k) > 5                  
        uma(k) = 5;                 
    elseif uma(k) < -5               
        uma(k) = -5;               
    end
    vma(k) = alpha_hat(1)*uma(k) + alpha_hat(2)*(uma(k)^2);
    xma(:,k+1) = A*xma(:,k) + B*vma(k) ; 
    yma(k) = C*xma(:,k) + w(k);
    
    if yma(k) < ycr*0.98 || yma(k) > ycr*1.02 
        tsma = k+1;
    end
    if yma(k) > ycr && yma(k) > ymma
        ymma = yma(k);
    end
    
    %--------------------------------------------------------------------------
    %             Fecha a malha manualmente para o sistema real
    %--------------------------------------------------------------------------
    vs_hata(k) = -K*xsa(:,k)*H  + Na(k)*yr(k)+pert(k);    
    vs_hata(k) = sigmaa(k)*vs_hata(k);  
    usa(k) = beta(1)*vs_hata(k) + beta(2)*vs_hata(k)^2;
    usa(k) = abs(usa(k));
    xsa(:,k+1) = dv_Tank(usa(k), xsa(:, k), Ts);
    ysysa(k) = xsa(2,k+1) + w(k);
    
    if ysysa(k) < ycr*0.98 || ysysa(k) > ycr*1.02 
        tssysa = k+1;
    end
    if ysysa(k) > ycr && ysysa(k)>ysysm
        ysysma = ysysa(k);
    end
        
    erroa(k+1) = ysysa(k) - ycr;
    erro(k+1) = ysys(k) - ycr;
    itsea(k+1) = k*(erroa(k+1).^2 + erroa(k).^2)/2;
    itse(k+1) = k*(erro(k+1).^2 + erro(k).^2)/2;

end
ymn(mc,:) = ym(1:end-1);
ysn(mc,:) = ysys(1:end-1);
yman(mc,:) = yma(1:end-1);
ysan(mc,:) = ysysa(1:end-1);

for k = 1:timefnn
    
    x1(:,k+1) = G1*x1(:,k) + minf*B*Ninf*ycr;  
    yMmin(k) = C*x1(:,k);
    
    x2(:,k+1) = G2*x2(:,k) + msup*B*Nsup*ycr;  
    yMmax(k) = C*x2(:,k);
    
    d1(k) = abs(ysysnn(k) - yMmin(k));          
    d2(k) = abs(ysysnn(k) - yMmax(k));          
    
    if (k >= 3)
        mu1(k) = d2(k)/(d1(k)+d2(k));         
        mu2(k) = 1 - mu1(k);                  
    else
        mu1(k) = 0.5;                         
        mu2(k) = 1 - mu1(k);
    end
    sigma(k) = minf*mu1(k) + msup*mu2(k);
    
   
    N(k) = ( 1/sigma(k) )/( C* ( (eye(size(G1)) - (mu1(k)*G1 + mu2(k)*G2)) \ B ) );
    
    %--------------------------------------------------------------------------
    %          Fecha a malha manualmente no modelo de Hammerstein
    %--------------------------------------------------------------------------
    vm_hatnn(k) = -K*xmnn(:,k) + N(k)*yr(k) +pert(k);        
    vm_hatnn(k) = sigma(k)*vm_hatnn(k);   
    %um(k) = -0.182796 + 0.0107527*sqrt(465000*vm_hat(k) + 289);
    umnn(k) = -0.18 + 0.01*sqrt(465000*vm_hatnn(k) + 289);
    if umnn(k) > 5                   
        umnn(k) = 5;                 
    elseif umnn(k) < -5               
        umnn(k) = -5;                 
    end
    vmnn(k) = alpha_hat(1)*umnn(k) + alpha_hat(2)*(umnn(k)^2);
    xmnn(:,k+1) = A*xmnn(:,k) + B*vmnn(k) ; % Fechando a malha
    ymnn(k) = C*xmnn(:,k) + w(k);
    
    

    
    %--------------------------------------------------------------------------
    %             Fecha a malha manualmente para o sistema real
    %--------------------------------------------------------------------------
    vs_hatnn(k) = -K*xsnn(:,k)*H  + N(k)*yr(k)+pert(k);   
    vs_hatnn(k) = sigma(k)*vs_hatnn(k);  
    %us(k) = -0.182796 + 0.0107527*sqrt(465000*vs_hat(k) + 289);
    usnn(k) = -0.18 + 0.01*sqrt(465000*vs_hatnn(k) + 289);
    usnn(k) = sign(real(usnn(k)))*abs(usnn(k));
    xsnn(:,k+1) = dv_Tank(usnn(k), xsnn(:, k), Ts);
    ysysnn(k) = xsnn(2,k) + w(k);
    
    d1a(k) = abs(ysysann(k) - yMmin(k));        
    d2a(k) = abs(ysysann(k) - yMmax(k));        
    
    if (k >= 3)
        mu1a(k) = d2a(k)/(d1a(k)+d2a(k));       
        mu2a(k) = 1 - mu1a(k);                  
    else
        mu1a(k) = 0.5;                          
        mu2a(k) = 1 - mu1a(k);
    end
    sigmaa(k) = minf*mu1a(k) + msup*mu2(k);
    
    Na(k) = ( 1/sigmaa(k) )/( C* ( (eye(size(G1)) - (mu1a(k)*G1 + mu2a(k)*G2)) \ B ) );
    
    %--------------------------------------------------------------------------
    %          Fecha a malha manualmente no modelo de Hammerstein
    %--------------------------------------------------------------------------
    vm_hatann(k) = -K*xmann(:,k) + Na(k)*yr(k)+pert(k);     
    vm_hatann(k) = sigmaa(k)*vm_hatann(k);   
    umann(k) =  beta(1)*vm_hatann(k) + beta(2)*vm_hatann(k)^2;
    if umann(k) > 5                  
        umann(k) = 5;                
    elseif umann(k) < -5              
        umann(k) = -5;                 
    end
    vmann(k) = alpha_hat(1)*umann(k) + alpha_hat(2)*(umann(k)^2);
    xmann(:,k+1) = A*xmann(:,k) + B*vmann(k) ; 
    ymann(k) = C*xmann(:,k) + w(k);
    

    
    %--------------------------------------------------------------------------
    %             Fecha a malha manualmente para o sistema real
    %--------------------------------------------------------------------------
    vs_hatann(k) = -K*xsann(:,k)*H  + Na(k)*yr(k)+pert(k);    
    vs_hatann(k) = sigmaa(k)*vs_hatann(k);  
    usann(k) = beta(1)*vs_hatann(k) + beta(2)*vs_hatann(k)^2;
    usann(k) = abs(usann(k));
    xsann(:,k+1) = dv_Tank(usann(k), xsann(:, k), Ts);
    ysysann(k) = xsann(2,k+1) + w(k);



end
ymnnn(mc,:) = ymnn(1:end-1);
ysnnn(mc,:) = ysysnn(1:end-1);
ymannn(mc,:) = ymann(1:end-1);
ysannn(mc,:) = ysysann(1:end-1);


end
sum(itse(1:200)) %criterio integral com inversa exata algebrica para degrau de 0,2 apenas
sum(itsea(1:200)) % criterio integral com inversa aproximada

erromodelo = ym(timef -1) - ycr;
erroreal = ym(timef - 1) - ycr;



figure
warning off;
stairs(usmfr(1:timef-1), 'k-', 'linewidth', 1.5)
hold on
stairs(us(1:timef-1), 'r--', 'linewidth',1.5)
stairs(pert(1:timef-1), 'b:', 'linewidth',1.5)
legend('$u_k$ proposed', '$u_k$ [19]','disturbance')
legend({},'Interpreter','latex','FontSize',16)
% title('Sinal controle linearização exata')
grid
xlabel('kT [s]','FontSize',16);
ylabel('Control signal','FontSize',16);
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);


figure
plot(1:timef-1, yr(1:timef-1), '-k', 'linewidth', 1.5);
hold on
config = stdshade(ysnnn,0,[128 128 128]/255,1:timef-1,'-k',1);
config.Color = [0 1 0];
hold on
config = stdshade(ymnnn,0,[128 128 128]/255,1:timef-1,'-k',1);
config.Color = [1 0 1];
hold on
stdshade(ysn,0.3,'r',1:timef-1,'--r',1.5);
hold on
stdshade(ymn,0.3,'b',1:timef-1,':b',1.5);
ylim([0 0.251])
hold off
legend('$y_r$','','$y_k$ real $x_0 \neq 0$','','$\hat{y}_k$ model $x_0 \neq 0$','', '$y_k$ real','','$\hat{y}_k$ model' , 'Location','SouthEast')
legend({},'FontSize',16,'Interpreter','latex')
grid
xlabel('kT - iterations');
ylabel('System output');
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);



figure
plot(1:timef-1, yr(1:timef-1), '-k', 'linewidth', 1.5);
hold on
config = stdshade(ysannn,0,[128 128 128]/255,1:timef-1,'-k',1);
config.Color = [0 1 0];
hold on
config = stdshade(ymannn,0,[128 128 128]/255,1:timef-1,'-k',1);
config.Color = [1 0 1];
hold on
stdshade(ysan,0.3,'r',1:timef-1,'--r',1.5);
hold on
stdshade(yman,0.3,'b',1:timef-1,':b',1.5);
ylim([0 0.52])
xlim([0 timef-1])
hold off
legend('$y_r$','','$y_k$ real $x_0 \neq 0$','','$\hat{y}_k$ model $x_0 \neq 0$','', '$y_k$ real','','$\hat{y}_k$ model' , 'Location','SouthEast')
legend({},'FontSize',16,'Interpreter','latex')
grid
% title('Resposta do sistema para referência de 0,2 inversa aprox')
xlabel('kT - iterations');
ylabel('System output');
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);
xticklabels('auto')
warning on;