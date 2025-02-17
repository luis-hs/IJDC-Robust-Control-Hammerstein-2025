% Luís Henrique dos Santos 16/02/2025
% Implementation of the numerical example from the article A new approach
% for robust control based on parametric Hammerstein models

clearvars
close all
clc

load('IdentTank.mat') % Loads the identified model of the two-tank cascade system
% into the variables Ao, Bo, Co, and alpha_hat

minf = 0.85;   % Maximum slope of the lower bound line for the intermediate signal
msup = 1.15;   % Minimum slope of the upper bound line for the intermediate signal
timef = 1001;  % Number of simulation iterations
ycr = 0.2;     % Desired output level for the second tank
eta = 1;       % Parameter related to the control law design
v0 = 750.1;    % Parameter related to the control law design
ts = 88;       % Desired settling time
%--------------------------------------------------------------------------
% reference signal
%--------------------------------------------------------------------------
yr = linspace(ycr,ycr,timef);
% the following three lines can be commented out to simulate just one step
yr(401:600) = 0.1*ones(size(yr(401:600)));
yr(601:800) = 0.15*ones(size(yr(601:800)));
yr(801:end) = 0.05*ones(size(yr(801:end)));


Co2 = Co(2,:); % Selects only the C matrix row corresponding to state 2 (the controlled state)
n = 1;         % Number of system outputs


%--------------------------------------------------------------------------
%           Expanding the system with the new state variable
%--------------------------------------------------------------------------
A = [Ao, zeros(size(Ao,1), size(Co2,1)); -Co2, eye(size(Co2,1), size(Co2,1))];
B = [Bo; zeros(n,size(Bo,2))];             % Input matrix
Br = [zeros(size(Bo)); eye(n,size(Bo,2))]; % Reference matrix
C = [Co2, zeros(size(Co2,1), n)];          % Output matrix

% Beta for the case where [u u^2], according to Rayouf, Ghorbel, and Braiek (2018)
beta = [1/alpha_hat(1) , -alpha_hat(2)/(alpha_hat(1)^3)];

%--------------------------------------------------------------------------
%        Compute the PI controller that stabilizes the model
%--------------------------------------------------------------------------

% Calls the function to design the controller that places the closed-loop poles
% in the region specified by the settling time ts, while also satisfying
% the other parameters
controller = dsolvelmi_region(A, B, minf, msup, ts, eta, v0);
%%
x0_zic = [0;0;0];          % Zero initial condition
x0_nic = [0.1;0.05;0.1];   % Nonzero initial condition
if controller.feas ~= 0
    % K = controller.S*pinv(controller.H); 
    K = -1*[6.9726  25.0562 -0.7936]; % use this to get the result of the article
    
    %----------------------------------------------------------------------
    %       Initialize all variables used in the system simulation
    %----------------------------------------------------------------------
    N_std = 1000; % Number of iterations to generate data for the standard deviation shadow
    % zic indicates null (zero) initial conditions
    % nic indicates non null initial condition
    % al indicates algebraic linearization
    % or indicates linearization by operation regions
    % std indicates that the variable is for plot shadow standart deviation
    % m indicates model
    % p indicates process
    % y indicates output
    % x indicates state variable
    % u indicates input of Hammerstein model
    % v indicates the intermediate signal
    % v_hat the control signal
    ym_zic_al_std = zeros(N_std,timef-1);
    yp_zic_al_std = zeros(N_std,timef-1);
    ym_zic_or_std = zeros(N_std,timef-1);
    yp_zic_or_std = zeros(N_std,timef-1);
    ym_nic_al_std = zeros(N_std,timef-1);
    yp_nic_al_std = zeros(N_std,timef-1);
    ym_nic_or_std = zeros(N_std,timef-1);
    yp_nic_or_std = zeros(N_std,timef-1);
    for it_std = 1:N_std
        % inicialization of variables of model with (zic) and (al)
        xm_zic_al = zeros(size(A, 2),timef);
        xm_zic_al(:,1) = x0_zic;
        ym_zic_al = zeros(1, timef);
        vm_hat_zic_al = zeros(1, timef);
        vm_zic_al = zeros(1, timef);
        um_zic_al = zeros(1, timef);
        
        % inicialization of variables of model with (zic) and (or)
        xm_zic_or = zeros(size(A, 2),timef);
        xm_zic_or(:,1) = x0_zic;
        ym_zic_or = zeros(1, timef);
        vm_hat_zic_or = zeros(1, timef);
        vm_zic_or = zeros(1, timef);
        um_zic_or = zeros(1, timef);
        
        % inicialization of variables of model with (nic) and (al)
        xm_nic_al = zeros(size(A, 2),timef);
        xm_nic_al(:,1) = x0_nic;
        ym_nic_al = zeros(1, timef);
        vm_hat_nic_al = zeros(1, timef);
        vm_nic_al = zeros(1, timef);
        um_nic_al = zeros(1, timef);
        
        % inicialization of variables of model with (nic) and (or)
        xm_nic_or = zeros(size(A, 2),timef);
        xm_nic_or(:,1) = x0_nic;
        ym_nic_or = zeros(1, timef);
        vm_hat_nic_or = zeros(1, timef);
        vm_nic_or = zeros(1, timef);
        um_nic_or = zeros(1, timef);
        
        % inicialization of variables of process with (zic) and (al)
        xp_zic_al = zeros(size(A, 2), timef);
        xp_zic_al(:,1) = x0_zic;
        yp_zic_al = zeros(1, timef);
        vp_hat_zic_al = zeros(1, timef);
        up_zic_al = zeros(1, timef);
        
        % inicialization of variables of process with (nic) and (al)
        xp_nic_al = zeros(size(A, 2), timef);
        xp_nic_al(:,1) = x0_nic;
        yp_nic_al = zeros(1, timef);
        vp_hat_nic_al = zeros(1, timef);
        up_nic_al = zeros(1, timef);
        
        % inicialization of variables of process with (zic) and (or)
        xp_zic_or = zeros(size(A, 2), timef);
        xp_zic_or(:,1) = x0_zic;
        yp_zic_or = zeros(1, timef);
        vp_hat_zic_or = zeros(1, timef);
        up_zic_or = zeros(1, timef);
        
        % inicialization of variables of process with (nic) and (or)
        xp_nic_or = zeros(size(A, 2), timef);
        xp_nic_or(:,1) = x0_nic;
        yp_nic_or = zeros(1, timef);
        vp_hat_nic_or = zeros(1, timef);
        up_nic_or = zeros(1, timef);
        
        
        mu1 = 0.3;   % manually selects the values of the polytope parameters
        mu2 = 0.7;
        sigma = minf*mu1 + msup*mu2;
        yp_al_max_val = ycr;
        ym_al_max_val = ycr;
        yp_or_max_val = ycr;
        ym_or_max_val = ycr;
        
        %------------------------------------------------------------------
        % noise: set the parameter 0.01 to 0 if you want to remove the noise
        %------------------------------------------------------------------
        w = randn(timef,1)*0.01;
        erro_or = zeros(1,timef);
        erro_al = zeros(1,timef);
        itse_or = zeros(1,timef);
        itse_al = zeros(1,timef);
        pert = zeros(1,timef);   % perturbation
        pert(250:end) = ones(size(pert(250:end)));
        %------------------------------------------------------------------
        %  Closes the loop using the same controller K for both systems in
        %  two methodologies: algebraic and operating region-based
        %------------------------------------------------------------------
        % The simulation was performed step by step, meaning that everything is
        % computed at each iteration, and the next iteration starts with updated values.
        
        %% For zero initial conditions
        for k = 1:timef
            %--------------------------------------------------------------
            %    Closes the loop in the Hammerstein model using the
            %    algebraic methodology (zero initial conditions)
            %--------------------------------------------------------------
            vm_hat_zic_al(k) = K*xm_zic_al(:,k)+pert(k);   % control law
            vm_hat_zic_al(k) = sigma*vm_hat_zic_al(k);
            % um_zic_al(k) = -0.182796 + 0.0107527*sqrt(465000*vm_hat_zic_al(k) + 289);
            um_zic_al(k) = -0.18 + 0.01*sqrt(465000*vm_hat_zic_al(k) + 289);
            if um_zic_al(k) > 5       % Saturates the control signal
                um_zic_al(k) = 5;
            elseif um_zic_al(k) < -5
                um_zic_al(k) = -5;
            end
            vm_zic_al(k) = alpha_hat(1)*um_zic_al(k) + alpha_hat(2)*(um_zic_al(k)^2); % Non-linearity
            xm_zic_al(:,k+1) = A*xm_zic_al(:,k) + B*vm_zic_al(k) + Br*yr(k); % Closing the loop
            ym_zic_al(k) = C*xm_zic_al(:,k+1) + w(k); % Computes the output
            
            if k < 201 && (ym_zic_al(k) < ycr*0.98 || ym_zic_al(k) > ycr*1.02)
                tsm_al = k+1; % Stores the settling time
            end
            if ym_zic_al(k) > ycr && ym_zic_al(k) > ym_al_max_val
                ym_al_max_val = ym_zic_al(k); % Stores the peak amplitude
            end
            
            %--------------------------------------------------------------
            %    Closes the loop in the simulated process using the
            %    algebraic methodology (zero initial conditions)
            %--------------------------------------------------------------
            vp_hat_zic_al(k) = K*xp_zic_al(:,k)+pert(k); % control law
            vp_hat_zic_al(k) = sigma*vp_hat_zic_al(k);
            % up_zic_al(k) = -0.182796 + 0.0107527*sqrt(465000*vp_hat_zic_al(k) + 289);
            up_zic_al(k) = -0.18 + 0.01*sqrt(465000*vp_hat_zic_al(k) + 289); % Pre-inverts the non-linearity
            % Applies the inverse result to the process to obtain level 2
            xp_zic_al(1:end-1,k+1) = dv_Tank(up_zic_al(k), xp_zic_al(1:end-1, k), Ts);
            yp_zic_al(k) = xp_zic_al(end-1,k) + w(k); % Computes the output
            % Computes the third state variable of the process (qk from the PI controller)
            xp_zic_al(end,k+1) = xp_zic_al(end,k) -Co(2,:)*xp_zic_al(1:end-1,k) + yr(k);
            
            if k < 201 && (yp_zic_al(k) < ycr*0.98 || yp_zic_al(k) > ycr*1.02)
                tsp_al = k+1; % Stores the settling time
            end
            if yp_zic_al(k) > ycr && yp_zic_al(k) > yp_al_max_val
                yp_al_max_val = yp_zic_al(k); % Stores the peak amplitude
            end
            
            %--------------------------------------------------------------
            %   Closes the loop in the Hammerstein model using the
            %   operating region-based methodology (zero initial conditions)
            %--------------------------------------------------------------
            vm_hat_zic_or(k) = K*xm_zic_or(:,k) + pert(k);
            vm_hat_zic_or(k) = sigma*vm_hat_zic_or(k);
            um_zic_or(k) =  beta(1)*vm_hat_zic_or(k) + beta(2)*vm_hat_zic_or(k)^2;
            if um_zic_or(k) > 5
                um_zic_or(k) = 5;
            elseif um_zic_or(k) < -5
                um_zic_or(k) = -5;
            end
            vm_zic_or(k) = alpha_hat(1)*um_zic_or(k) + alpha_hat(2)*(um_zic_or(k)^2);
            xm_zic_or(:,k+1) = A*xm_zic_or(:,k) + B*vm_zic_or(k) + Br*yr(k);
            ym_zic_or(k) = C*xm_zic_or(:,k+1) + w(k);
            
            if k < 201 && (ym_zic_or(k) < ycr*0.98 || ym_zic_or(k) > ycr*1.02)
                tsm_or = k+1;
            end
            if ym_zic_or(k) > ycr && ym_zic_or(k) > ym_or_max_val
                ym_or_max_val = ym_zic_or(k);
            end
            
            %--------------------------------------------------------------
            %   Closes the loop in the process using the
            %   operating region-based methodology (zero initial conditions)
            %--------------------------------------------------------------
            vp_hat_zic_or(k) = K*xp_zic_or(:,k)+pert(k);
            vp_hat_zic_or(k) = sigma*vp_hat_zic_or(k);
            up_zic_or(k) = beta(1)*vp_hat_zic_or(k) + beta(2)*vp_hat_zic_or(k)^2;
            up_zic_or(k) = abs(up_zic_or(k));
            xp_zic_or(1:end-1,k+1) = dv_Tank(up_zic_or(k), xp_zic_or(1:end-1, k), Ts);
            yp_zic_or(k) = xp_zic_or(end-1,k) + w(k);
            xp_zic_or(end,k+1) = xp_zic_or(end,k) -Co(2,:)*xp_zic_or(1:end-1,k) + yr(k);
            
            if k < 201 && (yp_zic_or(k) < ycr*0.98 || yp_zic_or(k) > ycr*1.02)
                tsp_or = k+1;
            end
            if yp_zic_or(k) > ycr && yp_zic_or(k) > yp_al_max_val
                yp_or_max_val = yp_zic_or(k);
            end
            
            erro_or(k+1) = yp_zic_or(k) - ycr;
            erro_al(k+1) = yp_zic_al(k) - ycr;
            itse_or(k+1) = k*(erro_or(k+1).^2 + erro_or(k).^2)/2;
            itse_al(k+1) = k*(erro_al(k+1).^2 + erro_al(k).^2)/2;
        end
        ym_zic_al_std(it_std,:) = ym_zic_al(1:end-1);
        yp_zic_al_std(it_std,:) = yp_zic_al(1:end-1);
        ym_zic_or_std(it_std,:) = ym_zic_or(1:end-1);
        yp_zic_or_std(it_std,:) = yp_zic_or(1:end-1);
        
        %% For non-zero initial conditions
        for k = 1:timef
            %--------------------------------------------------------------
            %    Closes the loop in the Hammerstein model using the
            %    algebraic methodology (non-zero initial conditions)
            %--------------------------------------------------------------
            vm_hat_nic_al(k) = K*xm_nic_al(:,k)+pert(k);
            vm_hat_nic_al(k) = sigma*vm_hat_nic_al(k);
            um_nic_al(k) = -0.18 + 0.01*sqrt(465000*vm_hat_nic_al(k) + 289);
            if um_nic_al(k) > 5
                um_nic_al(k) = 5;
            elseif um_nic_al(k) < -5
                um_nic_al(k) = -5;
            end
            
            vm_nic_al(k) = alpha_hat(1)*um_nic_al(k) + alpha_hat(2)*(um_nic_al(k)^2);
            xm_nic_al(:,k+1) = A*xm_nic_al(:,k) + B*vm_nic_al(k) + Br*yr(k);
            ym_nic_al(k) = C*xm_nic_al(:,k+1) + w(k);
            
            %--------------------------------------------------------------
            %    Closes the loop for the simulated process using the
            %    algebraic methodology (non-zero initial conditions)
            %--------------------------------------------------------------
            vp_hat_nic_al(k) = K*xp_nic_al(:,k)+pert(k);
            vp_hat_nic_al(k) = sigma*vp_hat_nic_al(k);
            % up_nic_al(k) = -0.182796 + 0.0107527*sqrt(465000*vp_hat_nic_al(k) + 289);
            up_nic_al(k) = -0.18 + 0.01*sqrt(465000*vp_hat_nic_al(k) + 289);
            xp_nic_al(1:end-1,k+1) = dv_Tank(up_nic_al(k), xp_nic_al(1:end-1, k), Ts);
            yp_nic_al(k) = xp_nic_al(end-1,k) + w(k);
            xp_nic_al(end,k+1) = xp_nic_al(end,k) -Co(2,:)*xp_nic_al(1:end-1,k) + yr(k);
            
            
            
            %--------------------------------------------------------------
            %   Closes the loop for the Hammerstein model using the
            %   operation regions methodology (non-zero initial conditions)
            %--------------------------------------------------------------
            vm_hat_nic_or(k) = K*xm_nic_or(:,k)+pert(k);
            vm_hat_nic_or(k) = sigma*vm_hat_nic_or(k);
            um_nic_or(k) =  beta(1)*vm_hat_nic_or(k) + beta(2)*vm_hat_nic_or(k)^2;
            if um_nic_or(k) > 5
                um_nic_or(k) = 5;
            elseif um_nic_or(k) < -5
                um_nic_or(k) = -5;
            end
            vm_nic_or(k) = alpha_hat(1)*um_nic_or(k) + alpha_hat(2)*(um_nic_or(k)^2);
            xm_nic_or(:,k+1) = A*xm_nic_or(:,k) + B*vm_nic_or(k) + Br*yr(k);
            ym_nic_or(k) = C*xm_nic_or(:,k+1) + w(k);
            
            
            %--------------------------------------------------------------
            %   Closes the loop for the process using the
            %   operation regions methodology (non-zero initial conditions)
            %--------------------------------------------------------------
            vp_hat_nic_or(k) = K*xp_nic_or(:,k)+pert(k);
            vp_hat_nic_or(k) = sigma*vp_hat_nic_or(k);
            up_nic_or(k) = beta(1)*vp_hat_nic_or(k) + beta(2)*vp_hat_nic_or(k)^2;
            up_nic_or(k) = abs(up_nic_or(k));
            xp_nic_or(1:end-1,k+1) = dv_Tank(up_nic_or(k), xp_nic_or(1:end-1, k), Ts);
            yp_nic_or(k) = xp_nic_or(end-1,k) + w(k);
            xp_nic_or(end,k+1) = xp_nic_or(end,k) -Co(2,:)*xp_nic_or(1:end-1,k) + yr(k);
            
        end
        ym_nic_al_std(it_std,:) = ym_nic_al(1:end-1);
        yp_nic_al_std(it_std,:) = yp_nic_al(1:end-1);
        ym_nic_or_std(it_std,:) = ym_nic_or(1:end-1);
        yp_nic_or_std(it_std,:) = yp_nic_or(1:end-1);
    end
    fprintf('ITSE with algebraic: %.4f\n', sum(itse_al(1:200)));% integral criterion with algebraic inverse for step of 0.2 only
    fprintf('ITSE with operation regions: %.4f\n', sum(itse_or(1:200))); % integral criterion with inverse by operation regions for step of 0.2 only
    disp('For the measurement below to be correct, set the noise w to zero.')
    fprintf('Ts with algebraic inversion: %i\n', tsp_al);% settling time
    fprintf('Ts with operation regions inversion: %i\n', tsp_or); % settling time
    
    
    %----------------------------------------------------------------------
    % Plots the results of the system and model simulation in closed loop
    %----------------------------------------------------------------------
    % test signal for inversion effectiveness
    vt_hat = 0:0.05:5;
    ut_al_ex = -0.182796 + 0.0107527*sqrt(465000*vt_hat + 289);
    ut_al_ap = -0.18 + 0.01*sqrt(465000*vt_hat + 289);
    ut_or = beta(1)*vt_hat + beta(2)*(vt_hat.^2);
    
    vt_ap = alpha_hat(1)*ut_al_ap + alpha_hat(2)*(ut_al_ap.^2);
    vt_ex = alpha_hat(1)*ut_al_ex + alpha_hat(2)*(ut_al_ex.^2);
    vt_or = alpha_hat(1)*ut_or + alpha_hat(2)*(ut_or.^2);
    
    figure
    warning off;
    plot(vt_hat, vt_hat,'k-')
    hold on
    plot(vt_hat, vt_ex,'r--')
    plot(vt_hat, vt_ap,'b-.')
    title('Inverse validation', 'FontSize',16)
    leg = legend({'linear input', 'exact inverse', 'inverse with parametric variation'});
    legend({},'FontSize',16,'Color','white')
    grid
    xlabel('$$\tilde{\hat{v}}_k$$','FontSize',16,'Interpreter','Latex');
    ylabel('amplitude','FontSize',16)
    set(gca,'FontName','Times','FontSize',14)
    set(gcf,'Color','w','Position',[50 200 500 350])
    
    
    figure
    stairs(up_zic_al(1:timef -1), 'k-', 'linewidth', 1.5);
    hold on
    stairs(um_zic_al(1:timef-1), 'r--', 'linewidth',1.5);
    stairs(pert(1:timef-1), 'b:', 'linewidth',1.5);
    legend('$\hat{v}_k$ real', '$\hat{v}_k$ model', 'perturbation')
    legend({},'Interpreter','latex','FontSize',16)
    title('Control signal using algebraic linearization')
    grid
    xlabel('kT [s]','FontSize',16);
    ylabel('amplitude','FontSize',16);
    set(gca,'FontName','Times','FontSize',16)
    set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);
    
    
    figure
    plot(1:timef-1, yr(1:timef-1), '-k', 'linewidth', 1.5);
    hold on
    config = stdshade(yp_nic_al_std,0,[128 128 128]/255,1:timef-1,'-k',1);
    config.Color = [0 1 0];
    hold on
    config = stdshade(ym_nic_al_std,0,[128 128 128]/255,1:timef-1,'-k',1);
    config.Color = [1 0 1];
    hold on
    stdshade(yp_zic_al_std,0.3,'r',1:timef-1,'--r',1.5);
    hold on
    stdshade(ym_zic_al_std,0.3,'b',1:timef-1,':b',1.5);
    ylim([0 0.251])
    hold off
    legend('$y_r$','','$y_k$ real $x_0 \neq 0$','','$\hat{y}_k$ model $x_0 \neq 0$','', '$y_k$ real','','$\hat{y}_k$ model' , 'Location','SouthEast')
    legend({},'FontSize',16,'Interpreter','latex')
    title('System response to a multi-reference with AL inverse')
    grid
    xlabel('kT - iterations');
    ylabel('System output');
    set(gca,'FontName','Times','FontSize',16)
    set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);
    
    
    figure
    plot(1:timef-1, yr(1:timef-1), '-k', 'linewidth', 1.5);
    hold on
    config = stdshade(yp_nic_or_std,0,[128 128 128]/255,1:timef-1,'-k',1);
    config.Color = [0 1 0];
    hold on
    config = stdshade(ym_nic_or_std,0,[128 128 128]/255,1:timef-1,'-k',1);
    config.Color = [1 0 1];
    hold on
    stdshade(yp_zic_or_std,0.3,'r',1:timef-1,'--r',1.5);
    hold on
    stdshade(ym_zic_or_std,0.3,'b',1:timef-1,':b',1.5);
    ylim([0 0.501])
    xlim([0 timef-1])
    hold off
    legend('$y_r$','','$y_k$ real $x_0 \neq 0$','','$\hat{y}_k$ model $x_0 \neq 0$','', '$y_k$ real','','$\hat{y}_k$ model' , 'Location','SouthEast')
    legend({},'FontSize',16,'Interpreter','latex')
    grid
    title('System response to a multi-reference with OR inverse')
    xlabel('kT - iterations');
    ylabel('System output');
    set(gca,'FontName','Times','FontSize',16)
    set(gcf,'Color', 'white','renderer','painters','Position',[50 100 800 300]);
    xticklabels('auto')
    warning on;
else
    tsp_al = 1000;
end
