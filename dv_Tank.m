function hddot = dv_Tank(u, hddot, Ts)
% Luís Henrique dos Santos 16/02/2025
% Implementation of the numerical example from the article A new approach 
% for robust control based on parametric Hammerstein models

% Outputs
h11 = hddot(1);
h22 = hddot(2);

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
S = 0.014;          % Cross-sectional area of tank 1 and tank 2
a = 5e-5;           % Outlet area of tank 1 and tank 2
eta = 2.4e-3;       % Pump flow rate constant
g = 9.81;           % Gravitational constant 

%--------------------------------------------------------------------------
%  saturation
%--------------------------------------------------------------------------

%saturation of h11
if h11 > 0.5
    sath11 = 0.5;
elseif h11 < 0    
    sath11 = 0;
else
    sath11 = h11;
end
 
%saturation of h22

if h22 > 0.5
    sath22 = 0.5;
elseif h22 < 0
    sath22 = 0;
else
    sath22 = h22;
end

%saturation of u

if u > 5
    satu = 5;
elseif u < -5
    satu = -5;
else
    satu = u;
end
  

%--------------------------------------------------------------------------
% diff equations (35) in Rayouf 2018
%--------------------------------------------------------------------------
hddot(1,:) = h11 + Ts*(eta*satu - (a/S)*(sqrt(2*g*sath11)));
hddot(2,:) = h22 + Ts*((a/S)*(sqrt(2*g*sath11)) - (a/S)*(sqrt(2*g*sath22)));


