%% Assignment 3
%% 
% 

% Reset Everything
close all
clear

% Constant
q_0 = 1.60217653e-19;                   % electron charge
m_0 = 9.10938215e-31;                   % electron mass
meff = 0.26*m_0;                        % electron effective mass
kb = 1.3806504e-23;                     % Boltzmann constant
tmn = 0.2e-12;                          % mean time between collisions

% Region Defining
L = 200e-9;
W = 100e-9;

% Current Condition and variables
num = 1e4;                              % Number of electrons
T = 300;                                % Temperature (Kelvin)
V = 0.1;                                % Voltage applied
vth_e = sqrt((2*kb*T)/(meff));          % Thermal velocity of an electron
vth_ex = (vth_e/sqrt(2))*randn(num, 1); % X-component of thermal velocity
vth_ey = (vth_e/sqrt(2))*randn(num, 1); % Y-component of thermal velocity
vthdis = sqrt(vth_ex.^2+vth_ey.^2);     % Distribution of electrons thermal velocity
vthav = mean(sqrt(vth_ex.^2+vth_ey.^2));% Average of thermal velocity 
MFP = vthav*tmn;                        % Mean free path of electrons
Ex = V/L;                               % Electric field on x-axis
F = Ex*q_0;                             % Force applied to electrons
accel = F/meff;                         % Acceleration of electrons

% Electrons Defining
Elec = zeros(num, 4);
Elec(:, 1) = L*rand(num, 1);
Elec(:, 2) = W*rand(num, 1);
Elec(:, 3) = vth_ex;
Elec(:, 4) = vth_ey;
previous = zeros(num, 4);
previous = Elec;

% Electron simulation
figure(1);
t = 1e-11;                          % Total Time
dt = 1e-14;                         % Time Step
Psat = 1 - exp(-dt/tmn);            % Exponential Scattering Probability
numplot = 5;                        % Number of electron plotted
color = hsv(numplot);               % Colour Setup
    
for n = 0:dt:t 
    
    Elec(:, 3) = Elec(:, 3)+ accel*dt;
    
    if Psat > rand()
        vth_ex = (vth_e/sqrt(2))*randn(num, 1); 
        vth_ey = (vth_e/sqrt(2))*randn(num, 1);
        Elec(:, 3) = vth_ex;
        Elec(:, 4) = vth_ey;
    end
    
    for p = 1:1:num
        previous(p, 1) = Elec(p, 1);
        previous(p, 2) = Elec(p, 2);
        Elec(p, 1) = Elec(p, 1)+ Elec(p, 3)*dt;
        Elec(p, 2) = Elec(p, 2)+ Elec(p, 4)*dt;
    end
    
    % Plotting limited amount of electrons
    figure(1)
    for q = 1:1:numplot
        title('Electrons movement');
        plot([previous(q, 1), Elec(q, 1)], [previous(q, 2), Elec(q,2)], 'color', color(q, :))
        xlim([0 L])
        ylim([0 W])
        hold on
    end
    
    % Setting up the boundaries
    for o = 1:1:num
        % Looping on x-axis
        if Elec(o, 1) > L                       
            Elec(o, 1) = Elec(o, 1) - L;
            previous = Elec;
        end
        if Elec(o, 1) < 0
            Elec(o, 1) = Elec(o, 1) + L;
            previous = Elec;
        end
        % Reflecting on y-axis
        if Elec(o, 2) > W || Elec(o, 2) < 0
            Elec(o, 4) = -1*Elec(o, 4);
        end
    end
    
    % Plotting average temperature
    vaver = mean(sqrt(Elec(:, 3).^2 + Elec(:, 4).^2)); % Average thermal velocity
    aveT = (0.5*meff*vaver^2)/kb;              % Average temperature
    figure(2)
    scatter(n, aveT, 'r.')
    axis tight
    title('Average temperature');
    hold on
    
    % Plotting Current density
    figure(3)
    I = vaver*num*Ex*q_0;                      % Drift current of electron
    scatter(n, I, 'g.')
    axis tight
    title('Current density of electrons');
    hold on
    %pause(0.01)
end