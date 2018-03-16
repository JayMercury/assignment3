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
vth_e = sqrt((kb*T)/(meff));            % Thermal velocity of an electron
vth_ex = (vth_e)*randn(num, 1);         % X-component of thermal velocity
vth_ey = (vth_e)*randn(num, 1);         % Y-component of thermal velocity
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
% previous = zeros(num, 4);
previous = Elec;

% Electron simulation
t = 1e-12;                          % Total Time
dt = 1e-14;                         % Time Step
Psat = 1 - exp(-dt/tmn);            % Exponential Scattering Probability
numplot = 5;                        % Number of electron plotted
color = hsv(numplot);               % Colour Setup
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
    
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
    set(0, 'CurrentFigure', f1)
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
%             previous = Elec;
        end
        if Elec(o, 1) < 0
            Elec(o, 1) = Elec(o, 1) + L;
%             previous = Elec;
        end
        % Reflecting on y-axis
        if Elec(o, 2) > W || Elec(o, 2) < 0
            Elec(o, 4) = -1*Elec(o, 4);
        end
    end
    
    % Plotting average temperature
    vaver = mean(sqrt(Elec(:, 3).^2 + Elec(:, 4).^2)); % Average thermal velocity
    aveT = (meff*vaver^2)/(kb);              % Average temperature
    set(0, 'CurrentFigure', f2)
    scatter(n, aveT, 'r.')
    axis tight
    title('Average temperature');
    hold on
    
    % Plotting Current density
    set(0, 'CurrentFigure', f3)
    I = vaver*num*Ex*q_0;                      % Drift current of electron
    scatter(n, I, 'g.')
    axis tight
    title('Current density of electrons');
    hold on
    pause(1e-7)
    
end

% Electron Density map
set(0, 'CurrentFigure', f4)
hist3(Elec(:, 1:2), [50 50]);
% Eden = hist3(Elec(:, 1:2), [50 50]);

% Temperature map
set(0, 'CurrentFigure', f5)
% Vend = sqrt(Elec(:, 3).^2 + Elec(:, 4).^2);
% Tend = (meff.*Vend.^2)./kb;
% hold on
% hist3(Elec(:, 1:2), [50 50]);
% Eden1 = Eden';
% Eden1(size(Eden, 1)+1, size(Eden, 2)+1) = 0;
% Xe = linspace(min(Elec(:, 1)), max(Elec(:, 1)), size(Eden, 1)+1);
% Ye = linspace(min(Elec(:, 2)), max(Elec(:, 2)), size(Eden, 1)+1);
% pcolor(Tend);
[binx, biny] = meshgrid(0:L/50:L, 0:W/50:W);
zcheck = zeros(51, 51);
tempcheck = zeros(51, 51);
counter = 0;
vtotal = 0;
for i = 1:50
    txmn = binx(1,i);
    txmx = binx(1, i+1);
    for r = 1:50
        tymn = biny(r, 1);
        tymx = biny(r+1, 1);
        for mm = 1:num
            if(Elec(mm,1)>txmn & Elec(mm,1)<txmx & Elec(mm,2)<tymx & Elec(mm,2)>tymn)
                counter = counter + 1;
                zcheck(i, r) = zcheck(i, r)+1;
                vtotal = vtotal + sqrt(Elec(mm, 3)^2+Elec(mm, 4)^2);
                if(counter ~= 0)
                    tempcheck(i,r) = meff*(vtotal^2)/(counter*kb);
                end
            end
        end
        vtotal = 0;
        counter = 0;
    end
end
surf(binx, biny,zcheck)
title("Temperature density of electrons")