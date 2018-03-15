%% Part 2(a)
% In part 2(a), we are now trying to model the current flow in the region
% with a bottle neck added.

% Reset Everything
close all
clear

% Setting variables
nx = 50;                % Length of the region
ny = nx*3/2;            % Width of the region, 3/2 of length
G = sparse(nx*ny);      % Initialize a G matrix
D = zeros(1, nx*ny);    % Initialize a matrix for G matrix operation
S = zeros(ny, nx);      % Initialize a matrix for sigma
sigma1 = .01;             % Setting up parameter of sigma in different region
sigma2 = 1;
box = [nx*2/5 nx*3/5 ny*2/5 ny*3/5]; % Setting up the bottle neck

%Sigma matrix setup
sigma = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        if i > box(1) && i < box(2) && (j < box(3)||j > box(4))
            sigma(i, j) = sigma1;
        else
            sigma(i, j) = sigma2;
        end
    end
end

% Implement the G matrix with the bottle neck condition in the region
for i = 1:nx
    for j = 1:ny
        
        n = j + (i-1)*ny;
        nip = j + (i+1-1)*ny;
        nim = j + (i-1-1)*ny;
        njp = j + 1 + (i-1)*ny;
        njm = j - 1 + (i-1)*ny;
        
        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            D(n) = 1;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
            D(n) = 0;
        elseif j == 1
            G(n, nip) = (sigma(i+1, j) + sigma(i,j))/2;
            G(n, nim) = (sigma(i-1, j) + sigma(i,j))/2;
            G(n, njp) = (sigma(i, j+1) + sigma(i,j))/2;            
            G(n, n) = -(G(n,nip)+G(n,nim)+G(n,njp));
        elseif j == ny
            G(n, nip) = (sigma(i+1, j) + sigma(i,j))/2;
            G(n, nim) = (sigma(i-1, j) + sigma(i,j))/2;
            G(n, njm) = (sigma(i, j-1) + sigma(i,j))/2;
            G(n, n) = -(G(n,nip)+G(n,nim)+G(n,njm));
        else
            G(n, nip) = (sigma(i+1, j) + sigma(i,j))/2;
            G(n, nim) = (sigma(i-1, j) + sigma(i,j))/2;
            G(n, njp) = (sigma(i, j+1) + sigma(i,j))/2;
            G(n, njm) = (sigma(i, j-1) + sigma(i,j))/2;
            G(n, n) = -(G(n,nip)+G(n,nim)+G(n,njp)+G(n,njm));
        end
    end
end

% % Creating a surface plot for sigma
% figure(1)
% surf(sigma);
% axis tight
% view([40 30]);
% title("Surface plot of sigma")

% Calculating the voltage
V = G\D';

% Inverting the G matrix
X = zeros(ny, nx, 1);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        X(j,i) = V(n);
    end
end

% Creating a surface plot for voltage
figure(1)
surf(X)
axis tight
view([40 30]);
title("Surface plot of voltage with bottle neck condition")

% Calculating the electric field from voltage
[Ex, Ey] = gradient(X);
% 
% % Creating surface plots for x and y component for electric field
% figure(3)
% surf(-Ex)
% axis tight
% view([40 30]);
% title("Surface plot of x-component of electric field")
% 
% figure(4)
% surf(-Ey)
% axis tight
% view([40 30]);
% title("Surface plot of y-component of electric field")
% 
% % Calculating the current density
% Jx = sigma'.*-Ex;
% Jy = sigma'.*-Ey;
% J = sqrt(Jx.^2 + Jy.^2);
% 
% % Creating a surface plot for the current density
% figure(5)
% surf(J)
% axis tight
% view([40 30]);
% title("Surface plot of current density")

% Vector plot of electric field
figure(2)
quiver(-Ex, -Ey);
axis tight