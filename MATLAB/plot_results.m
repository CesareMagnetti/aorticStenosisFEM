clear; close all; clc;

%load final results
load('results0_6gauss.mat');

%load specific results for value of k
%k=0.95;
%load(sprintf('jetflow1/results%0.2f.mat',k));

%render all figures with chosen height and width to have optimal images for
%the report
figure('Renderer', 'painters', 'Position', [0 0 650 400]);
% use function fem_surface_plot given on keats
fem_surface_plot(pres, pres.u);
%edit the plot for better visualisation
ax = gca;
ax.FontSize = 20;
title("pressure field");
xlabel("X1 (mm)");
ylabel("X2 (mm)");
c = colorbar;
c.Label.String = 'pressure (in g/mm^2)';
c.Label.FontSize = 20;
%fixed cbar limits to allow comparison between different results (else
%coloring won't be comparable)
c.Limits = [-250, 850];
%save figure to file
%saveas(gcf,sprintf('/Users/cesaremagnetti/Desktop/FEM_project/pressure/jetflow1/p_%0.2f.png',k))


%render all figures with chosen height and width to have optimal images for
%the report
figure('Renderer', 'painters', 'Position', [0 0 650 400]); 
% use function fem_surface_plot given on keats
fem_surface_plot(vel, vel.u(1:2:end));
%edit the plot for better visualisation
ax = gca;
ax.FontSize = 20;
title("velocity, X1 component");
xlabel("X1 (mm)");
ylabel("X2 (mm)");
c = colorbar;
c.Label.String = 'velocity (in mm/s)';
c.Label.FontSize = 20;
%fixed cbar limits to allow comparison between different results (else
%coloring won't be comparable)
c.Limits = [-1250, 400];
%save figure to file
%saveas(gcf,sprintf('/Users/cesaremagnetti/Desktop/FEM_project/pressure/jetflow1/v1_%0.2f.png',k))


%render all figures with chosen height and width to have optimal images for
%the report
figure('Renderer', 'painters', 'Position', [0 0 650 400]);
% use function fem_surface_plot given on keats
fem_surface_plot(vel, vel.u(2:2:end));
%edit the plot for better visualisation
ax = gca;
ax.FontSize = 20;
title("velocity, X2 component");
xlabel("X1 (mm)");
ylabel("X2 (mm)");
c = colorbar;
c.Label.String = 'velocity (in mm/s)';
c.Label.FontSize = 20;
%fixed cbar limits to allow comparison between different results (else
%coloring won't be comparable)
c.Limits = [-550, 450];
%save figure to file
%saveas(gcf,sprintf('/Users/cesaremagnetti/Desktop/FEM_project/pressure/jetflow1/v2_%0.2f.png',k))