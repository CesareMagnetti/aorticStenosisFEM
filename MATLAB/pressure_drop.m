% load('pres0.mat');
% load('vel0.mat');
load('results3.mat');
%find values of pres on boundary 1
pres_G1 = pres.u(unique(reshape(pres.b(pres.b(:,end) == 1,2:end-1),[],1)));
%find values of pres on boundary 2
pres_G2 = pres.u(unique(reshape(pres.b(pres.b(:,end) == 2,2:end-1),[],1)));
%find coordinates of points on boundary 1
X_G1 = pres.x(unique(reshape(pres.b(pres.b(:,end) == 1,2:end-1),[],1)),:);
%find coordinates of points on boundary2
X_G2 = pres.x(unique(reshape(pres.b(pres.b(:,end) == 2,2:end-1),[],1)),:);

% evaluate integral of pressure over boundary 1 (trapezoidal rule was applied)
G1 = 0;
h1 = 0;
for i = 1:size(X_G1, 1)-1
    h = norm(X_G1(i+1,:)-X_G1(i,:)); %evaluate segment length
    h1 = h1+ h; %add segment length to total boundary length
    G1 = G1 + h*(pres_G1(i)+pres_G1(i+1)); %apply trapezium rule
end
G1 = G1/2;

% evaluate integral of pressure over boundary 2 (trapezoidal rule was applied)
G2 = 0;
h2 = 0;
for i = 1:size(X_G2, 1)-1
    h = norm(X_G2(i+1,:)-X_G2(i,:)); %evaluate segment length
    h2 = h2+ h; %add segment length
    G2 = G2 + h*(pres_G2(i)+pres_G2(i+1)); %apply rule
end
G2 = G2/2;

%evaluate mean value by dividing the integral for the total length of the boundary
%then evaluate pressure drop
pressure_drop1 = G1/h1-G2/h2;

%evaluate pressure drop simply using the mean of the pressures without integrating
pressure_drop2 = mean(pres_G1) - mean(pres_G2);