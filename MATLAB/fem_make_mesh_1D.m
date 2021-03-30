function [u x T] = fem_make_mesh_1D(Omega, ho, p);
%  fem_make_mesh_1D -- (FEM Tutorials)
%
%    This function discretizes the domain Omega with mesh size ho (order p)
%
%    Mandatory Inputs  --
%        Omega  - a range [a,b] \subset R, a < b
%
%        ho     - an estimated element size (real)
%                 (note: this should be adjusted to match range)
%
%        p      - desired order of interpolation (integer)
%
%
%    Outputs --
%        u      - a zero column to store your solution 
%
%        x      - a column vector storing nodal points 
%
%        T      - an array (no elems x p+1) giving global indices 
%                 for element / local node pairings 
%
%  by Dr. D Nordsletten
%  Feb 2014
%

% setting defaults (none) ...

% Error checking ...
fem_err((p  < 1),'Supplied polynomial order is outside the range of schemes (p > 0)')
fem_err((p - floor(p)  > 0),'Supplied polynomial order must be an integer')
fem_err((ho <= 0),'Supplied mesh size must be greater than zero')
fem_err((length(Omega)~= 2),'Supplied Domain must be a vector [a b], where a < b')
fem_err(Omega(2) <= Omega(1),'Supplied Domain must satisfy a < b')

% Program beginning ...

% Creating the size of our interval ...
Interval = Omega(2) - Omega(1);

% Number of points ...
N = ceil(Interval / ho) + 1;

% Number of elements ...
Ne = N - 1;

% Returned mesh size ...
h = Interval / Ne;

% Initializing arrays ...
u = zeros(p * Ne + 1,1);
x = zeros(p * Ne + 1,1);
T = zeros(Ne, p+1);

% xi points for nodes within the element ...
xi = 0:1/p:1;

% Looping over the elements and nodes ...
n = 1; % counter for nodes ...
for e = 1:Ne % counter for elements ...
  x(n:n+p) = h * xi(1:p+1) + h*(e-1) + Omega(1); % computing the node points
  T(e,1:p+1) = n:n+p; % local node 1:p+1 --> global node map
  n = n + p; % incrementing node counter ...
end

