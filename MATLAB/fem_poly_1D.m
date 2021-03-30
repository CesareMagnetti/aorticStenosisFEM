function [val] = fem_poly_1D(xi, varargin);
%  fem_poly -- (FEM Tutorials)
%
%    This function generates a p^th order nodal Lagrange basis (default, p = 1), 
%    differentiates it diff times (default, diff=0), and returns the basis 
%    evaluated at all supplied xi-points.
%
%    Mandatory Inputs  --
%        xi     - the xi-points to compute at (singleton or column vector in [0,1])
%
%    Optional Inputs   --
%        p      - The nodal Lagrange polynomial order (integer > 0)
%
%        diff   - Number of times to differentiate the polynomial basis (integer >= 0)
%
%    Outputs --
%        val - a two-dimensional matrix (no xi-points by p+1)  with the nth-column giving the
%              evaluation of the n-th basis (ordered left to right)
%
%
%    Examples --
%
%        >>  y = fem_poly(0); 
%        % returns the 1-Nodal Lagrange basis functions evaluated at zero
%
%        >>  y = fem_poly([0; 0.1; 0.2], 2); 
%        % returns the 2nd order-Nodal Lagrange basis functions evaluated at [0; 0.1; 0.2]
%
%        >>  y = fem_poly([0; 1], 5, 2); 
%        % returns the 5th order-Nodal Lagrange basis functions, differentiated twice, evaluated at [0; 1]
%
%  by Dr. D Nordsletten
%  Jan 2015
%

% set defaults ...
p    = 1;
diff = 0;

% error checking ...
fem_err((max(xi) > 1) || (min(xi) < 0), 'Range of xi points supplied is outside the master element e_M = [0,1]');
fem_err(size(xi,2) > size(xi,1),'xi-points must be supplied as a column vector');
fem_err(nargin > 3, 'Number of arguments outside allowable range, see help');
if(nargin >= 2)
  p = varargin{1};
  fem_err(p < 1, 'Supplied polynomial order is outside the range allowed (p > 0)');
  fem_err(p - floor(p) > 0, 'Supplied polynomial order must be an integer');
  if(nargin == 3)
    diff = varargin{2};
    fem_err(diff < 0, 'Routine cannot integrate.');
    fem_err(diff - floor(diff) > 0, 'Supplied polynomial order must be an integer');
  end
end

% Initialize the returned array ...
val = zeros(size(xi,1), p+1);

% Getting the basis coefficients (a) and the power array (b) ...
[a  b] = load_coefficient_arrays(p);
  
% Computing the val - array ...
for base = 1:p+1 % computing for the basis (base) ...
   for pw = 1:p+1 % Evaluating each term and summing to the column val(:,base) ...
      ax = 1; % the coefficient due to differentiation (1 when diff = 0) ...
      for k = 0 : diff - 1 
         ax = ax * (b(pw) - k);
      end
      % Computing the pw^th term and adding to val(:,base) ...
      val(:,base) = val(:, base) + a(base,pw) * ax * (xi .^max(b(pw) - diff,0) );
   end
end
  
% Function which generates the coefficient / power arrays ..
function [a b] = load_coefficient_arrays(p);
  % constructing the node points ...
  x = [0:1/p:1]';

  % constructing the power array ...
  b = 0:p;

  % computing the vandermonde matrix ...
  for i = 1:p+1
    a(:,i) = x .^ b(i);
  end
  
  % computing the inverse transpose of the vandermonde matrix, i.e. the coefficient array ...
  a = inv(a)';


