function [val] = fem_poly_2D(xi, etype, varargin);
%  fem_poly_2D -- (FEM Tutorials)
%
%    This function generates a p^th order nodal Lagrange basis (default, p = 1), 
%    differentiates it diff times (default, diff=0), and returns the basis 
%    evaluated at all supplied xi-points.
%
%    Mandatory Inputs  --
%        xi     - the xi-points to compute at (single value or array)
%
%        etype  - the master element type (quadrilateral or triangle)
%
%    Optional Inputs   --
%        p      - The nodal Lagrange polynomial order (integer > 0)
%
%        diff   - Number of times to differentiate the polynomial basis (integer array >= 0)
%
%    Outputs --
%        val - a two-dimensional matrix (no. xi-points by no. basis functions)
%
%            - Quadrilateral:
%                  (p=1)             (p=2)
%               3 ------- 4       3 -- 9 -- 4
%               |         |       |         |
%               |         |       6    7    8
%               |         |       |         |
%               1 ------- 2       1 -- 5 -- 2
%
%            - Triangle:
%                  (p=1)             (p=2)
%               3                 3
%               |  '              |  '
%               |    '            5    6
%               |      '          |       '
%               1 ------- 2       1 -- 4 -- 2
%
%  by Dr. D Nordsletten
%  Jan 2015
%

  % set defaults ...
  p    = 1;
  diff = [0 0];

  % error checking ...
  fem_err(size(xi,2) ~= 2, 'Size of xi point array is not 2')
  for i = 1:2
    fem_err((max(xi(:,i)) > 1) || (min(xi(:,i)) < 0), 'Range of xi points supplied is outside the master element e_M = [0,1]');
  end
  fem_err(nargin > 4, 'Number of arguments outside allowable range, see help');
  if(nargin >= 3)
    p = varargin{1};
    fem_err((p < 1) || (p > 2), 'Supplied polynomial order is outside the range allowed (p > 0)');
    fem_err(p - floor(p) > 0, 'Supplied polynomial order must be an integer');
    if(nargin == 4)
      diff = varargin{2};
      fem_err(size(diff,2) ~= 2, 'Size of diff vector is not 2')
      for i = 1:2
        fem_err(diff(i) < 0, 'Routine cannot integrate.');
        fem_err(diff(i) - floor(diff(i)) > 0, 'Supplied derivative order must be an integer');
      end
    end
  end

  if(strcmpi(etype,'quadrilateral'))
    valx = fem_poly_1D(xi(:,1), p, diff(1));
    valy = fem_poly_1D(xi(:,2), p, diff(2));

    if(p == 1)
      oneD = [1 1; 2 1; 1 2; 2 2];
    else
      oneD = [1 1; 3 1; 1 3; 3 3; 2 1; 1 2; 2 2; 3 2; 2 3];
    end

    % Initialize the returned array ...
    val = zeros(size(xi,1), (p+1)^2);

    for i = 1:size(xi,1)
      for j = 1:(p+1)^2
        val(i,j) = valx(i,oneD(j,1)) * valy(i,oneD(j,2));
      end
    end
  elseif(strcmpi(etype,'triangle'))
    % Getting the 
    if(p == 1)
      N = [1 -1 -1; 0 1 0; 0 0 1];      
      k = 3;
    elseif(p == 2)
      N = [1    -3    -3     2     4     2;
           0    -1     0     2     0     0;
           0     0    -1     0     0     2;
           0     4     0    -4    -4     0;
           0     0     4     0    -4    -4;
           0     0     0     0     4     0];
      k = 6;
    else
      error('fem: unsupported polynomial order encountered');
    end
    for i = 1:size(xi,1)
      v = triangle_get_xi_vector(xi(i,1:2),diff);
      val(i,:) = (N * v(1:k))';
    end

  else
    error(['fem: unknown element type encountered' etype]);
  end


function [v] = triangle_get_xi_vector(xi, diff)
  if(diff(1) == 0) && (diff(2) ==0)
    v = [1; xi(1); xi(2); xi(1)^2; xi(1)*xi(2); xi(2)^2];
  elseif(diff(1) == 1) && (diff(2) == 0)
    v = [0;     1;     0; 2*xi(1);       xi(2);       0];
  elseif(diff(1) == 0) && (diff(2) == 1)
    v = [0;     0;     1;       0;       xi(1); 2*xi(2)];
  else
    error('fem: unrecognized case for requested derivatives')
  end
