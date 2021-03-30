function [x T Bndry] = fem_make_mesh_2D(NoElems1D, h, e_m);
%  fem_make_mesh_2D -- (FEM Tutorials)
%
%    This function creates a 2D mesh
%
%    Mandatory Inputs  --
%        NoElems1D  - Integer of the no of elements in each direction
%
%        h      - Real vector containing the size of the square domain in (x_1, x_2)
%
%        e_m    - the basis object for the desired mesh
%
%        -- Element Ordering --
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
%
%    Outputs --
%        x      - a column array (no global nodes x spatial dimension) storing 
%                 Global nodal coordinates for each node 
%
%        T      - an array (no elems x no local nodes) giving global node index 
%                 for element / local node index pairing 
%        
%        Bndry  - a 0 or 1 vector (no nodes x 1) labeling 
%                 interior nodes (0) and boundary nodes (1)
%
%  by Dr. D Nordsletten
%  Oct 2015
%

  % setting defaults (none) ...

  % Error checking ...
  fem_check_basis(e_m)
    % grabbing the polynomial order + no local basis functions...
    p = e_m.p;
    no_basis_func = size(e_m.y,2);

  fem_err((p  < 1),'Supplied polynomial order is outside the range of schemes (p > 0)')
  fem_err((p - floor(p)  > 0),'Supplied polynomial order must be an integer')
  for i = 1:2
    fem_err((h(i) <= 0),'Supplied mesh size must be greater than zero')
  end
  fem_err((NoElems1D < 2),'Supplied Domain must contain at least 2 elements in each coordinate direction')


  % Program beginning ...

  % Number of points ...
  Nnrow = NoElems1D * e_m.p + 1;
  N = Nnrow^2;

  % Number of elements ...
  Ne = NoElems1D^2;
  if(strcmpi(e_m.etype,'triangle'))
    Ne = 2 * Ne; % Each quadrilateral is split into two triangles ...
  end

  % Initializing arrays ...
  x = zeros(N,2);
  T = zeros(Ne, no_basis_func);
  Bndry = zeros(N,1);

  % xi points for nodes within the element ...
  xi = 0:(1/p):1;

  % Creating the Global nodes array and boundary array (x + Bndry) ...
  m = 0; % initializing counter ...
  Ni = h / p; % setting the nodal increment ...
  for j = 1:Nnrow % looping over the number of node rows ...
    for i = 1:Nnrow % looping over the nodes in each row ...
      m = m + 1; % incremementing counter ...
      x(m,1:2) = [ Ni(1) * (i-1); Ni(2) * (j-1) ]; % setting the coordinate values ...
      Bndry(m) = (j == 1) || (j == Nnrow) || (i == 1) || (i == Nnrow); % determining if the coordinate is on the boundary ...
    end
  end


  % Creating the Mesh Topology (T array) ...
  m = 0; % counter for nodes ...
  for ej = 1:NoElems1D
    for ei = 1:NoElems1D
      n = 1 + (ei-1) * p + (ej-1) * Nnrow * p; % incrementing node counter ...
      v = get_node_indexes(n, e_m, Nnrow); % local node 1:p+1 --> global node map
      for i = 1:size(v,1)
        m = m + 1;
        T(m,1:no_basis_func) = v(i,1:no_basis_func);
      end
    end
  end



function v = get_node_indexes(n, e_m, Nnrow);
  if(strcmpi(e_m.etype, 'quadrilateral'))
    v = zeros(1, size(e_m.y,2));
    if(e_m.p == 1)
      v(1,1:4) = [n  (n + 1)  (n+Nnrow) (n+Nnrow+1)];
    elseif(e_m.p == 2)
      v(1,1:9) = [n  (n + 2)  (n+2*Nnrow) (n+2*Nnrow+2) (n + 1) (n + Nnrow) (n+Nnrow+1) (n+Nnrow+2) (n+2*Nnrow+1)];
    else
      error('fem: unsupported polynomial order encountered')
    end

  elseif(strcmpi(e_m.etype, 'triangle'))
    v = zeros(2, size(e_m.y,2));
    if(e_m.p == 1)
      v(1,1:3) = [n  (n + 1)  (n+Nnrow)];
      v(2,1:3) = [(n+Nnrow+1) (n + 1)  (n+Nnrow)];
    elseif(e_m.p == 2)
      v(1,1:6) = [n  (n + 2)  (n+2*Nnrow) (n+1) (n+Nnrow) (n+Nnrow+1)];
      v(2,1:6) = [(n+2*Nnrow+2) (n + 2)  (n+2*Nnrow)  (n+Nnrow+2) (n+2*Nnrow+1) (n+Nnrow+1)];
    else
      error('fem: unsupported polynomial order encountered')
    end

  else
    error(['fem: unknown element type encountered' e_m.etype]);
  end

