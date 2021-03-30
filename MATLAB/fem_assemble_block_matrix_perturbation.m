function [A] = fem_assemble_block_matrix_perturbation(eleres, Omega, testsp, trialsp, vars);
%  fem_assemble_block_matrix_perturbation -- (FEM Tutorials)
%
%    This function assembles an FEM block matrix defined by eleres using
%    the test space (testsp) and trial space (trialsp);
%
%    Mandatory Inputs  --
%
%        eleres             - an element residual function with prototype
%                             eleres(e, teste, vars);
%                                 e       - the element index e (integer)
%                                 teste   - the test functions element object
%                                 vars    - a list of variables
%        
%        Omega              - the variable space to compute over
%
%        testsp             - the test space to use
%
%        trialsp            - the trial space to use
%
%        vars               - the list of all problem variables
%
%
%    Outputs --
%        A      - a block matrix of D R / dU
%
%  by David Nordsletten
%  Jan 2015
%

  % Error checking ...
  fem_check_space(Omega, Omega.dm)
  fem_check_space(testsp, Omega.dm)
  fem_check_space(trialsp, Omega.dm)

  fem_err( size(Omega.t,1) ~= size(testsp.t,1) , 'fem: the number of elements in the domain inconsistent with test space')
  fem_err( size(Omega.t,1) ~= size(trialsp.t,1), 'fem: the number of elements in the domain inconsistent with trial space')

  try
    names = fieldnames(vars);
    n_names = length(names);
    pertid = 0;
    for i = 1:n_names
      if(strcmpi(names{i},trialsp.name))
        pertid = i;
      end
    end
    fem_err(pertid == 0, 'fem: vas must contain the trialsp variable')
  catch
    error('fem: the structure vars must contain the trialsp variable.') 
  end

  %calculating the global matrix (A) size and element matrix (Ae) size ...
  n=testsp.dm * size(testsp.x,1);
  m=trialsp.dm * size(trialsp.x,1);

  ne=testsp.dm * size(testsp.t,2);
  me=trialsp.dm * size(trialsp.t,2);

  % Initializing the matrix entries ...
  A=sparse(n,m);
  Ae=zeros(ne,me);

  % Setting the perturbation parameter ...
  eps = 1e-4;

  % Looping over all the elements ...
  for e = 1:size(Omega.t,1)
    % mapping basis function quantities to the element ...
    [metrics] = fem_compute_metrics(Omega, e);
    [teste]  = map_element_quantities(metrics, testsp.e);

    % Mapping derivatives for different variables .**  and storing as .**e ...
    for i = 1:n_names
      eval(['vars.' names{i} 'e = map_element_quantities(metrics, vars.' names{i} '.e);']);
    end
    
    % Computing the element matrix ...
    for i = 1:size(trialsp.t,2) 
      for j = 1:trialsp.dm 
        nn = trialsp.dm * (trialsp.t(e,i)-1) + j;
        eval(['vars.' names{pertid} '.u(nn) = vars.' names{pertid} '.u(nn) + eps;']);
        Rp = eleres(e, testsp, teste, vars);
        eval(['vars.' names{pertid} '.u(nn) = vars.' names{pertid} '.u(nn) - 2 * eps;']);
        Rm = eleres(e, testsp, teste, vars);
        k = trialsp.dm * (i-1) + j;
        Ae(:,k) = (Rp - Rm)./(2*eps);
      end
    end

    % Adding element matrix contributions to the global matrix A ...
    for i = 1:testsp.dm 
      for j = 1:trialsp.dm  
        vi  = testsp.dm * (testsp.t(e,:)-1) + i;
        vj  = trialsp.dm * (trialsp.t(e,:)-1) + j;
        vei = 1:size(testsp.t,2); vei = testsp.dm * (vei' - 1) + i;
        vej = 1:size(trialsp.t,2); vej = trialsp.dm * (vej' - 1) + j;
        A(vi, vj) = A(vi, vj) + Ae(vei, vej);
      end
    end
  end

%complete  

function [metrics] = fem_compute_metrics(Omega, elem);
  % extracting the local weights for the element mapping of coordinate x and coordinate y ...
  x = Omega.x(Omega.t(elem,:),1); % x - a vector of length (no. local nodes) storing the weights for the space map, p.
  y = Omega.x(Omega.t(elem,:),2); % y - a vector of length (no. local nodes) storing the weights for the space map, p.

  % this step uses interpolation to compute the spatial coordinate of each Gauss point ...
  metrics.x = Omega.e.y(:,:) * x; % metrics.x - a vector of length (no. quad points) storing x position (i.e. x_g)
  metrics.y = Omega.e.y(:,:) * y; % metrics.y - a vector of length (no. quad points) storing y position (i.e. y_g)

  % this step uses interpolation to compute the derivative of spatial coordinate map at each Gauss point ...
  dxdxi1 = Omega.e.dy(:,:,1) * x; % dxdxi1 - a vector of length (no. quad points) storing dx / d xi_1 
  dxdxi2 = Omega.e.dy(:,:,2) * x; % dxdxi2 - a vector of length (no. quad points) storing dx / d xi_2 
  dydxi1 = Omega.e.dy(:,:,1) * y; % dydxi1 - a vector of length (no. quad points) storing dy / d xi_1 
  dydxi2 = Omega.e.dy(:,:,2) * y; % dydxi2 - a vector of length (no. quad points) storing dy / d xi_2 

  % computing the jacobian of the mapping ...
  metrics.j = abs(dxdxi1 .* dydxi2 - dxdxi2 .* dydxi1); % metrics.j - a vector of length (no. quad points) storing the jacobian

  % computing the transform s.t. Grad(phi,x) = [metrics.m11 metrics.m12; metrics.m21 metrics.m22] Grad(phi,xi) ...
  metrics.m11 =  dydxi2 ./ metrics.j; 
  metrics.m12 = -dxdxi2 ./ metrics.j;
  metrics.m21 = -dydxi1 ./ metrics.j;
  metrics.m22 =  dxdxi1 ./ metrics.j;


function [teste] = map_element_quantities(metrics, e);
  teste = e;
  teste.gw = e.gw .* metrics.j;
  
  for i = 1:size(e.y,2)
    teste.dy(:,i,1) = metrics.m11 .* e.dy(:,i,1) + metrics.m21 .* e.dy(:,i,2);
    teste.dy(:,i,2) = metrics.m12 .* e.dy(:,i,1) + metrics.m22 .* e.dy(:,i,2);
  end
