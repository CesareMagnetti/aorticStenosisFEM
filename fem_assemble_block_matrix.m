function [A] = fem_assemble_block_matrix(elemat, Omega, testsp, trialsp);
%  fem_assemble_block_matrix -- (FEM Tutorials)
%
%    This function assembles an FEM block matrix defined by elemat using
%    the test space (testsp) and trial space (trialsp);
%
%    Mandatory Inputs  --
%
%        elemat             - an element matrix function with prototype
%                             elemat(testsp, trialsp, teste, triale);
%                                 teste   - the test functions element object
%                                 triale  - the trial functions element object
%        
%        Omega              - the variable space to compute over
%
%        testsp             - the test space to use
%
%        trialsp            - the trial space to use
%
%
%    Outputs --
%        A      - a block matrix of elemat
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

  %calculating the global matrix (A) size and element matrix (Ae) size ...
  n=testsp.dm * size(testsp.x,1);
  m=trialsp.dm * size(trialsp.x,1);

  ne=testsp.dm * size(testsp.t,2);
  me=trialsp.dm * size(trialsp.t,2);

  % Initializing the matrix entries ...
  A=sparse(n,m);
  Ae=zeros(ne,me);

  % Looping over all the elements ...
  for e = 1:size(Omega.t,1)
    % mapping basis function quantities to the element ...
    [metrics] = fem_compute_metrics(Omega, e);
    [teste]  = map_element_quantities(metrics, testsp.e);
    [triale] = map_element_quantities(metrics, trialsp.e);

    % Computing the element matrix ...
    Ae = elemat(testsp, trialsp, teste, triale);

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
