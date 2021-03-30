function [e] = fem_get_basis(p, quad_rule, etype);
%  fem_get_basis -- (FEM Tutorials)
%
%    This function generates a p^th order nodal Lagrange basis (default, p = 1), 
%    differentiates it diff times (default, diff=0), and returns the basis 
%    evaluated at all supplied xi-points.
%
%    Mandatory Inputs  --
%        p          - The nodal Lagrange polynomial order (integer > 0)
%
%        quad_rule - The quadrature rule order required
%
%        etype - Character variable 
%
%    Outputs --
%        e          - the basis object
%
%  by Dr. D Nordsletten
%  Jan 2015
%

  % set defaults ...

  % error checking ...
  fem_err(p <= 0, 'Polynomial order must be positive')
  fem_err(p - floor(p) ~= 0, 'Polynomial order must be an integer')

  fem_err(quad_rule <= 0, 'Quadrature order must be positive')
  fem_err(quad_rule - floor(quad_rule) ~= 0, 'Quadrature order must be an integer')

  fem_err(1-strcmpi(etype,'triangle')-strcmpi(etype,'quadrilateral'), 'Element type must be triangle or quadrilateral')
  
  e.etype = etype
  e.p = p;
  [e.gp, e.gw]=fem_get_quadrule_2D(quad_rule, etype); 
  [e.y]=fem_poly_2D(e.gp, etype, p);
  [e.dy(:,:,1)]=fem_poly_2D(e.gp, etype, p, [1 0]);
  [e.dy(:,:,2)]=fem_poly_2D(e.gp, etype, p, [0 1]);

  [gp gw] = fem_get_quadrule_1D(5);
  if(strcmpi(etype,'triangle'))
    % order of faces : face 1 - between nodes 1 & 2, 
    %                  face 2 - between nodes 1 & 3,
    %                  face 3 - between nodes 2 & 3,
    id = [ 1 2; 1 3; 2 3];
    e.normal = [-1  0  sqrt(2);
                 0 -1  sqrt(2)];
    e.x = [0 0; 1 0; 0 1; 0.5 0; 0 0.5; 0.5 0.5];
  else
    % order of faces : face 1 - between nodes 1 & 2, 
    %                  face 2 - between nodes 1 & 3,
    %                  face 3 - between nodes 2 & 4,
    %                  face 4 - between nodes 3 & 4,
    id = [ 1 2; 1 3; 2 4; 3 4];
    e.normal = [-1  0  0  1;
                 0 -1  1  0];
    e.x = [0 0; 1 0; 0 1; 1 1; 0.5 0; 0 0.5; 0.5 0.5; 1 0.5; 0.5 1];
  end
  
  % making face array ...
  for i = 1:size(id,1)
     % Getting Gauss weights and scaling by the line length ...
     e.f{i}.gw=gw * norm(e.x(id(i,1),:) - e.x(id(i,2),:),2);
     % Getting Gauss points and scaling by the line length ...
     e.f{i}.gp(:,1) = gp * (e.x(id(i,2),1) - e.x(id(i,1),1)) + e.x(id(i,1),1);
     e.f{i}.gp(:,2) = gp * (e.x(id(i,2),2) - e.x(id(i,1),2)) + e.x(id(i,1),2);
     % Getting basis functions + basis functions derivatives at the face i ...
    [e.f{i}.y]=fem_poly_2D(e.f{i}.gp, etype, p);
    [e.f{i}.dy(:,:,1)]=fem_poly_2D(e.f{i}.gp, etype, p, [1 0]);
    [e.f{i}.dy(:,:,2)]=fem_poly_2D(e.f{i}.gp, etype, p, [0 1]);
  end

  fem_check_basis(e)

