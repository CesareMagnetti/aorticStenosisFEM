function [Re] = darcy_momentum_eqn(e,testsp, teste, vars);
% darcy_momentum_eqn -- (FEM Tutorials)
%
% This function computes the element matrix for the 2D discrete Laplacian
%
% Inputs --
%
%        e       - the current element (integer)
%
%        testsp  - the space of test functions ...
%
%        teste   - the element / basis object for test functions 
%                  evaluated for the local element ...
%
%        vars    - a structure containing variables ...
%                  contains vars.[name]e - which contains basis functions 
%                                          evaluated for the local element
%
%
% Outputs --
%
%        Re - the element residual for Darcy momentum equation ...
%
% by David Nordsletten
% Jan 2015
%

  % Problem parameters
  K_H = 10;
  K_L = 1;
  % Getting the local nodal values of velocity and pressure, note that the
  % two components of velocity are interleaved
  vl(:,1) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
  vl(:,2) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 2);
  pl(:,1) = vars.pres.u (vars.pres.dm * (vars.pres.t(e,:)-1)+1);

  % assemble velocity and pressure arrays @ quadrature points
  % (weighted sum of weights and basis functions)
  v = [vars.vele.y(:,:) * vl(:,1) vars.vele.y(:,:) * vl(:,2)];
  p = vars.prese.y(:,:) * pl(:,1);


  % Computing the nonlinear permeability tensor @ quadrature points ...
  k = 0.5*(K_H*(1+atan(p)) + K_L*(1-atan(p)));

  % Getting local row / local column sizes ...
  ne=size(testsp.t,2);
  Re = zeros(testsp.dm * ne, 1);

  for i = 1:ne 
    % getting the local element residual indices ordered v_1^{n=1}, v_2^{n=1}, v_1^{n=2}, v_2^{n=2} ...
    vei = (testsp.dm * (i - 1) + 1):(testsp.dm * i);
    % computing k v . w - p div w ...
    Re(vei) = [dot(teste.gw, k.*v(:,1).*teste.y(:,i) - p.*(teste.dy(:,i,1)))...
               dot(teste.gw, k.*v(:,2).*teste.y(:,i) - p.*(teste.dy(:,i,2)))];
  end


