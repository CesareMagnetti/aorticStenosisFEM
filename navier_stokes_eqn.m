function [Re] = navier_stokes_eqn(e,testsp, teste, vars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% navier_stokes_eqn: steady flow --> 
%       STRONG FORM: row*V*grad(V) - mu*div(grad(V)) + grad(P)
%       WEAK FORM: row*V*grad(V)*w + mu*grad(V)**grad(w) - P*div(w)
%
% This function computes residual vector for navier stokes with steady flow
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
% by Cesare Magnetti
% Mar 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Problem parameters
  mu = 4e-2;
  ro = 1e-3;
  % Getting the local nodal values of velocity and pressure, note that the
  % two components of velocity are interleaved
  vl(:,1) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
  vl(:,2) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 2);
  pl(:,1) = vars.pres.u (vars.pres.dm * (vars.pres.t(e,:)-1)+1);

  % assemble velocity and pressure arrays @ quadrature points
  % (weighted sum of weights and basis functions)
  v = [vars.vele.y(:,:) * vl(:,1) vars.vele.y(:,:) * vl(:,2)];
  p = vars.prese.y(:,:) * pl(:,1);
  % assemble velocity gradient tensor @ quadrature points
  % (weighted sum of weights and basis functions)
  
  %this code should work no matter the dimensionality of the problem
  gradv = zeros(size(vars.vele.dy(:,:,1),1),size(vl,2),size(vl,2));
  for a = 1:size(v,2)
      for b = 1:size(v,2)
         gradv(:,a,b) = vars.vele.dy(:,:,b) * vl(:,a);
      end
  end
 
  % Getting local row / local column sizes ...
  ne=size(testsp.t,2);
  Re = zeros(testsp.dm * ne, 1);

  for i = 1:ne 
    % getting the local element residual indices ordered v_1^{n=1}, v_2^{n=1}, v_1^{n=2}, v_2^{n=2} ...
    vei = (testsp.dm * (i - 1) + 1):(testsp.dm * i);
%     % computing row*V*grad(V)*w + mu*grad(V)**grad(w) - P*div(w) ...
%     Re(vei) = [dot(teste.gw, ro*v(:,1).*(gradv(:,1,1)+gradv(:,1,2)).* teste.y(:,i) +...
%                              mu*(gradv(:,1,1).*teste.dy(:,i,1) + gradv(:,1,2).*teste.dy(:,i,2)) -...
%                              p.*(teste.dy(:,i,1))),...
%                dot(teste.gw, ro*v(:,2).*(gradv(:,2,1)+gradv(:,2,2)).* teste.y(:,i) +...
%                              mu*(gradv(:,2,1).*teste.dy(:,i,1) + gradv(:,2,2).*teste.dy(:,i,2)) -...
%                              p.*(teste.dy(:,i,2)))];
%                          
    Re(vei) = [dot(teste.gw, ro*(v(:,1).*gradv(:,1,1)+v(:,2).*gradv(:,1,2)).* teste.y(:,i) +...
                             mu*(gradv(:,1,1).*teste.dy(:,i,1) + gradv(:,1,2).*teste.dy(:,i,2)) -...
                             p.*(teste.dy(:,i,1))),...
               dot(teste.gw, ro*(v(:,1).*gradv(:,2,1)+v(:,2).*gradv(:,2,2)).* teste.y(:,i) +...
                             mu*(gradv(:,2,1).*teste.dy(:,i,1) + gradv(:,2,2).*teste.dy(:,i,2)) -...
                             p.*(teste.dy(:,i,2)))];
  end


