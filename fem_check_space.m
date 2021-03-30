function [] = fem_check_space(space, Omegadm)
%  fem_check_space -- (FEM Tutorials)
%
%    This function checks that a space object is well defined.  
%
%    Mandatory Inputs  --
%
%        space              - the space object to check
%
%        Omegadm            - the spatial dimension of our domain
%
%    Optional Inputs   --
%
%    Outputs --
%
%  by David Nordsletten
%  May 2014
%

% Checking the space contains all necessary components ...
  fem_err(isfield(space,'name') == 0,'fem: there is no .name defined for space')
  fem_err(isfield(space,'dm')   == 0,['fem: there is no .dm (dimension) defined for space ' space.name])
  fem_err(isfield(space,'x')    == 0,['fem: there is no .x (mesh coordinates) defined for space ' space.name])
  fem_err(isfield(space,'t')    == 0,['fem: there is no .t (mesh topology) defined for space ' space.name])
  fem_err(isfield(space,'p')    == 0,['fem: there is no .p (polynomial order) defined for space ' space.name])
  fem_err(isfield(space,'e')    == 0,['fem: there is no .e (basis object) defined for space ' space.name])
  fem_check_basis(space.e);

% Checking to make sure things are sensible ...
  fem_err((1 > Omegadm)                , ['fem: spatial dimension of domain has dimension less than 1'])
  fem_err(Omegadm - floor(Omegadm) ~= 0, ['fem: spatial dimension of domain non-integer dimension'])

  fem_err((1 > space.dm)                 , ['fem: space ' space.name ' has dimension less than 1'])
  fem_err(space.dm - floor(space.dm) ~= 0, ['fem: space ' space.name ' has non-integer dimension'])
  
  fem_err(size(space.t,2) ~= size(space.e.y,2), ['fem: polynomial order of ' space.name ' inconsistent with topology'])

  nop= zeros(size(space.x,1),1);   
  for i = 1:size(space.t,1)
    for j = 1:size(space.t,2)
      fem_err(space.t(i,j) < 1, ['fem: topology of ' space.name ' has zero nodes'])
      fem_err(space.t(i,j) > size(space.x,1), ['fem: topology of ' space.name ' has nodes exceeding .x array'])
      nop(space.t(i,j)) = 1;
    end
  end
  fem_err(min(nop) == 0, ['fem: topology of ' space.name ' has nodes which do not appear in any elements'])

