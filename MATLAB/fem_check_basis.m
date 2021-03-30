function [] = check_element_object_for_errors(e);
%  check_element_object_for_errors -- (FEM Tutorials)
%
%    This function checks that the element object (e) is contains the requisite 
%    information.
%
%    Mandatory Inputs  --
%        e      - the element object
%
%    Outputs --
%        the function flags an error if the object is missing components or they 
%        are inappropriately formed.
%
%
%  by Dr. D Nordsletten
%  Jan 2015
%

% Checking that all relevant objects are present ...
fem_err(isfield(e,'gp') == 0,'fem: element object contains no gauss point array')
fem_err(isfield(e,'gw') == 0,'fem: element object contains no gauss weight array')
fem_err(isfield(e,'y') == 0 ,'fem: element object contains no precomputed basis')
fem_err(isfield(e,'dy') == 0,'fem: element object contains no precomputed basis derivatives')
fem_err(isfield(e,'p') == 0,'fem: element object contains polynomial degree of basis')

% Checking that the length of quadrature arrays is correct ...
fem_err(size(e.gp,1) ~= size(e.gw,1), 'fem: element object has incompatible quadrature arrays')
fem_err(size(e.gp,1) ~= size(e.y,1), 'fem: element object has incompatible quadrature arrays')
fem_err(size(e.gp,1) ~= size(e.dy,1), 'fem: element object has incompatible quadrature arrays')

if(strcmpi(e.etype,'quadrilateral')) || (strcmpi(e.etype,'triangle'))
else
  error('fem: element object is of unknown type (must be triangle or quadrilateral')
end


