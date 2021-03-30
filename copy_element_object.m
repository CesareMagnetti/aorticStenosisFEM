function [e_copy] = copy_element_object(e);
%  copy_element_object -- (FEM Tutorials)
%
%    This function copies the element object and makes a duplicate.
%
%    Mandatory Inputs  --
%        e      - the element object
%
%    Outputs --
%        e_copy - replicate of e.
%
%
%  by Dr. D Nordsletten
%  Jan 2015
%

fem_check_basis(e);

% copying object ...
e_copy.gp = e.gp;
e_copy.gw = e.gw;
e_copy.y  = e.y;
e_copy.dy  = e.dy;


