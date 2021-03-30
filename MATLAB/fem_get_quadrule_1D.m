function [gpt gw] = fem_get_quadrule_1D(order,varargin);
%  fem_get_quadrule_1D -- (FEM Tutorials)
%
%    This function returns a Gaussian quadrature rule
%
%    Mandatory Inputs  --
%        order  - the order [1,5] of the returned quadrature scheme
%
%    Optional Inputs   --
%        range  - A vector giving the desired range of integration [a,b]
%
%
%    Outputs --
%        gpt    - a column vector with the quadrature points 
%
%        gw     - a column vector with the quadrature weights 
%
%  by Dr. D Nordsletten
%  Feb 2014
%

% setting defaults ..
InputRange = 0;

% Error checking ...
fem_err( (order < 1) || (order > 5),'Supplied order is outside the range of schemes [1,5]')
fem_err( nargin > 2, 'Number of supplied inputs is inconsistent with that allowed, see help.')
if(nargin == 2)
  range = varargin{1}; InputRange = 1;
  fem_err(length(range) ~= 2, 'Supplied range exceeds required dimensions (2)')
end

% Checking the order, and assigning the appropriate quadrature scheme ...
if(order == 1)
  gpt = [ 0.5 ];
  gw  = [ 1   ];
elseif(order == 2)
  a = 0.5*sqrt(1/3);
  gpt = [ (0.5 - a)  (0.5 + a) ];
  gw  = [ 1      1 ] / 2;
elseif(order == 3)
  a = 0.5 * sqrt(3/5);
  gpt = [ 0.5   (0.5-a)  (0.5+a) ];
  gw  = [ 8      5      5] / 18;
elseif(order == 4)
  a = 0.5 * sqrt( (1/7)*(3 - 2 *sqrt(6/5)));
  b = 0.5 * sqrt( (1/7)*(3 + 2 *sqrt(6/5)));
  gpt = [ (0.5 - a)  (0.5 + a)  (0.5 - b)  (0.5 + b) ];
  gw  = [ 18+sqrt(30)      18+sqrt(30)      18-sqrt(30)      18-sqrt(30)] / 72;
elseif(order == 5)
  a = (1/6) * sqrt(5 - 2 * sqrt(10/7));
  b = (1/6) * sqrt(5 + 2 * sqrt(10/7));
  gpt = [ (0.5 - b)    (0.5-a)  0.5   (0.5 + a) (0.5 + b) ];
  gw  = [ (322-13*sqrt(70))  (322+13*sqrt(70))  (115200/225)   (322+13*sqrt(70))  (322-13*sqrt(70))  ] / 1800;  
end

% adjusting to range ...
if(InputRange)
  gpt = (range(2) - range(1)) * gpt + range(1);  
  gw  = (range(2) - range(1)) * gw;
end

% return column vectors ...
gpt = gpt'; gw = gw';

