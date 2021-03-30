function [gpt gw] = fem_get_quadrule_2D(rule, etype);
%  fem_get_quadrule_2D -- (FEM Tutorials)
%
%    This function returns a Quadrature points for 2D master elements of type (etype) 
%
%    Mandatory Inputs  --
%        rule  - the rule of the returned quadrature scheme
%                for quadrilateral : rule = 1, 2, ... or 5
%                                     routine returns rule^2 Gaussian Quadrature points
%
%                for triangle : rule = 0, 1, 2, ... or 21
%                               routine returns the (rule)^{th} - Lyness Quadrature point scheme
%
%        etype  - the master element type (quadrilateral or triangle)
%
%    Outputs --
%        gpt    - an array (no. quad points by spatial dimension) with the Quadrature points 
%
%        gw     - an array (no. quad points) with the Quadrature weights 
%
%  by Dr. D Nordsletten
%  Oct 2014
%

  if(strcmpi(etype,'quadrilateral'))
    % Error checking ...
    fem_err( (rule < 1) || (rule > 5),'Supplied rule is outside the range of schemes [1,5]')

    % Getting 1D quadrature rule ...
    [gpt1d gw1d] = fem_get_quadrule_1D(rule);
    n = length(gw1d);

    gpt = zeros(n*n,2); gw = zeros(n*n,1);
    m = 0;
    for i = 1:n
      for j = 1:n
        m = m + 1;
        gpt(m,1:2) = [gpt1d(j), gpt1d(i)];
        gw(m) = gw1d(i) * gw1d(j);
      end
    end
    disp(['The selected rule provides accuracy up to O(' num2str(2*rule-1) ')'])

  elseif(strcmpi(etype,'triangle'))
    % Error checking ...
    fem_err( (rule < 0) || (rule > 21),'Supplied rule is outside the range of schemes [0,21]')

    % Getting the rule and rule numbers ...
    addpath('triangle_quadrature/')
    % Initializing rule index and number of points + finding the number of Lyness rules ...
    [gw gpt] = lyness_rule(rule, lyness_order(rule));
    gpt = gpt';
    gw =  0.5 * gw';
    disp(['The selected rule provides accuracy up to O(' num2str(lyness_precision(rule)) ')'])

  else
    error(['fem: unknown element type ' etype])
  end



