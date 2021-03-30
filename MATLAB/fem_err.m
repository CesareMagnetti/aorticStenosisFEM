function [] = fem_err(err,msg)
%  fem_err -- (FEM Tutorials)
%
%    This function checks an error code (err), and throws an error message (msg)
%
%    Mandatory Inputs  --
%        err     - the error code (0 - no error, 1 - an error)
%        msg     - a message to return to prompt
%
%  by Dr. D Nordsletten
%  Feb 2014
%
if(err)
  error(msg)
end