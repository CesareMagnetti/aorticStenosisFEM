function triangle_lyness_rule_test ( )

%*****************************************************************************80
%
%% TRIANGLE_LYNESS_RULE_TEST tests the LYNESS_RULE library.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 September 2010
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TRIANGLE_LYNESS_RULE_TEST:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Test the TRIANGLE_LYNESS_RULE library.\n' );

  triangle_lyness_rule_test01 ( );
  triangle_lyness_rule_test02 ( );
  triangle_lyness_rule_test03 ( );
  triangle_lyness_rule_test04 ( );
  triangle_lyness_rule_test05 ( );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TRIANGLE_LYNESS_RULE_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );

  return
end