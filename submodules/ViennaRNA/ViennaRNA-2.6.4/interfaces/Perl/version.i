%{
/** @file version.i.in
 * @brief Set $RNA::VERSION to the bindings version
 */
%}

%perlcode %{
our $VERSION = '2.6.4';
sub VERSION () { $VERSION };
%}

