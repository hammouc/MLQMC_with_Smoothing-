function k = i4_xor ( i, j )

%*****************************************************************************80
%
%% i4_xor() calculates the exclusive OR of two integers.
%
%  Discussion:
%
%    MATLAB provides a BITXOR ( i, j ) function which should be used instead!
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2005
%
%  Author:
%
%   John Burkardt
%
%  Input:
%
%    integer I, J, two values whose exclusive OR is needed.
%
%  Output:
%
%    integer K, the exclusive OR of I and J.
%
  k = 0;
  l = 1;

  i = floor ( i );
  j = floor ( j );

  while ( i ~= 0 || j ~= 0 )
%
%  Check the current right-hand bits of I and J.
%  If they differ, set the appropriate bit of K.
%
    i2 = floor ( i / 2 );
    j2 = floor ( j / 2 );

    if ( ...
      ( ( i == 2 * i2 ) && ( j ~= 2 * j2 ) ) || ...
      ( ( i ~= 2 * i2 ) && ( j == 2 * j2 ) ) )
      k = k + l;
    end

    i = i2;
    j = j2;
    l = 2 * l;

  end

  return
end
