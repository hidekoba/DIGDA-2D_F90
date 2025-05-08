module module_zchem
  implicit none
contains

! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

  function fk0(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk0

  fk0 = &
 &  exp ( &
 &     - 162.8301d0 &
 &     + 218.2968d0 * 1.d2 / t &
 &     + 90.9241d0  * log( t / 1.d2 ) &
 &     - 1.47696d0  * ( t / 1.d2 )**2 &
 &     + s * (   0.025695d0 &
 &             - 0.025225d0  * ( t / 1.d2 ) &
 &             + 0.0049867d0 * ( t / 1.d2 )**2 ) &
 &  )

  end function fk0

  function fk1(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk1

  fk1 = &
 &  10.d0 ** ( -( &
 &     + 3670.7d0 / t &
 &     - 62.008d0 &
 &     + 9.7944d0 * log( t ) &
 &     - 0.0118d0 * s &
 &     + 0.000116d0 * s * s &
 &  ))

  end function fk1

  function fk2(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk2

  fk2 = &
 &  10.d0 ** ( -( &
 &     + 1394.7d0 / t &
 &     + 4.777d0 &
 &     - 0.0184d0 * s &
 &     + 0.000118d0 * s * s &
 &  ))

  end function fk2

  function fkb(t, s)

  double precision, intent(in) :: t, s
  double precision :: fkb

  fkb = &
 &  exp ( &
 &     + ( - 8966.90d0 &
 &         - 2890.53d0 * sQRt( s ) &
 &         - 77.942d0  * s &
 &         + 1.728d0   * s * sQRt( s ) &
 &         - 0.0996d0  * s * s ) / t &
 &     + ( + 148.0248d0 &
 &         + 137.1942d0 * sQRt( s ) &
 &         + 1.62142d0  * s ) &
 &     + ( - 24.4344d0 &
 &         - 25.085d0 * sQRt( s ) &
 &         - 0.2474d0 * s ) * log( t ) &
 &     + 0.053105 * sQRt( s ) * t &
 &  )

  end function fkb

  function fkw(t, s)

  double precision, intent(in) :: t, s
  double precision :: fkw

  fkw = &
 &  exp ( &
 &     - 13847.26d0 / t &
 &     + 148.9652d0 &
 &     - 23.6521d0 * log( t ) &
 &     + ( + 118.67d0 / t &
 &         - 5.977d0 &
 &         + 1.0495d0 * log( t ) ) * sQRt(s) &
 &     - 0.01615d0 * s &
 &  )

  end function fkw

  function fk1p(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk1p

  fk1p = &
 &  exp ( &
 &     - 4576.752d0 / t &
 &     + 115.525d0 &
 &     - 18.453d0 * log( t ) &
 &     + ( - 106.736d0 / t &
 &         + 0.69171       ) / sQRt(s) &
 &     + ( - 0.65643d0 / t &
 &         - 0.01844       ) * s &
 &  )

  end function fk1p

  function fk2p(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk2p

  fk2p = &
 &  exp ( &
 &     - 8814.15d0 / t &
 &     + 172.0883d0 &
 &     - 27.927d0 * log(t) &
 &     + ( - 160.340d0 / t &
 &         + 1.3566d0      ) / sQRt(s) &
 &     + ( + 0.37335d0 / t &
 &         - 0.05778d0     ) * s &
 &  )

  end function fk2p

  function fk3p(t, s)

  double precision, intent(in) :: t, s
  double precision :: fk3p

  fk3p = &
 &  exp ( &
 &     - 3070.75d0 / t &
 &     - 18.141d0 &
 &     + ( + 17.27039d0 / t &
 &         + 2.81197d0      ) / sQRt(s) &
 &     + ( + 0.37335d0 / t &
 &         - 0.05778       ) * s &
 &  )

  end function fk3p

  function fksi(t, s)

  double precision, intent(in) :: t, s
  double precision :: fksi

  fksi = &
 &  exp ( &
 &     - 8904.2d0 / t &
 &     + 117.385d0 &
 &     - 19.334d0 * log(t) &
 &     + ( - 458.79d0 / t &
 &         + 3.5913d0     ) / sQRt(s) &
 &     + ( + 188.74d0 / t &
 &         - 1.5998d0     ) * &
 &       ( 19.924d0 * s &
 &         / ( 1000.d0 - 1.005d0 * s )) &
 &     + ( - 12.1652 / t &
 &         + 0.07871d0   ) * &
 &       ( 19.924d0 * s &
 &         / ( 1000.d0 - 1.005d0 * s ))**2 &
 &     + log( 1.0 - 0.001005d0 * s ) &
 &  )

  end function fksi

  function fkso4(t, s)

  double precision, intent(in) :: t, s
  double precision :: fkso4

  fkso4 = &
 &  exp ( &
 &     - 4276.1d0 / t &
 &     + 141.328d0 &
 &     - 23.093d0 * log(t) &
 &     + ( - 13856.0d0 / t &
 &         + 324.57d0 &
 &         - 47.986d0 * log(t) ) * &
 &       sQRt( 19.924d0 * s / ( 1000.d0 - 1.005d0 * s )) &
 &     + ( + 35474.0d0 / t &
 &         - 771.54d0 &
 &         + 114.7231 * log(t) ) * &
 &       ( 19.924d0 * s / ( 1000.d0 - 1.005d0 * s )) &
 &     - 2698.0d0 / t * &
 &       ( 19.924d0 * s / ( 1000.d0 - 1.005d0 * s ))**1.5d0 &
 &     + 1776.0d0 / t * &
 &       ( 19.924d0 * s / ( 1000.d0 - 1.005d0 * s ))**2.0d0 &
 &     + log( 1.0d0 - 0.001005d0 * s ) &
 &  )

  end function fkso4

  function fkf(t, s, xks)

  double precision, intent(in) :: t, s, xks
  double precision :: fkf

  fkf = &
 &  exp ( &
 &     + 1590.2d0 / t &
 &     - 12.641d0 &
 &     + 1.525d0 / sQRt(s) &
 &     + log( 1.0d0 - 0.001005d0 * s ) &
 &     + log( 1.0d0 + ( 0.1400d0 / 96.062d0 )) * &
 &       ( s / 1.86655d0 / xks ) &
 &  )

  end function fkf

  function fbt(s)

  double precision, intent(in) :: s
  double precision :: fbt

  fbt = &
 &     + 0.000232d0 / 10.811d0 &
 &       * s / 1.80655d0

  end function fbt

  function fso4t(s)

  double precision, intent(in) :: s
  double precision :: fso4t

  fso4t = &
 &     + 0.14d0 / 96.062d0 &
 &       * s / 1.80655d0

  end function fso4t

  function fflt(s)

  double precision, intent(in) :: s
  double precision :: fflt

  fflt = &
 &     + 0.000067d0 / 18.9984d0 &
 &       * s / 1.80655d0

  end function fflt

  function falk( &
 &          xh, x1, x2, xb, xw, &
 &          x1p, x2p,  x3p, &
 &          xsi, xso4, xf, &
 &          c, b, p, si, so4, fl )

  double precision, intent(in) :: xh, x1, x2, xb, xw, x1p, x2p, x3p
  double precision, intent(in) :: xsi, xso4, xf
  double precision, intent(in) :: c, b, p, si, so4, fl
  double precision :: falk

  falk &
 &   =   c * ( xh*x1 + 2.d0*x1*x2 ) &
 &         / ( xh*xh + xh*x1 + x1*x2 ) &
 &     + b * xb / ( xh + xb ) &
 &     + xw / xh &
 &     - xh &
 &     + p * (  xh*x1p*x2p + 2.d0*x1p*x2p*x3p - xh*xh*xh ) &
 &         / (  xh*xh*xh + x1p*xh*xh + x1p*x2p*xh &
 &            + x1p*x2p*x3p ) &
 &     + si * xsi / ( xh + xsi ) &
 &     + xh / ( 1.0d0 + so4 / xso4 ) &
 &     - so4 / (  1.0d0 &
 &              + xso4 / ( xh / ( 1.0d0 + so4 / xso4 ))) &
 &     - fl / ( 1.0d0 + xf / xh )

  end function falk

end module module_zchem
