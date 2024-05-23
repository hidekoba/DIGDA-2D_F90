module module_zchem
  implicit none
contains

! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

  function FK0(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK0

  FK0 = &
 &  EXP ( &
 &     - 162.8301D0 &
 &     + 218.2968D0 * 1.D2 / T &
 &     + 90.9241D0  * LOG( T / 1.D2 ) &
 &     - 1.47696D0  * ( T / 1.D2 )**2 &
 &     + S * (   0.025695D0 &
 &             - 0.025225D0  * ( T / 1.D2 ) &
 &             + 0.0049867D0 * ( T / 1.D2 )**2 ) &
 &  )

  end function FK0

  function FK1(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK1

  FK1 = &
 &  10.D0 ** ( -( &
 &     + 3670.7D0 / T &
 &     - 62.008D0 &
 &     + 9.7944D0 * LOG( T ) &
 &     - 0.0118D0 * S &
 &     + 0.000116D0 * S * S &
 &  ))

  end function FK1

  function FK2(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK2

  FK2 = &
 &  10.D0 ** ( -( &
 &     + 1394.7D0 / T &
 &     + 4.777D0 &
 &     - 0.0184D0 * S &
 &     + 0.000118D0 * S * S &
 &  ))

  end function FK2

  function FKB(T, S)

  double precision, intent(in) :: T, S
  double precision :: FKB

  FKB = &
 &  EXP ( &
 &     + ( - 8966.90D0 &
 &         - 2890.53D0 * SQRT( S ) &
 &         - 77.942D0  * S &
 &         + 1.728D0   * S * SQRT( S ) &
 &         - 0.0996D0  * S * S ) / T &
 &     + ( + 148.0248D0 &
 &         + 137.1942D0 * SQRT( S ) &
 &         + 1.62142D0  * S ) &
 &     + ( - 24.4344D0 &
 &         - 25.085D0 * SQRT( S ) &
 &         - 0.2474D0 * S ) * LOG( T ) &
 &     + 0.053105 * SQRT( S ) * T &
 &  )

  end function FKB

  function FKW(T, S)

  double precision, intent(in) :: T, S
  double precision :: FKW

  FKW = &
 &  EXP ( &
 &     - 13847.26D0 / T &
 &     + 148.9652D0 &
 &     - 23.6521D0 * LOG( T ) &
 &     + ( + 118.67D0 / T &
 &         - 5.977D0 &
 &         + 1.0495D0 * LOG( T ) ) * SQRT(S) &
 &     - 0.01615D0 * S &
 &  )

  end function FKW

  function FK1P(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK1P

  FK1P = &
 &  EXP ( &
 &     - 4576.752D0 / T &
 &     + 115.525D0 &
 &     - 18.453D0 * LOG( T ) &
 &     + ( - 106.736D0 / T &
 &         + 0.69171       ) / SQRT(S) &
 &     + ( - 0.65643D0 / T &
 &         - 0.01844       ) * S &
 &  )

  end function FK1P

  function FK2P(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK2P

  FK2P = &
 &  EXP ( &
 &     - 8814.15D0 / T &
 &     + 172.0883D0 &
 &     - 27.927D0 * LOG(T) &
 &     + ( - 160.340D0 / T &
 &         + 1.3566D0      ) / SQRT(S) &
 &     + ( + 0.37335D0 / T &
 &         - 0.05778D0     ) * S &
 &  )

  end function FK2P

  function FK3P(T, S)

  double precision, intent(in) :: T, S
  double precision :: FK3P

  FK3P = &
 &  EXP ( &
 &     - 3070.75D0 / T &
 &     - 18.141D0 &
 &     + ( + 17.27039D0 / T &
 &         + 2.81197D0      ) / SQRT(S) &
 &     + ( + 0.37335D0 / T &
 &         - 0.05778       ) * S &
 &  )

  end function FK3P

  function FKSI(T, S)

  double precision, intent(in) :: T, S
  double precision :: FKSI

  FKSI = &
 &  EXP ( &
 &     - 8904.2D0 / T &
 &     + 117.385D0 &
 &     - 19.334D0 * LOG(T) &
 &     + ( - 458.79D0 / T &
 &         + 3.5913D0     ) / SQRT(S) &
 &     + ( + 188.74D0 / T &
 &         - 1.5998D0     ) * &
 &       ( 19.924D0 * S &
 &         / ( 1000.D0 - 1.005D0 * S )) &
 &     + ( - 12.1652 / T &
 &         + 0.07871D0   ) * &
 &       ( 19.924D0 * S &
 &         / ( 1000.D0 - 1.005D0 * S ))**2 &
 &     + LOG( 1.0 - 0.001005D0 * S ) &
 &  )

  end function FKSI

  function FKSO4(T, S)

  double precision, intent(in) :: T, S
  double precision :: FKSO4

  FKSO4 = &
 &  EXP ( &
 &     - 4276.1D0 / T &
 &     + 141.328D0 &
 &     - 23.093D0 * LOG(T) &
 &     + ( - 13856.0D0 / T &
 &         + 324.57D0 &
 &         - 47.986D0 * LOG(T) ) * &
 &       SQRT( 19.924D0 * S / ( 1000.D0 - 1.005D0 * S )) &
 &     + ( + 35474.0D0 / T &
 &         - 771.54D0 &
 &         + 114.7231 * LOG(T) ) * &
 &       ( 19.924D0 * S / ( 1000.D0 - 1.005D0 * S )) &
 &     - 2698.0D0 / T * &
 &       ( 19.924D0 * S / ( 1000.D0 - 1.005D0 * S ))**1.5D0 &
 &     + 1776.0D0 / T * &
 &       ( 19.924D0 * S / ( 1000.D0 - 1.005D0 * S ))**2.0D0 &
 &     + LOG( 1.0D0 - 0.001005D0 * S ) &
 &  )

  end function FKSO4

  function FKF(T, S, XKS)

  double precision, intent(in) :: T, S, XKS
  double precision :: FKF

  FKF = &
 &  EXP ( &
 &     + 1590.2D0 / T &
 &     - 12.641D0 &
 &     + 1.525D0 / SQRT(S) &
 &     + LOG( 1.0D0 - 0.001005D0 * S ) &
 &     + LOG( 1.0D0 + ( 0.1400D0 / 96.062D0 )) * &
 &       ( S / 1.86655D0 / XKS ) &
 &  )

  end function FKF

  function FBT(S)

  double precision, intent(in) :: S
  double precision :: FBT

  FBT = &
 &     + 0.000232D0 / 10.811D0 &
 &       * S / 1.80655D0

  end function FBT

  function FSO4T(S)

  double precision, intent(in) :: S
  double precision :: FSO4T

  FSO4T = &
 &     + 0.14D0 / 96.062D0 &
 &       * S / 1.80655D0

  end function FSO4T

  function FFLT(S)

  double precision, intent(in) :: S
  double precision :: FFLT

  FFLT = &
 &     + 0.000067D0 / 18.9984D0 &
 &       * S / 1.80655D0

  end function FFLT

  function FALK( &
 &          XH, X1, X2, XB, XW, &
 &          X1P, X2P,  X3P, &
 &          XSI, XSO4, XF, &
 &          C, B, P, SI, SO4, FL )

  double precision, intent(in) :: XH, X1, X2, XB, XW, X1P, X2P, X3P
  double precision, intent(in) :: XSI, XSO4, XF
  double precision, intent(in) :: C, B, P, SI, SO4, FL
  double precision :: FALK

  FALK &
 &   =   C * ( XH*X1 + 2.D0*X1*X2 ) &
 &         / ( XH*XH + XH*X1 + X1*X2 ) &
 &     + B * XB / ( XH + XB ) &
 &     + XW / XH &
 &     - XH &
 &     + P * (  XH*X1P*X2P + 2.D0*X1P*X2P*X3P - XH*XH*XH ) &
 &         / (  XH*XH*XH + X1P*XH*XH + X1P*X2P*XH &
 &            + X1P*X2P*X3P ) &
 &     + SI * XSI / ( XH + XSI ) &
 &     + XH / ( 1.0D0 + SO4 / XSO4 ) &
 &     - SO4 / (  1.0D0 &
 &              + XSO4 / ( XH / ( 1.0D0 + SO4 / XSO4 ))) &
 &     - FL / ( 1.0D0 + XF / XH )

  end function FALK

end module module_zchem
