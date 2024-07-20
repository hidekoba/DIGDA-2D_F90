!######################################################################
module module_digda_parameter
  use module_digda_grid
  implicit none

  ! Density and Molar mass
  ! Density [g/cm^3], Molar mass [g/Mol], Molar mass / Density [cm^3/Mol]

  double precision, parameter :: rhocla = 2.60d0    !! 2.6 [g/cm^3]
  double precision, parameter :: rhopom = 2.60d0    !! 2.6 [g/cm^3]
  double precision, parameter :: rhocar = 2.60d0    !! 2.6 [g/cm^3]
  double precision, parameter :: rhoopa = 2.60d0    !! 2.6 [g/cm^3]

  double precision, parameter :: molcla = 430.0d0   !! 430 [g/Mol]
  double precision, parameter :: molpom = 12.0d0    !!  12 [g/Mol]
  double precision, parameter :: molcar = 100.0d0   !! 100 [g/Mol]
  double precision, parameter :: molopa = 60.0d0    !!  60 [g/Mol]

  double precision :: mcla, mpom, mcar, mopa

  double precision, parameter :: rho13pom = 2.60d0  !! 2.6 [g/cm^3]
  double precision, parameter :: rho14pom = 2.60d0  !! 2.6 [g/cm^3]
  double precision, parameter :: rho13car = 2.60d0  !! 2.6 [g/cm^3]
  double precision, parameter :: rho14car = 2.60d0  !! 2.6 [g/cm^3]

  double precision, parameter :: mol13pom = 13.0d0  !!  13 [g/Mol]
  double precision, parameter :: mol14pom = 14.0d0  !!  14 [g/Mol]
  double precision, parameter :: mol13car = 101.0d0 !! 101 [g/Mol]
  double precision, parameter :: mol14car = 102.0d0 !! 102 [g/Mol]

  double precision :: m13pom, m14pom, m13car, m14car

  ! Redfield ratio
  ! Fe:P, P:P, N:P, C:P, O:P
  double precision, parameter :: fp = 4.70d-4
  double precision, parameter :: pp = 1.0d0
  double precision, parameter :: np = 16.0d0
  double precision, parameter :: cp = 117.0d0
  double precision, parameter :: op = 170.0d0

  ! 13C, 14C
  double precision, parameter :: pdb = 0.0112372d0  !! standard 13C/12C
  double precision, parameter :: dstd = 1.1892d-12  !! standard 14C/12C
  double precision, parameter :: eoph = -22.0, eosh = 0.0  !! fractionation factor
  double precision, parameter :: eps = 1.0d-10
  double precision, parameter :: thalfyr = 5370.d0  !! half decay time [yr], decay constant

  ! Porosity
  double precision, parameter :: phi0 = 0.9d0  !! 0.9
  double precision :: z0 = 1.823d0  !! 1.823 [m]

  ! Bioturbation diffusivity
  double precision :: db0 = 3.0d-5 !! 3.0*10-5 [m^2/yr]
  double precision :: zb = 1.0d-1  !! 0.10 [m]

  ! Porewater effective molecular diffusivity
  double precision :: dsw = 4.59d-10      !! 4.59*10^-10 [m^2/s]
  double precision :: dswnit = 5.0d-10    !! 5.00*10^-10 [m^2/s]
  double precision :: dswdic = 6.4d-10    !! 6.40*10^-10 [m^2/s]
  double precision :: dswalk = 6.4d-10    !! 6.40*10^-10 [m^2/s]
  double precision :: dswdo2 = 1.2d-9     !! 1.20*10^-9  [m^2/s]
  double precision :: dswsil = 5.0d-10    !! 5.00*10^-10 [m^2/s]
  double precision :: dsw13dic = 6.4d-10  !! 6.40*10^-10 [m^2/s]
  double precision :: dsw14dic = 6.4d-10  !! 6.40*10^-10 [m^2/s]
  double precision, parameter :: nopa = 2.50d0

  ! Porosity
  double precision :: phi(sdz)
  double precision :: phii(sdz+1)

  ! Bioturbation diffusivity
  double precision :: dbi(sdz+1)

  ! Porewater diffusivity
  double precision :: deinit(sdz+1)
  double precision :: deidic(sdz+1)
  double precision :: deialk(sdz+1)
  double precision :: deido2(sdz+1)
  double precision :: deisil(sdz+1)
  double precision :: dei13dic(sdz+1)
  double precision :: dei14dic(sdz+1)

contains

  subroutine set_digda_parameter()
  end subroutine set_digda_parameter

end module module_digda_parameter
