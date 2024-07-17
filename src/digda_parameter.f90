!######################################################################
module module_digda_parameter
  use module_digda_grid
  implicit none

  ! Density and Molar mass
  ! Density [g/cm^3], Molar mass [g/Mol], Molar mass / Density [cm^3/Mol]

  double precision, parameter :: RHOCLA = 2.60D0    !! 2.6 [g/cm^3]
  double precision, parameter :: RHOPOM = 2.60D0    !! 2.6 [g/cm^3]
  double precision, parameter :: RHOCAR = 2.60D0    !! 2.6 [g/cm^3]
  double precision, parameter :: RHOOPA = 2.60D0    !! 2.6 [g/cm^3]

  double precision, parameter :: MOLCLA = 430.0D0   !! 430 [g/Mol]
  double precision, parameter :: MOLPOM = 12.0D0    !!  12 [g/Mol]
  double precision, parameter :: MOLCAR = 100.0D0   !! 100 [g/Mol]
  double precision, parameter :: MOLOPA = 60.0D0    !!  60 [g/Mol]

  double precision :: MCLA, MPOM, MCAR, MOPA

  double precision, parameter :: RHO13POM = 2.60D0  !! 2.6 [g/cm^3]
  double precision, parameter :: RHO14POM = 2.60D0  !! 2.6 [g/cm^3]
  double precision, parameter :: RHO13CAR = 2.60D0  !! 2.6 [g/cm^3]
  double precision, parameter :: RHO14CAR = 2.60D0  !! 2.6 [g/cm^3]

  double precision, parameter :: MOL13POM = 13.0D0  !!  13 [g/Mol]
  double precision, parameter :: MOL14POM = 14.0D0  !!  14 [g/Mol]
  double precision, parameter :: MOL13CAR = 101.0D0 !! 101 [g/Mol]
  double precision, parameter :: MOL14CAR = 102.0D0 !! 102 [g/Mol]

  double precision :: M13POM, M14POM, M13CAR, M14CAR

  ! Redfield ratio
  ! Fe:P, P:P, N:P, C:P, O:P
  double precision, parameter :: FP = 4.70D-4
  double precision, parameter :: PP = 1.0D0
  double precision, parameter :: NP = 16.0D0
  double precision, parameter :: CP = 117.0D0
  double precision, parameter :: OP = 170.0D0

  ! 13C, 14C
  double precision, parameter :: PDB = 0.0112372D0  !! standard 13C/12C
  double precision, parameter :: DSTD = 1.1892D-12  !! standard 14C/12C
  double precision, parameter :: EOPH = -22.0, EOSH = 0.0  !! fractionation factor
  double precision, parameter :: EPS = 1.0D-10
  double precision, parameter :: THALFYR = 5370.D0  !! half decay time [yr], decay constant

  ! Porosity
  double precision, parameter :: PHI0 = 0.9D0  !! 0.9
  double precision :: Z0 = 1.823D0  !! 1.823 [m]

  ! Bioturbation diffusivity
  double precision :: DB0 = 3.0D-5 !! 3.0*10-5 [m^2/yr]
  double precision :: ZB = 1.0D-1  !! 0.10 [m]

  ! Porewater effective molecular diffusivity
  double precision :: DSW = 4.59D-10      !! 4.59*10^-10 [m^2/s]
  double precision :: DSWNIT = 5.0D-10    !! 5.00*10^-10 [m^2/s]
  double precision :: DSWDIC = 6.4D-10    !! 6.40*10^-10 [m^2/s]
  double precision :: DSWALK = 6.4D-10    !! 6.40*10^-10 [m^2/s]
  double precision :: DSWDO2 = 1.2D-9     !! 1.20*10^-9  [m^2/s]
  double precision :: DSWSIL = 5.0D-10    !! 5.00*10^-10 [m^2/s]
  double precision :: DSW13DIC = 6.4D-10  !! 6.40*10^-10 [m^2/s]
  double precision :: DSW14DIC = 6.4D-10  !! 6.40*10^-10 [m^2/s]
  double precision, parameter :: NOPA = 2.50D0

  ! Porosity
  double precision :: PHI(SDZ)
  double precision :: PHII(SDZ+1)

  ! Bioturbation diffusivity
  double precision :: DBI(SDZ+1)

  ! Porewater diffusivity
  double precision :: DEINIT(SDZ+1)
  double precision :: DEIDIC(SDZ+1)
  double precision :: DEIALK(SDZ+1)
  double precision :: DEIDO2(SDZ+1)
  double precision :: DEISIL(SDZ+1)
  double precision :: DEI13DIC(SDZ+1)
  double precision :: DEI14DIC(SDZ+1)

contains

  subroutine set_digda_parameter()
  end subroutine set_digda_parameter

end module module_digda_parameter
