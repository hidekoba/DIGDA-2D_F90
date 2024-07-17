!######################################################################
module module_zcphys
  implicit none

  double precision :: er, grav
  parameter(er = 6.37d6, grav = 9.8d0)
! double precision :: cp, rair
! parameter(cp = 1004.6, rair = 287.04)
! double precision :: el, emelt
! parameter(el = 2.5d6, emelt = 3.4d5)
! double precision :: cpvap, rvap
! parameter(cpvap = 1810., rvap = 461.)
! double precision :: dwatr
! parameter(dwatr = 1.d3)
! double precision :: es0
! parameter(es0 = 611.)
! double precision :: stb
! parameter(stb = 5.67d-8)
! double precision :: fkarm
! parameter(fkarm = 0.4)
! double precision :: tmelt
! parameter(tmelt = 273.15d0)
! double precision :: cpwatr, cpice
! parameter(cpwatr = 4200., cpice = 2000.)
! double precision :: tfrzs, tqice
! parameter(tfrzs = 271.35, tqice = 273.15)
! double precision :: epsv, epsvt
! parameter(epsv = rair/rvap, epsvt = 1.d0/epsv-1.d0)

! double precision ::    gravit
! parameter(gravit = 9.8d+2)
! double precision ::    rhoo, rhoi, rhos, rhow
! parameter(rhoo = 1.d0, rhoi = 0.9d0, rhos = 0.33d0, rhow = 1.d0)
! double precision ::    ckarm
! parameter(ckarm = 0.4d0)
  double precision :: kelvin
  parameter(kelvin = 273.15d0)
! double precision :: hfus
! parameter(hfus = emelt * 1.d4)
! double precision :: cpo, cpi
! parameter(cpo = 3.990d7, cpi= 2.093d7)
! double precision :: cdi, cds
! parameter(cdi = 2.04d5, cds = 0.31d5)
! double precision :: dtds, dtdz
! parameter(dtds = -0.0543d0, dtdz = -7.59d-6)

  double precision :: pi
  parameter(pi = 3.1415926535d0)

end module module_zcphys
