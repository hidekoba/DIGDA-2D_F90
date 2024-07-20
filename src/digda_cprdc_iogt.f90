!######################################################################

! --- information -----------------------------------------------------
!
!  Three-dimensional sediment model
!
!  HISTORY
!     '14.08.01  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_cprdc
      use module_zcfile
      use module_digda_grid
      implicit none

! ==================================================
!   Time
! ==================================================

      INTEGER   tstr, tend, tendyr
      PARAMETER (tstr   =      1)
      PARAMETER (tendyr =      1)  !!      10 [yr]
!     PARAMETER (tendyr =     10)  !!      10 [yr]
!     PARAMETER (tendyr =     20)  !!      20 [yr]
!     PARAMETER (tendyr =     30)  !!      30 [yr]
!     PARAMETER (tendyr =     40)  !!      40 [yr]
!     PARAMETER (tendyr =     50)  !!      50 [yr]
!     PARAMETER (tendyr =    100)  !!     100 [yr]
!     PARAMETER (tendyr =    200)  !!     200 [yr]
!     PARAMETER (tendyr =    300)  !!     300 [yr]
!     PARAMETER (tendyr =    400)  !!     400 [yr]
!     PARAMETER (tendyr =    500)  !!     500 [yr]
!     PARAMETER (tendyr =   1000)  !!   1,000 [yr]
!     PARAMETER (tendyr =   2000)  !!   2,000 [yr]
!     PARAMETER (tendyr =   3000)  !!   3,000 [yr]
!     PARAMETER (tendyr =   4000)  !!   4,000 [yr]
!     PARAMETER (tendyr =   5000)  !!   5,000 [yr]
!     PARAMETER (tendyr =  10000)  !!  10,000 [yr]
!     PARAMETER (tendyr =  20000)  !!  20,000 [yr]
!     PARAMETER (tendyr =  30000)  !!  30,000 [yr]
!     PARAMETER (tendyr =  40000)  !!  40,000 [yr]
!     PARAMETER (tendyr =  50000)  !!  50,000 [yr]
!     PARAMETER (tendyr = 100000)  !! 100,000 [yr]
!     PARAMETER (tendyr = 200000)  !! 200,000 [yr]
!     PARAMETER (tendyr = 300000)  !! 300,000 [yr]
!     PARAMETER (tendyr = 400000)  !! 400,000 [yr]
!     PARAMETER (tendyr = 500000)  !! 500,000 [yr]
      PARAMETER (tend    = tendyr / tsyr )
!     PARAMETER (tend    =      1 )

      INTEGER   tout, toutyr
      PARAMETER (toutyr  =     120000)
      PARAMETER (tout    = toutyr / tsyr )

! ==================================================
!   Sediments
! ==================================================

      REAL*8 SDTX(nxydim, sdzdim, sedtdm)
      SAVE   SDTX
      REAL*8 SDTXBUF(nxdim, nydim, sdzdim, sedtdm)

! --------------------------------------------------
!   Particles
!      Clay, POM, CaCO3, Opal, agetracer
! --------------------------------------------------

!-- Concentration :
!      clay [g/m^3], pom, caco3, opal [mMol/m^3], agetracer [yr]
      REAL*8 scla(nxy, sdz),   sncla(nxy, sdz)
      REAL*8 spom(nxy, sdz),   snpom(nxy, sdz)
      REAL*8 scar(nxy, sdz),   sncar(nxy, sdz)
      REAL*8 sopa(nxy, sdz),   snopa(nxy, sdz)
      REAL*8 sage(nxy, sdz),   snage(nxy, sdz)
      REAL*8 s13pom(nxy, sdz), sn13pom(nxy, sdz)
      REAL*8 s14pom(nxy, sdz), sn14pom(nxy, sdz)
      REAL*8 s13car(nxy, sdz), sn13car(nxy, sdz)
      REAL*8 s14car(nxy, sdz), sn14car(nxy, sdz)
      REAL*8 sclag(nxydim, sdzdim)
      REAL*8 spomg(nxydim, sdzdim)
      REAL*8 scarg(nxydim, sdzdim)
      REAL*8 sopag(nxydim, sdzdim)
      REAL*8 sageg(nxydim, sdzdim)
      REAL*8 s13pomg(nxydim, sdzdim)
      REAL*8 s14pomg(nxydim, sdzdim)
      REAL*8 s13carg(nxydim, sdzdim)
      REAL*8 s14carg(nxydim, sdzdim)

!-- Weight percent [%] :
      REAL*8 pcla(nxy, sdz)
      REAL*8 ppom(nxy, sdz)
      REAL*8 pcar(nxy, sdz)
      REAL*8 popa(nxy, sdz)
      REAL*8 wt(nxy, sdz)
      REAL*8 pclag(nxydim, sdzdim)
      REAL*8 ppomg(nxydim, sdzdim)
      REAL*8 pcarg(nxydim, sdzdim)
      REAL*8 popag(nxydim, sdzdim)
      REAL*8 wtg(nxydim, sdzdim)

!-- Dissolution rate [1/s], Dissolution term [mMol/m^3/s] :
      REAL*8 rcla(nxy, sdz),   trcla(nxy, sdz)
      REAL*8 rpom(nxy, sdz),   trpom(nxy, sdz)
      REAL*8 rcar(nxy, sdz),   trcar(nxy, sdz)
      REAL*8 ropa(nxy, sdz),   tropa(nxy, sdz)
      REAL*8 r13pom(nxy, sdz), tr13pom(nxy, sdz)
      REAL*8 r14pom(nxy, sdz), tr14pom(nxy, sdz)
      REAL*8 r13car(nxy, sdz), tr13car(nxy, sdz)
      REAL*8 r14car(nxy, sdz), tr14car(nxy, sdz)

! --------------------------------------------------
!   Dissolved material
!      Phosphate, Dissolved Inorganic Carbon,
!      Alkalinity, Dissolved Oxygen, Silicate
! --------------------------------------------------

!-- Concentration [mMol/m^3] :
      REAL*8 cnit(nxy, sdz),   cnnit(nxy, sdz)
      REAL*8 cdic(nxy, sdz),   cndic(nxy, sdz)
      REAL*8 calk(nxy, sdz),   cnalk(nxy, sdz)
      REAL*8 cdo2(nxy, sdz),   cndo2(nxy, sdz)
      REAL*8 csil(nxy, sdz),   cnsil(nxy, sdz)
      REAL*8 c13dic(nxy, sdz), cn13dic(nxy, sdz)
      REAL*8 c14dic(nxy, sdz), cn14dic(nxy, sdz)
      REAL*8 cdc13(nxy, sdz),  cdc14(nxy, sdz)
      REAL*8 cco32(nxy, sdz)
      REAL*8 satco32(nxy)
      REAL*8 satsil(nxy,sdz)
      REAL*8 cnitg(nxydim, sdzdim)
      REAL*8 cdicg(nxydim, sdzdim)
      REAL*8 calkg(nxydim, sdzdim)
      REAL*8 cdo2g(nxydim, sdzdim)
      REAL*8 csilg(nxydim, sdzdim)
      REAL*8 c13dicg(nxydim, sdzdim)
      REAL*8 c14dicg(nxydim, sdzdim)
      REAL*8 cdc13g(nxydim, sdzdim)
      REAL*8 cdc14g(nxydim, sdzdim)
      REAL*8 cco32g(nxydim, sdzdim)
      REAL*8 satco32g(nxydim)
      REAL*8 satsilg(nxydim,sdzdim)

!-- Concentration (bottom water) :
!      temp [C], salinity [psu], the others [mMol/m^3]
      REAL*8 bwtem(nxy)
      REAL*8 bwsal(nxy)
      REAL*8 bwnit(nxy)
      REAL*8 bwdic(nxy)
      REAL*8 bwalk(nxy)
      REAL*8 bwdo2(nxy)
      REAL*8 bwsil(nxy)
      REAL*8 bw13dic(nxy)
      REAL*8 bw14dic(nxy)
      REAL*8 bwtemg(nxydim)
      REAL*8 bwsalg(nxydim)
      REAL*8 bwnitg(nxydim)
      REAL*8 bwdicg(nxydim)
      REAL*8 bwalkg(nxydim)
      REAL*8 bwdo2g(nxydim)
      REAL*8 bwsilg(nxydim)
      REAL*8 bw13dicg(nxydim)
      REAL*8 bw14dicg(nxydim)

!-- Concentration (total) :
!      total phosphorous [mMol/m^3], total silicon [mMol/m^3]
      REAL*8 PT(nxy)
      REAL*8 SIT(nxy)

!-- Dissolution term [mMol/m^3/s] :
      REAL*8 trnit(nxy, sdz)
      REAL*8 trdic(nxy, sdz)
      REAL*8 tralk(nxy, sdz)
      REAL*8 trdo2(nxy, sdz)
      REAL*8 trsil(nxy, sdz)
      REAL*8 tr13dic(nxy, sdz)
      REAL*8 tr14dic(nxy, sdz)

!-- Sea-bottom Flux :
!      clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
      REAL*8 ficla(nxy),   focla(nxy)
      REAL*8 fipom(nxy),   fopom(nxy)
      REAL*8 ficar(nxy),   focar(nxy)
      REAL*8 fiopa(nxy),   foopa(nxy)
      REAL*8 fiage(nxy)
      REAL*8 fi13pom(nxy), fo13pom(nxy)
      REAL*8 fi14pom(nxy), fo14pom(nxy)
      REAL*8 fi13car(nxy), fo13car(nxy)
      REAL*8 fi14car(nxy), fo14car(nxy)
      REAL*8 ficlag(nxydim), foclag(nxydim)
      REAL*8 fipomg(nxydim), fopomg(nxydim)
      REAL*8 ficarg(nxydim), focarg(nxydim)
      REAL*8 fiopag(nxydim), foopag(nxydim)
      REAL*8 fiageg(nxydim), foageg(nxydim)
      REAL*8 fi13pomg(nxydim), fo13pomg(nxydim)
      REAL*8 fi14pomg(nxydim), fo14pomg(nxydim)
      REAL*8 fi13carg(nxydim), fo13carg(nxydim)
      REAL*8 fi14carg(nxydim), fo14carg(nxydim)

!-- Sea-bottom Flux [mMol/m^3*cm/s] :
      REAL*8 finit(nxy),   fonit(nxy)
      REAL*8 fidic(nxy),   fodic(nxy)
      REAL*8 fialk(nxy),   foalk(nxy)
      REAL*8 fido2(nxy),   fodo2(nxy)
      REAL*8 fisil(nxy),   fosil(nxy)
      REAL*8 fi13dic(nxy), fo13dic(nxy)
      REAL*8 fi14dic(nxy), fo14dic(nxy)

!-- Sediment-bottom Flux :
!      clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
      REAL*8 bfcla(nxy)
      REAL*8 bfpom(nxy)
      REAL*8 bfcar(nxy)
      REAL*8 bfopa(nxy)
      REAL*8 bf13pom(nxy)
      REAL*8 bf14pom(nxy)
      REAL*8 bf13car(nxy)
      REAL*8 bf14car(nxy)
      REAL*8 bfclag(nxydim)
      REAL*8 bfpomg(nxydim)
      REAL*8 bfcarg(nxydim)
      REAL*8 bfopag(nxydim)
      REAL*8 bf13pomg(nxydim)
      REAL*8 bf14pomg(nxydim)
      REAL*8 bf13carg(nxydim)
      REAL*8 bf14carg(nxydim)

!-- Density, Molar mass :
!      Density [g/cm^3], Molar mass [g/Mol],
!      Molar mass / Density [cm^3/Mol]
      REAL*8 rhocla, rhopom, rhocar, rhoopa
      REAL*8 molcla, molpom, molcar, molopa
      REAL*8 mcla,   mpom,   mcar,   mopa
      SAVE   rhocla, rhopom, rhocar, rhoopa
      SAVE   molcla, molpom, molcar, molopa
      SAVE   mcla,   mpom,   mcar,   mopa
      DATA   rhocla, rhopom, rhocar, rhoopa &
     &    / 2.60d0, 2.60d0, 2.60d0, 2.60d0 /
          !! 2.6, 2.6, 2.6, 2.6 [g/cm^3]
      DATA   molcla, molpom, molcar, molopa &
     &    / 430.0d0, 12.0d0, 100.0d0,  60.0d0 /
          !! 430, 12, 100, 60 [g/Mol]
      REAL*8 rho13pom, rho14pom, rho13car, rho14car
      REAL*8 mol13pom, mol14pom, mol13car, mol14car
      REAL*8 m13pom,   m14pom,   m13car,   m14car
      SAVE   rho13pom, rho14pom, rho13car, rho14car
      SAVE   mol13pom, mol14pom, mol13car, mol14car
      SAVE   m13pom,   m14pom,   m13car,   m14car
      DATA   rho13pom, rho14pom, rho13car, rho14car &
     &    / 2.60d0, 2.60d0, 2.60d0, 2.60d0 /
          !! 2.6, 2.6, 2.6, 2.6 [g/cm^3]
      DATA   mol13pom, mol14pom, mol13car, mol14car &
     &    /  13.0d0, 14.0d0, 101.0d0, 102.0d0 /
          !! 13, 14, 101, 102 [g/Mol]

!-- Redfield ratio :
      REAL*8 FP, PP, NP, CP, OP
          !! Fe:P, P:P, N:P, C:P, O:P
      SAVE   FP, PP, NP, CP, OP
      DATA   FP, PP, NP, CP, OP &
     &    / 4.70d-4, 1.0d0, 16.0d0, 117.0d0, 170.0d0 /
          !! 0.00047, 1, 16, 117, 170

!-- Advection, Diffusion :
!      Sedimentation rate
      REAL*8 w(nxy, sdz+1)
      REAL*8 WG(nxydim, sdzdim+1)
!      Porosity
      REAL*8 phi0, z0
      SAVE   phi0, z0
      DATA   phi0, z0 &
     &    / 0.9d0, 1.823d0 /
          !! 0.9, 1.823 [m]
      REAL*8 phi(sdz), phii(sdz+1)
      SAVE   phi, phii
!      Bioturbation diffusivity
      REAL*8 db0, zb
      SAVE   db0, zb
      DATA   db0, zb &
     &    / 3.0d-5, 1.0d-1 /
          !! 3.0*10-5 [m^2/yr], 0.10 [m]
      REAL*8 dbi(sdz+1)
      SAVE   dbi
!      Pw Effective molecular diffusivity
      REAL*8 dsw
      REAL*8 dswnit, dswdic, dswalk, dswdo2, dswsil
      REAL*8 dsw13dic, dsw14dic
      REAL*8 nopa
      SAVE   dsw
      SAVE   dswnit, dswdic, dswalk, dswdo2, dswsil
      SAVE   dsw13dic, dsw14dic
      SAVE   nopa
      DATA   dsw &
     &    / 4.59d-10 /
          !! 4.59*10-10 [m^2/s]
      DATA   dswnit, dswdic, dswalk, dswdo2, dswsil &
     &    / 5.0d-10, 6.4d-10, 6.4d-10, 1.2d-9, 5.0d-10 /
          !! 5.0*10-10, 6.4*10-10, 6.4*10-10, 1.2*10-9, 5.0*10-10  [m^2/s]
      DATA   dsw13dic, dsw14dic &
     &    / 6.4d-10, 6.4d-10 /
          !! 6.4*10-10, 6.4*10-10 [m^2/s]
      DATA   nopa &
     &    / 2.50d0 /
      REAL*8 deinit(sdz+1)
      REAL*8 deidic(sdz+1)
      REAL*8 deialk(sdz+1)
      REAL*8 deido2(sdz+1)
      REAL*8 deisil(sdz+1)
      REAL*8 dei13dic(sdz+1)
      REAL*8 dei14dic(sdz+1)
      SAVE   deinit, deidic, deialk, deido2, deisil
      SAVE   dei13dic, dei14dic

!-- 13C, 14C :
      REAL*8 pdb              !! standard 13C/12C
      REAL*8 dstd             !! standard 14C/12C
      REAL*8 eoph, eosh, eps  !! fractionaton factor
      REAL*8 feoph, feosh, r13, r14
      REAL*8 thalf, adcay     !! half decay time [yr], decay constant
      REAL*8 dltc13, dltc14
      SAVE   pdb
      SAVE   dstd
      SAVE   eoph, eosh
      SAVE   feoph, feosh
      SAVE   thalf, adcay
      DATA   pdb &
     &    / 0.0112372d0 /
      DATA   dstd &
     &    / 1.1892d-12 /
      DATA   eoph, eosh &
     &    / -22.0, 0.0 /
      DATA   eps &
     &    / 1.0d-10 /
      DATA   thalf &
     &    / 5370.d0 /

!-- MASK :
      REAL*8 ZW(nz)
      SAVE ZW
      REAL*8 bdep(nxy), bdepg(nxydim) !! ocean depth [m]
      SAVE bdep, bdepg

!-- Parallel computing :
      INTEGER, ALLOCATABLE :: nrcvct2d(:),  ndispt2d(:)
      INTEGER, ALLOCATABLE :: nrcvct3d(:),  ndispt3d(:)
      INTEGER, ALLOCATABLE :: nrcvcw3d(:),  ndispw3d(:)
      REAL*8,  ALLOCATABLE :: DUMMY01(:), DUMMY02(:), DUMMY03(:)
      REAL*8,  ALLOCATABLE :: DUMMY04(:), DUMMY05(:), DUMMY06(:)
      REAL*8,  ALLOCATABLE :: DUMMY07(:), DUMMY08(:), DUMMY09(:)
      REAL*8,  ALLOCATABLE :: DUMMY11(:), DUMMY12(:), DUMMY13(:)
      REAL*8,  ALLOCATABLE :: DUMMY14(:), DUMMY15(:), DUMMY16(:)
      REAL*8,  ALLOCATABLE :: DUMMY17(:), DUMMY18(:), DUMMY19(:)
      REAL*8,  ALLOCATABLE :: DUMMY20(:)
      REAL*8,  ALLOCATABLE :: DUMMY21(:), DUMMY22(:), DUMMY23(:)
      REAL*8,  ALLOCATABLE :: DUMMY31(:), DUMMY32(:), DUMMY33(:)
      REAL*8,  ALLOCATABLE :: DUMMY34(:), DUMMY35(:), DUMMY36(:)
      REAL*8,  ALLOCATABLE :: DUMMY37(:), DUMMY38(:)
      REAL*8,  ALLOCATABLE :: DUMMY41(:), DUMMY42(:), DUMMY43(:)
      REAL*8,  ALLOCATABLE :: DUMMY44(:), DUMMY45(:), DUMMY46(:)
      REAL*8,  ALLOCATABLE :: DUMMY47(:), DUMMY48(:)
      REAL*8,  ALLOCATABLE :: DUMMY51(:), DUMMY52(:), DUMMY53(:)
      REAL*8,  ALLOCATABLE :: DUMMY54(:), DUMMY55(:), DUMMY56(:)
      REAL*8,  ALLOCATABLE :: DUMMY57(:), DUMMY58(:)

!-- Option :
      LOGICAL SDFIrst
      SAVE SDFIrst
      DATA SDFIrst / .true. /
      LOGICAL SDofirst
      SAVE SDofirst
      DATA SDofirst / .true.  /

      LOGICAL DEIARC
      SAVE DEIARC
      namelist /NMSDDEAR/ DEIARC
      DATA DEIARC / .true. /  !! Archer's DEI formulation ?

      LOGICAL  CALc13, CALc14
      SAVE     CALc13, CALc14
      namelist /NMSDISO/ CALc13, CALc14
      DATA CALc13, CALc14 / .true., .true. /  !! carbon isotopes ?

!-- etc :
      INTEGER I, J, ij, k, L, IC
      CHARACTER SEDTNM(sedtdm)*16
      SAVE SEDTNM
      INTEGER TT

      LOGICAL  wotstp
      SAVE     wotstp
      DATA     wotstp / .true. /  !! writeout time step ?
      INTEGER     nfsots
      PARAMETER ( nfsots = 999 )

      CHARACTER CFILE*128, CTITLE*128, CUNIT*128

end module globals_digda_cprdc

! =====================================================================

module module_digda_cprdc
      use globals_digda_cprdc
      use module_digda_clay
      use module_digda_pom
      use module_digda_carbonate
      use module_digda_silicate
      use module_digda_step_solid
      use module_digda_step_liquid
      use module_digda_chemeq
      implicit none
contains

      subroutine digda_cprdc()
    ! subroutine digda_cprdc( &
    !&     bdepg, &
    !&     bfrcg, &
    !&     sdtx )

!-- SEDIMENT TRACERS : setup --

if (use_mpi) then
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(mpi_comm_world, nprocs, IERR)
      call MPI_COMM_RANK(mpi_comm_world, myrank, IERR)

      if ( nprocs /= ncpu ) then
         if ( myrank == 0 ) then
            PRINT *, '=== ERROR ==='
            call MPI_FINALIZE(IERR)
            STOP
         end if
      end if

      do irank = 0, nprocs-1
         call set_parallel_range( &
     &                    1, nxydim, nprocs, irank, &
     &                    jjstr(irank), jjend(irank) &
     &                  )
      end do
end if

      OPEN(nfparm, file="PARAMET")
      OPEN(nfstdo, file="STDOUT")

      write(nfstdo, *) '-----------------------------'
      write(nfstdo, *) '-- Setup of Sediment Model --'
      write(nfstdo, *) '-----------------------------'
      rewind(nfparm)
      read(nfparm, NMSDDEAR, end=99)
99 continue
      write(nfstdo, NMSDDEAR)
      rewind(nfparm)
      read(nfparm, NMSDISO, end=98)
98 continue
      write(nfstdo, NMSDISO)

!     ### for Sediment Model ###
      SEDTNM(1)  = 'scla'
      SEDTNM(2)  = 'spom'
      SEDTNM(3)  = 'scar'
      SEDTNM(4)  = 'sopa'
      SEDTNM(5)  = 'sage'
      SEDTNM(6)  = 'cnit'
      SEDTNM(7)  = 'cdic'
      SEDTNM(8)  = 'calk'
      SEDTNM(9)  = 'csil'
      SEDTNM(10) = 'cdo2'
      SEDTNM(11) = 's13pom'
      SEDTNM(12) = 's13car'
      SEDTNM(13) = 'c13dic'
      SEDTNM(14) = 's14pom'
      SEDTNM(15) = 's14car'
      SEDTNM(16) = 'c14dic'

!CCCCCcall SDFILESET
      call set_digda_clay
      call set_digda_pom
      call set_digda_carbonate
      call set_digda_silicate
      call set_digda_chemeq
      call set_digda_step_solid
      call set_digda_step_liquid

!-- SEDIMENT TRACERS : setup (end) --

!     -- Ocean depth --
      CFILE="bound/bdep"
      call read_gtool( CFILE, nxydim,  bdepg )
if (use_mpi) then
      do ij = 1, nxy
         bdep(ij) = max( bdepg(myrank*nxy+ij), 0.0d0 )
      end do
else
      do ij = 1, nxy
         bdep(ij) = max( bdepg(ij), 0.0d0 )
      end do
end if

!     -- Convert unit --

      z0  = z0  * 1.0d+2                       !! [m] --> [cm]

      db0 = db0 * 1.0d+4 / 3.1536d+7           !! [m^2/yr] --> [cm^2/s]
      zb  = zb  * 1.0d+2                       !! [m] --> [cm]

      do k = skst, sken
         lz(k)   = lz(k)   * 1.0d+2            !! [m] --> [cm]
         dlz(k)  = dlz(k)  * 1.0d+2            !! [m] --> [cm]
      end do
      do k = skst, sken+1
         lzi(k)  = lzi(k)  * 1.0d+2            !! [m] --> [cm]
         dlzi(k) = dlzi(k) * 1.0d+2            !! [m] --> [cm]
      end do

      dswnit   = dswnit   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      dswdic   = dswdic   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      dswalk   = dswalk   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      dswdo2   = dswdo2   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      dswsil   = dswsil   * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      if ( CALc13 ) then
      dsw13dic = dsw13dic * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      end if
      if ( CALc14 ) then
      dsw14dic = dsw14dic * 1.0d+4             !! [m^2/s] --> [cm^2/s]
      end if
      dsw    = dsw        * 1.0d+4             !! [m^2/s] --> [cm^2/s]

      mcla   = 1.0d-6 /  rhocla                !! [m^3/g]
      mpom   = 1.0d-9 / ( rhopom / molpom )    !! [m^3/mMol]
      mcar   = 1.0d-9 / ( rhocar / molcar )    !! [m^3/mMol]
      mopa   = 1.0d-9 / ( rhoopa / molopa )    !! [m^3/mMol]
      if ( CALc13 ) then
      m13pom = 1.0d-9 / ( rhopom / mol13pom )  !! [m^3/mMol]
      m13car = 1.0d-9 / ( rhocar / mol13car )  !! [m^3/mMol]
      end if
      if ( CALc14 ) then
      m14pom = 1.0d-9 / ( rhopom / mol14pom )  !! [m^3/mMol]
      m14car = 1.0d-9 / ( rhocar / mol14car )  !! [m^3/mMol]
      thalf  = thalf * 3.1536d+7               !! [yr] --> [s]
      adcay  = log( 2.0d0 ) / thalf / 365.0d0 / 8.64d+4
      end if
      if ( CALc13 .or. CALc14 ) then
      feoph  = 1.0d0 + eoph * 1.0d-3
      feosh  = 1.0d0 + eosh * 1.0d-3
      end if

      do ij = 1, nxy
         bdep(ij) = bdep(ij) * 1.0d+2          !! [m] --> [cm]
      end do

!     -- Porosity, Diffusivity --

      do k = skst, sken
         phi(k) = phi0 * exp( - lz(k)  / z0 )
      end do

      do k = skst, sken+1
         phii(k) = phi0 * exp( - lzi(k) / z0 )
         dbi(k)  = ( db0 / ( 1.0d0 + exp( 2.0d0 * ( lzi(k) - zb ) ) ) )
      end do

      do k = skst, sken+1
         if (DEIARC) then
            deinit(k)   = dswnit   * ( phii(k) ** 3.0d0 )
            deidic(k)   = dswdic   * ( phii(k) ** 3.0d0 )
            deialk(k)   = dswalk   * ( phii(k) ** 3.0d0 )
            deido2(k)   = dswdo2   * ( phii(k) ** 3.0d0 )
            deisil(k)   = dswsil   * ( phii(k) ** 3.0d0 )
            if ( CALc13 ) then
            dei13dic(k) = dsw13dic * ( phii(k) ** 3.0d0 )
            end if
            if ( CALc14 ) then
            dei14dic(k) = dsw14dic * ( phii(k) ** 3.0d0 )
            end if
         else
            deinit(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            deidic(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            deialk(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            deido2(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            deisil(k)   = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            if ( CALc13 ) then
            dei13dic(k) = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            end if
            if ( CALc14 ) then
            dei14dic(k) = ( phii(k) ** ( nopa - 1.0d0 ) ) * dsw
            end if
         end if
      end do

!     -- Set boundary conditions --

      CFILE="bound/ficla"
      call read_gtool( CFILE, nxydim,  ficlag )
      CFILE="bound/fipom"
      call read_gtool( CFILE, nxydim,  fipomg )
      CFILE="bound/ficar"
      call read_gtool( CFILE, nxydim,  ficarg )
      CFILE="bound/fiopa"
      call read_gtool( CFILE, nxydim,  fiopag )
if (use_mpi) then
      do ij = 1, nxy
         ficla(ij) = max( ficlag(myrank*nxy+ij), 0.0d0 )
         fipom(ij) = max( fipomg(myrank*nxy+ij), 0.0d0 )
         ficar(ij) = max( ficarg(myrank*nxy+ij), 0.0d0 )
         fiopa(ij) = max( fiopag(myrank*nxy+ij), 0.0d0 )
         fiage(ij) = 0.0d0
      end do
else
      do ij = 1, nxy
         ficla(ij) = max( ficlag(ij), 0.0d0 )
         fipom(ij) = max( fipomg(ij), 0.0d0 )
         ficar(ij) = max( ficarg(ij), 0.0d0 )
         fiopa(ij) = max( fiopag(ij), 0.0d0 )
         fiage(ij) = 0.0d0
      end do
end if
      CFILE="bound/btem"
      call read_gtool( CFILE, nxydim,  bwtemg )
      CFILE="bound/bsal"
      call read_gtool( CFILE, nxydim,  bwsalg )
      CFILE="bound/bnit"
      call read_gtool( CFILE, nxydim,  bwnitg )
      CFILE="bound/bdic"
      call read_gtool( CFILE, nxydim,  bwdicg )
      CFILE="bound/balk"
      call read_gtool( CFILE, nxydim,  bwalkg )
      CFILE="bound/bdo2"
      call read_gtool( CFILE, nxydim,  bwdo2g )
      CFILE="bound/bsil"
      call read_gtool( CFILE, nxydim,  bwsilg )
if (use_mpi) then
      do ij = 1, nxy
         bwtem(ij)   = bwtemg(myrank*nxy+ij)
         bwsal(ij)   = bwsalg(myrank*nxy+ij)
         bwnit(ij)   = max( bwnitg(myrank*nxy+ij), 0.0d0 )
         bwdic(ij)   = max( bwdicg(myrank*nxy+ij), 0.0d0 )
         bwalk(ij)   = max( bwalkg(myrank*nxy+ij), 0.0d0 )
         bwdo2(ij)   = max( bwdo2g(myrank*nxy+ij), 0.0d0 )
         bwsil(ij)   = max( bwsilg(myrank*nxy+ij), 0.0d0 )
      end do
else
      do ij = 1, nxy
         bwtem(ij)   = bwtemg(ij)
         bwsal(ij)   = bwsalg(ij)
         bwnit(ij)   = max( bwnitg(ij), 0.0d0 )
         bwdic(ij)   = max( bwdicg(ij), 0.0d0 )
         bwalk(ij)   = max( bwalkg(ij), 0.0d0 )
         bwdo2(ij)   = max( bwdo2g(ij), 0.0d0 )
         bwsil(ij)   = max( bwsilg(ij), 0.0d0 )
      end do
end if

      if ( CALc13 ) then
      do ij = ijtstr, ijtend
         fi13pom(ij) = feoph * pdb  * fipom(ij)
         fi13car(ij) = feosh * pdb  * ficar(ij)
      end do
      CFILE="bound/bc13"
      call read_gtool( CFILE, nxydim,  bw13dicg )
if (use_mpi) then
      do ij = 1, nxy
         bw13dic(ij) = max( bw13dicg(myrank*nxy+ij), 0.0d0 )
      end do
else
      do ij = 1, nxy
         bw13dic(ij) = max( bw13dicg(ij), 0.0d0 )
      end do
end if
      end if

      if ( CALc14 ) then
      do ij = ijtstr, ijtend
         fi14pom(ij) = feoph * dstd * fipom(ij)
         fi14car(ij) = feosh * dstd * ficar(ij)
      end do
      CFILE="bound/bc14"
      call read_gtool( CFILE, nxydim,  bw14dicg )
if (use_mpi) then
      do ij = 1, nxy
         bw14dic(ij) = max( bw14dicg(myrank*nxy+ij), 0.0d0 )
      end do
else
      do ij = 1, nxy
         bw14dic(ij) = max( bw14dicg(ij), 0.0d0 )
      end do
end if
      end if

      do ij = 1, nxy
         SIT(ij) = 0.51d0
         PT(ij)  = 7.70d0
      end do

!     do l = 1, sedtdm
!        do k = 1, sdzdim
!           DO J = 1, nydim
!              DO I = 1, nxdim
!                 ij = i + (j - 1) * nxdim
!                 SDTXBUF(i, j, k, l) = SDTX(ij, k, l)
!              end do
!           end do
!        end do
!     end do
!     do l = 1, sedtdm
!        do k = skst, sken
!           DO J = jstr, jend
!              DO I = 1, istr-1
!                 SDTXBUF(i, j, k, l) = SDTXBUF(iend-istr+1+I, j, k, l)
!              end do
!              DO I = iend+1, nxdim
!                 SDTXBUF(i, j, k, l) = SDTXBUF(istr-1+I-iend, j, k, l)
!              end do
!           end do
!        end do
!     end do
!     do l = 1, sedtdm
!        do k = 1, sdzdim
!           DO J = 1, nydim
!              DO I = 1, nxdim
!                 ij = i + (j - 1) * nxdim
!                 SDTX(ij, k, l) = SDTXBUF(i, j, k, l)
!              end do
!           end do
!        end do
!     end do
!     do l = 1, sedtdm
!        do k = 1, skst-1
!           do ij = 1, nxy
!              SDTX(ij, k, l) = SDTX(ij, skst, l)
!           end do
!        end do
!     end do

!     -- Time Integration --

      do tt = tstr, tend

         YY = int( tsyr * TT )

      if (SDFIrst) then
         do k = skst, sken
            do ij = ijtstr, ijtend
               scla(ij, k) = rhocla * 1.0d+6
               spom(ij, k) = 0.0d0
               scar(ij, k) = 0.0d0
               sopa(ij, k) = 0.0d0
               sage(ij, k) = 0.0d0
               cnit(ij, k) = bwnit(ij)
               cdic(ij, k) = bwdic(ij)
               calk(ij, k) = bwalk(ij)
               cdo2(ij, k) = bwdo2(ij)
               csil(ij, k) = bwsil(ij)
               if ( CALc13 ) then
               s13pom(ij, k) = 0.0d0
               s13car(ij, k) = 0.0d0
               c13dic(ij, k) = bw13dic(ij)
               end if
               if ( CALc14 ) then
               s14pom(ij, k) = 0.0d0
               s14car(ij, k) = 0.0d0
               c14dic(ij, k) = bw14dic(ij)
               end if
               pcla(ij, k) = 1.00+2
               ppom(ij, k) = 0.0d0
               pcar(ij, k) = 0.0d0
               popa(ij, k) = 0.0d0
            end do
         end do
         SDFIrst = .false.
      else
if (use_mpi) then
         do k = skst, sken
            do ij = 1, nxy
               scla(ij, k) = SDTX(myrank*nxy+ij, k,  1)
               spom(ij, k) = SDTX(myrank*nxy+ij, k,  2)
               scar(ij, k) = SDTX(myrank*nxy+ij, k,  3)
               sopa(ij, k) = SDTX(myrank*nxy+ij, k,  4)
               sage(ij, k) = SDTX(myrank*nxy+ij, k,  5)
               cnit(ij, k) = SDTX(myrank*nxy+ij, k,  6)
               cdic(ij, k) = SDTX(myrank*nxy+ij, k,  7)
               calk(ij, k) = SDTX(myrank*nxy+ij, k,  8)
               csil(ij, k) = SDTX(myrank*nxy+ij, k,  9)
               cdo2(ij, k) = SDTX(myrank*nxy+ij, k, 10)
            end do
         end do
         if ( CALc13 ) then
         do k = skst, sken
            do ij = 1, nxy
               s13pom(ij, k) = SDTX(myrank*nxy+ij, k, 11)
               s13car(ij, k) = SDTX(myrank*nxy+ij, k, 12)
               c13dic(ij, k) = SDTX(myrank*nxy+ij, k, 13)
            end do
         end do
         end if
         if ( CALc14 ) then
         do k = skst, sken
            do ij = 1, nxy
               s14pom(ij, k) = SDTX(myrank*nxy+ij, k, 14)
               s14car(ij, k) = SDTX(myrank*nxy+ij, k, 15)
               c14dic(ij, k) = SDTX(myrank*nxy+ij, k, 16)
            end do
         end do
         end if
else
         do k = skst, sken
            do ij = 1, nxy
               scla(ij, k) = SDTX(ij, k,  1)
               spom(ij, k) = SDTX(ij, k,  2)
               scar(ij, k) = SDTX(ij, k,  3)
               sopa(ij, k) = SDTX(ij, k,  4)
               sage(ij, k) = SDTX(ij, k,  5)
               cnit(ij, k) = SDTX(ij, k,  6)
               cdic(ij, k) = SDTX(ij, k,  7)
               calk(ij, k) = SDTX(ij, k,  8)
               csil(ij, k) = SDTX(ij, k,  9)
               cdo2(ij, k) = SDTX(ij, k, 10)
            end do
         end do
         if ( CALc13 ) then
         do k = skst, sken
            do ij = 1, nxy
               s13pom(ij, k) = SDTX(ij, k, 11)
               s13car(ij, k) = SDTX(ij, k, 12)
               c13dic(ij, k) = SDTX(ij, k, 13)
            end do
         end do
         end if
         if ( CALc14 ) then
         do k = skst, sken
            do ij = 1, nxy
               s14pom(ij, k) = SDTX(ij, k, 14)
               s14car(ij, k) = SDTX(ij, k, 15)
               c14dic(ij, k) = SDTX(ij, k, 16)
            end do
         end do
         end if
end if
      end if

!     -- Saturated carbonate ion concentration --

      call calc_digda_chemeq( &
     &    bwtem, bwsal, cdic, calk, &
     &    SIT, PT, &
     &    cco32 )

      call calc_carbonate_saturation( &
     &    bwtem, bwsal, bdep, &
     &    satco32 )

!     -- Dissolution term --

      call calc_clay_dissolution( &
     &    scla, &
     &    rcla )

      call calc_pom_dissolution( &
     &    spom, cdo2, CP, OP, phi, &
     &    rpom )

      call calc_carbonate_dissolution( &
     &    scar, cco32, satco32, phi, &
     &    rcar )

      call calc_silicate_dissolution( &
     &    sopa, csil, pcla, popa, bwtem, &
     &    ropa, satsil )

!     -- Dissolution rate --
!        CLA, POM, CAR, OPA

      do k = skst, sken
         do ij = ijtstr, ijtend
            trcla(ij, k)   = ( 1.0d0 - phi(k) ) * rcla(ij, k) &
     &                     * scla(ij, k)
            trpom(ij, k)   = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * spom(ij, k)
            trcar(ij, k)   = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * scar(ij, k)
            tropa(ij, k)   = ( 1.0d0 - phi(k) ) * ropa(ij, k) &
     &                     * sopa(ij, k)
         end do
      end do

      if ( CALc13 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr13pom(ij, k) = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * s13pom(ij, k)
            tr13car(ij, k) = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * s13car(ij, k)
         end do
      end do
      end if

      if ( CALc14 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr14pom(ij, k) = ( 1.0d0 - phi(k) ) * rpom(ij, k) &
     &                     * s14pom(ij, k)
            tr14car(ij, k) = ( 1.0d0 - phi(k) ) * rcar(ij, k) &
     &                     * s14car(ij, k)
         end do
      end do
      end if

!     -- Sedimentation rate --

      do ij = ijtstr, ijtend
         w(ij, skst) = ( ( ( ficla(ij) * mcla ) &
     &                   + ( fipom(ij) * mpom ) &
     &                   + ( ficar(ij) * mcar ) &
     &                   + ( fiopa(ij) * mopa ) ) &
     &               / ( 1.0d0 - phii(skst) ) )
         do k = skst, sken
            w(ij, k+1) = ( ( 1.0d0 - phii(k) ) / ( 1.0d0 - phii(k+1) ) &
     &                   * w(ij, k) ) &
     &                 - ( dlz(k) * ( ( trpom(ij, k) * mpom ) &
     &                              + ( trcar(ij, k) * mcar ) &
     &                              + ( tropa(ij, k) * mopa ) ) &
     &                   / ( 1.0d0 - phii(k+1) ) )
         end do
      end do

!     -- Particles (t --> t+dt) --
!        Clay, POM, CaCO3, Opal

      call calc_step_solid_clay( &
     &    ficla, trcla, w, phi, phii, dbi, scla, rhocla, &
     &    sncla, bfcla )

      call calc_step_solid( &
     &    fipom, trpom, w, phi, phii, dbi, spom, &
     &    snpom, bfpom )

      call calc_step_solid( &
     &    ficar, trcar, w, phi, phii, dbi, scar, &
     &    sncar, bfcar )

      call calc_step_solid( &
     &    fiopa, tropa, w, phi, phii, dbi, sopa, &
     &    snopa, bfopa )

      call calc_step_solid_age( &
     &    fiage, w, phi, phii, dbi, sage, &
     &    snage )

      if ( CALc13 ) then
      call calc_step_solid( &
     &    fi13pom, tr13pom, w, phi, phii, dbi, s13pom, &
     &    sn13pom, bf13pom )

      call calc_step_solid( &
     &    fi13car, tr13car, w, phi, phii, dbi, s13car, &
     &    sn13car, bf13car )
      end if

      if ( CALc14 ) then
      call calc_step_solid( &
     &    fi14pom, tr14pom, w, phi, phii, dbi, s14pom, &
     &    sn14pom, bf14pom )

      call calc_step_solid( &
     &    fi14car, tr14car, w, phi, phii, dbi, s14car, &
     &    sn14car, bf14car )
      do k = skst, sken
         do ij = ijtstr, ijtend
            sn14pom(ij, k) = sn14pom(ij, k) &
     &                     - adcay * ts * s14pom(ij, k)
            sn14car(ij, k) = sn14car(ij, k) &
     &                     - adcay * ts * s14car(ij, k)
         end do
      end do
      end if

!     -- Dissolution rate --
!        Phosphate, DIC, Alkalinity, Oxygen, Silicate

      do k = skst, sken
         do ij = ijtstr, ijtend
            trnit(ij, k)   = ( ( PP / CP ) * trpom(ij, k) )
            trdic(ij, k)   = ( trpom(ij, k) + trcar(ij, k) )
            tralk(ij, k)   = ( - ( NP / CP ) * trpom(ij, k) &
     &                         + ( 2.0d0 * trcar(ij, k) ) )
            trdo2(ij, k)   = - ( ( OP / CP ) * trpom(ij, k) )
            if ( cdo2(ij, k) <= 0.0d0 ) then
               trdo2(ij, k) = 0.0d0
            end if
            trsil(ij, k)   = tropa(ij, k)
         end do
      end do

      if ( CALc13 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr13dic(ij, k) = ( tr13pom(ij, k) + tr13car(ij, k) )
         end do
      end do
      end if

      if ( CALc14 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            tr14dic(ij, k) = ( tr14pom(ij, k) + tr14car(ij, k) )
         end do
      end do
      end if

!     -- Dissolved material (t --> t+dt) --
!        Phosphate, DIC, Alkalinity, Oxygen, Silicate

      call calc_step_liquid( &
     &    bwnit, trnit, phi, phii, deinit, cnit, &
     &    cnnit )

      call calc_step_liquid( &
     &    bwdic, trdic, phi, phii, deidic, cdic, &
     &    cndic )

      call calc_step_liquid( &
     &    bwalk, tralk, phi, phii, deialk, calk, &
     &    cnalk )

      call calc_step_liquid_o2( &
     &    bwdo2, trdo2, phi, phii, deido2, cdo2, &
     &    cndo2 )

      call calc_step_liquid( &
     &    bwsil, trsil, phi, phii, deisil, csil, &
     &    cnsil )

      if ( CALc13 ) then
      call calc_step_liquid( &
     &    bw13dic, tr13dic, phi, phii, dei13dic, c13dic, &
     &    cn13dic )
      end if

      if ( CALc14 ) then
      call calc_step_liquid( &
     &    bw14dic, tr14dic, phi, phii, dei14dic, c14dic, &
     &    cn14dic )
      end if
      do k = skst, sken
         do ij = ijtstr, ijtend
            cn14dic(ij, k) = cn14dic(ij, k) &
     &                     - adcay * ts * c14dic(ij, k)
         end do
      end do

!     -- Particles weight percent --

      do k = skst, sken
         do ij = ijtstr, ijtend
            wt(ij, k) = sncla(ij, k) &
     &                + ( snpom(ij, k) * molpom * 1.0d-3 ) &
     &                + ( sncar(ij, k) * molcar * 1.0d-3 ) &
     &                + ( snopa(ij, k) * molopa * 1.0d-3 )
            if ( wt(ij, k) == 0.0d0 ) then
               pcla(ij, k) = 0.0d0
               ppom(ij, k) = 0.0d0
               pcar(ij, k) = 0.0d0
               popa(ij, k) = 0.0d0
            else
               pcla(ij, k) = &
     &            ( sncla(ij, k) / wt(ij, k) ) &
     &            * 1.0d+2
               ppom(ij, k) = &
     &            ( snpom(ij, k) * molpom * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
               pcar(ij, k) = &
     &            ( sncar(ij, k) * molcar * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
               popa(ij, k) = &
     &            ( snopa(ij, k) * molopa * 1.0d-3 / wt(ij, k) ) &
     &            * 1.0d+2
            end if
         end do
      end do

!     -- Isotopic fractionation --

      if ( CALc13 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            if ( cndic(ij, k) /= 0.0d0 ) then
               r13 = cn13dic(ij, k) / max( cndic(ij, k), eps )
               cdc13(ij, k) = &
     &            ( 1.0d0 / 0.0112372d0 &
     &              * r13 / ( 1.0d0 - r13 ) - 1.0d0 ) * 1.0d+3
            else
               cdc13(ij, k) = 0.0d0
            end if
         end do
      end do
      end if

      if ( CALc14 ) then
      do k = skst, sken
         do ij = ijtstr, ijtend
            if ( cndic(ij, k) /= 0.0d0 ) then
               r14 = cn14dic(ij, k) / max( cndic(ij, k), eps )
               cdc14(ij, k) = &
     &            ( 1.0d0 / 1.1892d-12 &
     &              * r14 / ( 1.0d0 - r14 ) - 1.0d0 ) * 1.0d+3
            else
               cdc14(ij, k) = 0.0d0
            end if
         end do
      end do
      do k = skst, sken
         do ij = ijtstr, ijtend
            dltc13 = cdc13(ij, k) * 1.0d-3  !! [permil]
            dltc14 = cdc14(ij, k) * 1.0d-3  !! [permil]
            cdc14(ij, k) = &
     &         ( dltc14 - 2.0d0 * ( dltc13 + 0.025d0 ) &
     &                          * ( 1.0d0 + dltc14 ) &
     &         ) * 1.0d+3   !! [permil]
         end do
      end do
      end if

!     -- Sea-bottom flux --

      do ij = ijtstr, ijtend
         fonit(ij) = - ( phii(skst) * deinit(skst) &
     &      * ( bwnit(ij) - cnnit(ij, skst) ) / dlzi(skst) )
         fodic(ij) = - ( phii(skst) * deidic(skst) &
     &      * ( bwdic(ij) - cndic(ij, skst) ) / dlzi(skst) )
         foalk(ij) = - ( phii(skst) * deialk(skst) &
     &      * ( bwalk(ij) - cnalk(ij, skst) ) / dlzi(skst) )
         fodo2(ij) = - ( phii(skst) * deido2(skst) &
     &      * ( bwdo2(ij) - cndo2(ij, skst) ) / dlzi(skst) )
         fosil(ij) = - ( phii(skst) * deisil(skst) &
     &      * ( bwsil(ij) - cnsil(ij, skst) ) / dlzi(skst) )
         finit(ij)  = fipom(ij) * ( PP / CP )
         fidic(ij)  = fipom(ij) + ficar(ij)
         fialk(ij)  = - fipom(ij) * ( NP / CP ) + 2.0d0 * ficar(ij)
         fido2(ij)  = fipom(ij) * ( OP / CP )
         fisil(ij)  = fiopa(ij)
         if ( CALc13 ) then
         fo13dic(ij) = - ( phii(skst) * dei13dic(skst) &
     &      * ( bw13dic(ij) - cn13dic(ij, skst) ) / dlzi(skst) )
         fi13dic(ij)  = fi13pom(ij) + fi13car(ij)
         end if
         if ( CALc14 ) then
         fo14dic(ij) = - ( phii(skst) * dei14dic(skst) &
     &      * ( bw14dic(ij) - cn14dic(ij, skst) ) / dlzi(skst) )
         fi14dic(ij)  = fi14pom(ij) + fi14car(ij)
         end if
      end do

      do ij = ijtstr, ijtend
         fopom(ij) = CP * fonit(ij)
         focar(ij) = 0.5d0 * ( + NP * fonit(ij) + foalk(ij) )
!         fopom(ij) = ( CP / ( NP + 2.0d0 * CP ) ) &
!     &      * ( 2.0d0 * fodic(ij) - foalk(ij) )
!         focar(ij) = ( CP / ( NP + 2.0d0 * CP ) ) &
!     &      * ( ( NP / CP ) * fodic(ij) + foalk(ij) )
         foopa(ij) = fosil(ij)
      end do

!     -- OUTPUT --
if (use_mpi) then
      do k = skst, sken
         do ij = 1, nxy
            SDTX(myrank*nxy+ij, k,  1) = sncla(ij, k)
            SDTX(myrank*nxy+ij, k,  2) = snpom(ij, k)
            SDTX(myrank*nxy+ij, k,  3) = sncar(ij, k)
            SDTX(myrank*nxy+ij, k,  4) = snopa(ij, k)
            SDTX(myrank*nxy+ij, k,  5) = snage(ij, k)
            SDTX(myrank*nxy+ij, k,  6) = cnnit(ij, k)
            SDTX(myrank*nxy+ij, k,  7) = cndic(ij, k)
            SDTX(myrank*nxy+ij, k,  8) = cnalk(ij, k)
            SDTX(myrank*nxy+ij, k,  9) = cnsil(ij, k)
            SDTX(myrank*nxy+ij, k, 10) = cndo2(ij, k)
         end do
      end do
      if ( CALc13 ) then
      do k = skst, sken
         do ij = 1, nxy
            SDTX(myrank*nxy+ij, k, 11) = sn13pom(ij, k)
            SDTX(myrank*nxy+ij, k, 12) = sn13car(ij, k)
            SDTX(myrank*nxy+ij, k, 13) = cn13dic(ij, k)
         end do
      end do
      end if
      if ( CALc14 ) then
      do k = skst, sken
         do ij = 1, nxy
            SDTX(myrank*nxy+ij, k, 14) = sn14pom(ij, k)
            SDTX(myrank*nxy+ij, k, 15) = sn14car(ij, k)
            SDTX(myrank*nxy+ij, k, 16) = cn14dic(ij, k)
         end do
      end do
      end if
else
      do k = skst, sken
         do ij = 1, nxy
            SDTX(ij, k,  1) = sncla(ij, k)
            SDTX(ij, k,  2) = snpom(ij, k)
            SDTX(ij, k,  3) = sncar(ij, k)
            SDTX(ij, k,  4) = snopa(ij, k)
            SDTX(ij, k,  5) = snage(ij, k)
            SDTX(ij, k,  6) = cnnit(ij, k)
            SDTX(ij, k,  7) = cndic(ij, k)
            SDTX(ij, k,  8) = cnalk(ij, k)
            SDTX(ij, k,  9) = cnsil(ij, k)
            SDTX(ij, k, 10) = cndo2(ij, k)
         end do
      end do
      if ( CALc13 ) then
      do k = skst, sken
         do ij = 1, nxy
            SDTX(ij, k, 11) = sn13pom(ij, k)
            SDTX(ij, k, 12) = sn13car(ij, k)
            SDTX(ij, k, 13) = cn13dic(ij, k)
         end do
      end do
      end if
      if ( CALc14 ) then
      do k = skst, sken
         do ij = 1, nxy
            SDTX(ij, k, 14) = sn14pom(ij, k)
            SDTX(ij, k, 15) = sn14car(ij, k)
            SDTX(ij, k, 16) = cn14dic(ij, k)
         end do
      end do
      end if
end if

!     -- Time step --
      if ( wotstp ) then
         write(nfsots, *), 'Time Step :', TT, 'step', YY, 'year'
      end if

!     -- Conclude --

!CCCCCif ( TT == tend ) then
      if ( TT == tend .OR. mod(TT, tout) == 0 ) then

if (use_mpi) then

      if (SDofirst) then
      allocate(nrcvct2d(0:nprocs-1))
      allocate(ndispt2d(0:nprocs-1))
      allocate(nrcvct3d(0:nprocs-1))
      allocate(ndispt3d(0:nprocs-1))
      allocate(nrcvcw3d(0:nprocs-1))
      allocate(ndispw3d(0:nprocs-1))

      allocate(DUMMY01(nxyzdm+nprocs))
      allocate(DUMMY02(nxyzdm+nprocs))
      allocate(DUMMY03(nxyzdm+nprocs))
      allocate(DUMMY04(nxyzdm+nprocs))
      allocate(DUMMY05(nxyzdm+nprocs))
      allocate(DUMMY11(nxyzdm+nprocs))
      allocate(DUMMY12(nxyzdm+nprocs))
      allocate(DUMMY13(nxyzdm+nprocs))
      allocate(DUMMY14(nxyzdm+nprocs))
      allocate(DUMMY15(nxyzdm+nprocs))
      allocate(DUMMY16(nxyzdm+nprocs))
      allocate(DUMMY21(nxyzdm+nprocs))
      allocate(DUMMY22(nxydim+nprocs))
      allocate(DUMMY23(nxyzdm+nxydim+nprocs))
      allocate(DUMMY31(nxydim+nprocs))
      allocate(DUMMY32(nxydim+nprocs))
      allocate(DUMMY33(nxydim+nprocs))
      allocate(DUMMY34(nxydim+nprocs))
      allocate(DUMMY41(nxydim+nprocs))
      allocate(DUMMY42(nxydim+nprocs))
      allocate(DUMMY43(nxydim+nprocs))
      allocate(DUMMY44(nxydim+nprocs))
      allocate(DUMMY51(nxydim+nprocs))
      allocate(DUMMY52(nxydim+nprocs))
      allocate(DUMMY53(nxydim+nprocs))
      allocate(DUMMY54(nxydim+nprocs))
      if ( CALc13 ) then
      allocate(DUMMY06(nxyzdm+nprocs))
      allocate(DUMMY08(nxyzdm+nprocs))
      allocate(DUMMY17(nxyzdm+nprocs))
      allocate(DUMMY19(nxyzdm+nprocs))
      allocate(DUMMY35(nxydim+nprocs))
      allocate(DUMMY37(nxydim+nprocs))
      allocate(DUMMY45(nxydim+nprocs))
      allocate(DUMMY47(nxydim+nprocs))
      allocate(DUMMY55(nxydim+nprocs))
      allocate(DUMMY57(nxydim+nprocs))
      end if
      if ( CALc14 ) then
      allocate(DUMMY07(nxyzdm+nprocs))
      allocate(DUMMY09(nxyzdm+nprocs))
      allocate(DUMMY18(nxyzdm+nprocs))
      allocate(DUMMY20(nxyzdm+nprocs))
      allocate(DUMMY36(nxydim+nprocs))
      allocate(DUMMY38(nxydim+nprocs))
      allocate(DUMMY46(nxydim+nprocs))
      allocate(DUMMY48(nxydim+nprocs))
      allocate(DUMMY56(nxydim+nprocs))
      allocate(DUMMY58(nxydim+nprocs))
      end if

      nsndct2d = nxy
      do irank = 0, nprocs-1
         nrcvct2d(irank) = nxy
         ndispt2d(irank) = nxy * irank
      end do
      nsndct3d = nxyz
      do irank = 0, nprocs-1
         nrcvct3d(irank) = nxyz
         ndispt3d(irank) = nxyz * irank
      end do
      nsndcw3d = ( nxyz + nxy )
      do irank = 0, nprocs-1
         nrcvcw3d(irank) = ( nxyz + nxy )
         ndispw3d(irank) = ( nxyz + nxy ) * irank
      end do

      SDofirst = .false.
      end if

      call mpi_gatherv( scla,    nsndct3d, mpi_real8, &
     &                  DUMMY01, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( spom,    nsndct3d, mpi_real8, &
     &                  DUMMY02, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( scar,    nsndct3d, mpi_real8, &
     &                  DUMMY03, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( sopa,    nsndct3d, mpi_real8, &
     &                  DUMMY04, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( sage,    nsndct3d, mpi_real8, &
     &                  DUMMY05, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cnit,    nsndct3d, mpi_real8, &
     &                  DUMMY11, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cdic,    nsndct3d, mpi_real8, &
     &                  DUMMY12, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( calk,    nsndct3d, mpi_real8, &
     &                  DUMMY13, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cdo2,    nsndct3d, mpi_real8, &
     &                  DUMMY14, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( csil,    nsndct3d, mpi_real8, &
     &                  DUMMY15, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cco32,   nsndct3d, mpi_real8, &
     &                  DUMMY16, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( satsil,  nsndct3d, mpi_real8, &
     &                  DUMMY21, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( satco32, nsndct2d, mpi_real8, &
     &                  DUMMY22, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( W,       nsndcw3d, mpi_real8, &
     &                  DUMMY23, nrcvcw3d, ndispw3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( ficla,   nsndct2d, mpi_real8, &
     &                  DUMMY31, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fipom,   nsndct2d, mpi_real8, &
     &                  DUMMY32, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( ficar,   nsndct2d, mpi_real8, &
     &                  DUMMY33, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fiopa,   nsndct2d, mpi_real8, &
     &                  DUMMY34, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( focla,   nsndct2d, mpi_real8, &
     &                  DUMMY41, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fopom,   nsndct2d, mpi_real8, &
     &                  DUMMY42, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( focar,   nsndct2d, mpi_real8, &
     &                  DUMMY43, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( foopa,   nsndct2d, mpi_real8, &
     &                  DUMMY44, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bfcla,   nsndct2d, mpi_real8, &
     &                  DUMMY51, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bfpom,   nsndct2d, mpi_real8, &
     &                  DUMMY52, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bfcar,   nsndct2d, mpi_real8, &
     &                  DUMMY53, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bfopa,   nsndct2d, mpi_real8, &
     &                  DUMMY54, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )

      if ( CALc13 ) then
      call mpi_gatherv( s13pom,  nsndct3d, mpi_real8, &
     &                  DUMMY06, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( s13car,  nsndct3d, mpi_real8, &
     &                  DUMMY08, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( c13dic,  nsndct3d, mpi_real8, &
     &                  DUMMY17, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cdc13,   nsndct3d, mpi_real8, &
     &                  DUMMY19, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fi13pom, nsndct2d, mpi_real8, &
     &                  DUMMY35, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fi13car, nsndct2d, mpi_real8, &
     &                  DUMMY37, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fo13pom, nsndct2d, mpi_real8, &
     &                  DUMMY45, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fo13car, nsndct2d, mpi_real8, &
     &                  DUMMY47, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bf13pom, nsndct2d, mpi_real8, &
     &                  DUMMY55, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bf13car, nsndct2d, mpi_real8, &
     &                  DUMMY57, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      end if

      if ( CALc14 ) then
      call mpi_gatherv( s14pom,  nsndct3d, mpi_real8, &
     &                  DUMMY07, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( s14car,  nsndct3d, mpi_real8, &
     &                  DUMMY09, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( c14dic,  nsndct3d, mpi_real8, &
     &                  DUMMY18, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( cdc14,   nsndct3d, mpi_real8, &
     &                  DUMMY20, nrcvct3d, ndispt3d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fi14pom, nsndct2d, mpi_real8, &
     &                  DUMMY36, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fo14pom, nsndct2d, mpi_real8, &
     &                  DUMMY46, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( fo14car, nsndct2d, mpi_real8, &
     &                  DUMMY48, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bf14pom, nsndct2d, mpi_real8, &
     &                  DUMMY56, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      call mpi_gatherv( bf14car, nsndct2d, mpi_real8, &
     &                  DUMMY58, nrcvct2d, ndispt2d, mpi_real8, &
     &                  IROOT, mpi_COMM_WORLD, IERR     )
      end if

      if ( myrank .eq. IROOT ) then
         do irank = 0, nprocs-1
            ijkg = ndispt3d(irank)
            do k = skst, sken
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  sclag(ij, k)   = DUMMY01(ijkg)
                  spomg(ij, k)   = DUMMY02(ijkg)
                  scarg(ij, k)   = DUMMY03(ijkg)
                  sopag(ij, k)   = DUMMY04(ijkg)
                  sageg(ij, k)   = DUMMY05(ijkg)
                  cnitg(ij, k)   = DUMMY11(ijkg)
                  cdicg(ij, k)   = DUMMY12(ijkg)
                  calkg(ij, k)   = DUMMY13(ijkg)
                  cdo2g(ij, k)   = DUMMY14(ijkg)
                  csilg(ij, k)   = DUMMY15(ijkg)
                  cco32g(ij, k)  = DUMMY16(ijkg)
               if ( CALc13 ) then
                  s13pomg(ij, k) = DUMMY06(ijkg)
                  s13carg(ij, k) = DUMMY08(ijkg)
                  c13dicg(ij, k) = DUMMY17(ijkg)
                  cdc13g(ij, k)  = DUMMY19(ijkg)
               end if
               if ( CALc14 ) then
                  s14pomg(ij, k) = DUMMY07(ijkg)
                  s14carg(ij, k) = DUMMY09(ijkg)
                  c14dicg(ij, k) = DUMMY18(ijkg)
                  cdc14g(ij, k)  = DUMMY20(ijkg)
               end if
                  satsilg(ij, k) = DUMMY21(ijkg)
               end do
            end do
            ijkg = ndispw3d(irank)
            do k = skst, sken+1
               do ij = jjstr(irank), jjend(irank)
                  ijkg = ijkg + 1
                  WG(ij, k)      = DUMMY23(ijkg)
               end do
            end do
            do ij = jjstr(irank), jjend(irank)
               satco32g(ij) = DUMMY22(ij)
               ficlag(ij)   = DUMMY31(ij)
               fipomg(ij)   = DUMMY32(ij)
               ficarg(ij)   = DUMMY33(ij)
               fiopag(ij)   = DUMMY34(ij)
               foclag(ij)   = DUMMY41(ij)
               fopomg(ij)   = DUMMY42(ij)
               focarg(ij)   = DUMMY43(ij)
               foopag(ij)   = DUMMY44(ij)
               bfclag(ij)   = DUMMY51(ij)
               bfpomg(ij)   = DUMMY52(ij)
               bfcarg(ij)   = DUMMY53(ij)
               bfopag(ij)   = DUMMY54(ij)
               if ( CALc13 ) then
               fi13pomg(ij) = DUMMY35(ij)
               fi13carg(ij) = DUMMY37(ij)
               fo13pomg(ij) = DUMMY45(ij)
               fo13carg(ij) = DUMMY47(ij)
               bf13pomg(ij) = DUMMY55(ij)
               bf13carg(ij) = DUMMY57(ij)
               end if
               if ( CALc14 ) then
               fi14pomg(ij) = DUMMY36(ij)
               fi14carg(ij) = DUMMY38(ij)
               fo14pomg(ij) = DUMMY46(ij)
               fo14carg(ij) = DUMMY48(ij)
               bf14pomg(ij) = DUMMY56(ij)
               bf14carg(ij) = DUMMY58(ij)
               end if
            end do
         end do

         do k = skst, sken
            do ij = 1, nxydim
               wtg(ij, k)   = sclag(ij, k) &
     &                      + ( spomg(ij, k) * molpom * 1.0d-3 ) &
     &                      + ( scarg(ij, k) * molcar * 1.0d-3 ) &
     &                      + ( sopag(ij, k) * molopa * 1.0d-3 )
               pclag(ij, k) = sclag(ij, k) &
     &                      / wtg(ij, k) * 1.0d+2
               ppomg(ij, k) = spomg(ij, k) * molpom * 1.0d-3 &
     &                      / wtg(ij, k) * 1.0d+2
               pcarg(ij, k) = scarg(ij, k) * molcar * 1.0d-3 &
     &                      / wtg(ij, k) * 1.0d+2
               popag(ij, k) = sopag(ij, k) * molopa * 1.0d-3 &
     &                      / wtg(ij, k) * 1.0d+2
            end do
         end do

      end if
      call mpi_barrier(mpi_comm_world, IERR)

else

         do k = skst, sken
            do ij = 1, nxydim
               sclag(ij, k) = scla(ij, k)
               spomg(ij, k) = spom(ij, k)
               scarg(ij, k) = scar(ij, k)
               sopag(ij, k) = sopa(ij, k)
               sageg(ij, k) = sage(ij, k)
            end do
         end do

         do k = skst, sken
            do ij = 1, nxydim
               pclag(ij, k) = pcla(ij, k)
               ppomg(ij, k) = ppom(ij, k)
               pcarg(ij, k) = pcar(ij, k)
               popag(ij, k) = popa(ij, k)
            end do
         end do

         do k = skst, sken
            do ij = 1, nxydim
                cnitg(ij, k)   = cnit(ij, k)
                cdicg(ij, k)   = cdic(ij, k)
                calkg(ij, k)   = calk(ij, k)
                cdo2g(ij, k)   = cdo2(ij, k)
                csilg(ij, k)   = csil(ij, k)
                cco32g(ij, k)  = cco32(ij, k)
            end do
         end do

         do k = skst, sken
            do ij = 1, nxydim
                satsilg(ij, k) = satsil(ij, k)
            end do
         end do

         do ij = 1, nxydim
            satco32g(ij) = satco32(ij)
         end do

         do k = skst, sken
            do ij = 1, nxydim
                WG(ij, k) = w(ij, k)
            end do
         end do

         if ( CALc13 ) then
         do k = skst, sken
            do ij = 1, nxydim
                s13pomg(ij, k) = s13pom(ij, k)
                s13carg(ij, k) = s13car(ij, k)
                c13dicg(ij, k) = c13dic(ij, k)
                cdc13g(ij, k)  = cdc13(ij, k)
            end do
         end do
         end if

         if ( CALc14 ) then
         do k = skst, sken
            do ij = 1, nxydim
                s14pomg(ij, k) = s14pom(ij, k)
                s14carg(ij, k) = s14car(ij, k)
                c14dicg(ij, k) = c14dic(ij, k)
                cdc14g(ij, k)  = cdc14(ij, k)
            end do
         end do
         end if

         do ij = 1, nxydim
            ficlag(ij)     = ficla(ij)
            fipomg(ij)     = fipom(ij)
            ficarg(ij)     = ficar(ij)
            fiopag(ij)     = fiopa(ij)
            foclag(ij)     = focla(ij)
            fopomg(ij)     = fopom(ij)
            focarg(ij)     = focar(ij)
            foopag(ij)     = foopa(ij)
            bfclag(ij)     = bfcla(ij)
            bfpomg(ij)     = bfpom(ij)
            bfcarg(ij)     = bfcar(ij)
            bfopag(ij)     = bfopa(ij)
         end do

         if ( CALc13 ) then
         do ij = 1, nxydim
            fi13pomg(ij)   = fi13pom(ij)
            fi13carg(ij)   = fi13car(ij)
            fo13pomg(ij)   = fo13pom(ij)
            fo13carg(ij)   = fo13car(ij)
            bf13pomg(ij)   = bf13pom(ij)
            bf13carg(ij)   = bf13car(ij)
         end do
         end if

         if ( CALc14 ) then
         do ij = 1, nxydim
            fi14pomg(ij)   = fi14pom(ij)
            fi14carg(ij)   = fi14car(ij)
            fo14pomg(ij)   = fo14pom(ij)
            fo14carg(ij)   = fo14car(ij)
            bf14pomg(ij)   = bf14pom(ij)
            bf14carg(ij)   = bf14car(ij)
         end do
         end if

end if

if (use_mpi) then
      if ( myrank .eq. IROOT ) then
         CFILE  = 'scla'
         CTITLE = 'Clay'
         CUNIT  = 'g/m^3'
         call write_gtoolfile(      YY,   sclag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'spom'
         CTITLE = 'OrgM'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   spomg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'scar'
         CTITLE = 'CaCO3'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   scarg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'sopa'
         CTITLE = 'SiO2'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   sopag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'sage'
         CTITLE = 'Age tracer'
         CUNIT  = 'ND'
         call write_gtoolfile(      YY,   sageg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'cnit'
         CTITLE = 'Phosphate'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   cnitg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'cdic'
         CTITLE = 'dic'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   cdicg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'calk'
         CTITLE = 'Alkalinity'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   calkg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'cdo2'
         CTITLE = 'Oxygen'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   cdo2g,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'csil'
         CTITLE = 'Silicate'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   csilg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'ficla'
         CTITLE = 'Clay Bottom Flux'
         CUNIT  = 'g/m^3*cm/s'
         call write_gtoolfile(      YY,  ficlag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'fipom'
         CTITLE = 'OrgM Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  fipomg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'ficar'
         CTITLE = 'CaCO3 Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  ficarg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'fiopa'
         CTITLE = 'SiO2 Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  fiopag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'fopom'
         CTITLE = 'OrgM Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  fopomg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'focar'
         CTITLE = 'CaCO3 Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  focarg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'foopa'
         CTITLE = 'SiO2 Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  foopag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'bfcla'
         CTITLE = 'Clay Burial'
         CUNIT  = 'g/m^3*cm/s'
         call write_gtoolfile(      YY,  bfclag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'bfpom'
         CTITLE = 'OrgM Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  bfpomg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'bfcar'
         CTITLE = 'CaCO3 Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  bfcarg,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'bfopa'
         CTITLE = 'SiO2 Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  bfopag,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,       1  )
         CFILE  = 'w'
         CTITLE = 'Sediment Rate'
         CUNIT  = 'cm/s'
         call write_gtoolfile(      YY,      WG,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim, sdzdim+1 )
         CFILE  = 'cdc13'
         CTITLE = 'delta carbon-13'
         CUNIT  = 'permil'
         call write_gtoolfile(      YY,  cdc13g,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
         CFILE  = 'cdc14'
         CTITLE = 'Delta carbon-14'
         CUNIT  = 'permil'
         call write_gtoolfile(      YY,  cdc14g,   CFILE,  CTITLE,   CUNIT, &
     &                   nxdim,   nydim,  sdzdim  )
      end if
      call mpi_barrier(mpi_comm_world, IERR)
end if

      end if ! -- End of Output

!     call SDwrite_gtoolfile(
!    I     SEDTNM,
!    I     SDTX )

      end do ! -- End of Time integration

      end subroutine digda_cprdc

end module module_digda_cprdc

! =====================================================================
      subroutine read_gtoolfile( FLX, NSIZE, CFILE)

      IMPLICIT NONE

      INTEGER NSIZE
      REAL*8  flx(NSIZE)
      CHARACTER CFILE*128
      CHARACTER CHEAD(64)*16

      OPEN(11, FILE=CFILE, FORM='UNFORMATTED')
      READ(11) CHEAD
      READ(11) FLX
      CLOSE(11)

      end subroutine read_gtoolfile

! =====================================================================
      subroutine write_gtoolfile( &
     &           YY,    DAT3,  CFILE, CTITLE,  CUNIT, &
     &           nx,    ny,    nz     )

      IMPLICIT NONE

      INTEGER I, J, K
      INTEGER nx, ny, nz
      REAL*8  DAT3(nx,ny,nz)
      CHARACTER CFILE*128, CTITLE*128, CUNIT*128
      CHARACTER CHEAD(64)*16
      CHARACTER YEAR*8
      INTEGER YY

      DATA CHEAD / 64*'' /
      write(YEAR, 900) YY
 900  FORMAT(I8.8)

      CHEAD(1) = '9010'    !!  1 : IDFM  : I10   :   GTOOL3 version                 : 9010
      CHEAD(2) = ''        !!  2 : DSET  : A16   :   Dataset name                   :
      CHEAD(3) = CTITLE    !!  3 : ITEM  : A16   :   Identification (variable) name :
      CHEAD(4) = ''        !!  4 : EDIT1 : A16   :   Edit abbreviation 01           :
      CHEAD(5) = ''        !!  5 : EDIT2 : A16   :   Edit abbreviation 02           :
      CHEAD(6) = ''        !!  6 : EDIT3 : A16   :   Edit abbreviation 03           :
      CHEAD(7) = ''        !!  7 : EDIT4 : A16   :   Edit abbreviation 04           :
      CHEAD(8) = ''        !!  8 : EDIT5 : A16   :   Edit abbreviation 05           :
      CHEAD(9) = ''        !!  9 : EDIT6 : A16   :   Edit abbreviation 06           :
      CHEAD(10) = ''       !! 10 : EDIT7 : A16   :   Edit abbreviation 07           :
      CHEAD(11) = ''       !! 11 : EDIT8 : A16   :   Edit abbreviation 08           :
      CHEAD(12) = ''       !! 12 : FNUM  : I10   :   File number                    : 0
      CHEAD(13) = ''       !! 13 : DNUM  : I10   :   Data number                    : 0
      CHEAD(14) = CTITLE   !! 14 : TITL1 : A16   :   Title 1                        :
      CHEAD(15) = ''       !! 15 : TITL2 : A16   :   Title 2                        :
      CHEAD(16) = CUNIT    !! 16 : UNIT  : A16   :   Unit                           :
      CHEAD(17) = ''       !! 17 : ETTL1 : A16   :   Edit title 01                  :
      CHEAD(18) = ''       !! 18 : ETTL2 : A16   :   Edit title 02                  :
      CHEAD(19) = ''       !! 19 : ETTL3 : A16   :   Edit title 03                  :
      CHEAD(20) = ''       !! 20 : ETTL4 : A16   :   Edit title 04                  :
      CHEAD(21) = ''       !! 21 : ETTL5 : A16   :   Edit title 05                  :
      CHEAD(22) = ''       !! 22 : ETTL6 : A16   :   Edit title 06                  :
      CHEAD(23) = ''       !! 23 : ETTL7 : A16   :   Edit title 07                  :
      CHEAD(24) = ''       !! 24 : ETTL8 : A16   :   Edit title 08                  :
      CHEAD(25) = ''       !! 25 : TIME  : I10   :   Total time                     : 0
      CHEAD(26) = ''       !! 26 : UTIM  : A16   :   Time unit                      :
!     CHEAD(27) = ''       !! 27 : DATE  : A16   :   Time (yyyymmddhhmmss)          :
      CHEAD(28) = ''       !! 28 : TDUR  : I10   :   Data representative time       : 0
!     CHEAD(29) = ''       !! 29 : AITM1 : A16   :   Axis grid identification name  :
!     CHEAD(30) = ''       !! 30 : AStr1 : I10   :   Start of axis grid number      : 1
!     CHEAD(31) = ''       !! 31 : Aend1 : I10   :   End of axis grid number        :
!     CHEAD(32) = ''       !! 32 : AITM2 : A16   :   Axis grid identification name  :
!     CHEAD(33) = ''       !! 33 : ASTR2 : I10   :   Start of axis grid number      : 1
!     CHEAD(34) = ''       !! 34 : Aend2 : I10   :   End of axis grid number        :
!     CHEAD(35) = ''       !! 35 : AITM3 : A16   :   Axis grid identification name  :
!     CHEAD(36) = ''       !! 36 : ASTR3 : I10   :   Start of axis grid number      : 1
!     CHEAD(37) = ''       !! 37 : Aend+3 : I10   :   End of axis grid number        :
!     CHEAD(38) = ''       !! 38 : DFMT  : A16   :   Data format                    : UR4
!     CHEAD(39) = ''       !! 39 : MISS  : E15.7 :   Missing value                  : -999.
!     CHEAD(40) = ''       !! 40 : DMIN  : E15.7 :   Range (minimum)                : -999.
!     CHEAD(41) = ''       !! 41 : DMAX  : E15.7 :   Range (maximum)                : -999.
!     CHEAD(42) = ''       !! 42 : DIVS  : E15.7 :   Interval (small)               : -999.
!     CHEAD(43) = ''       !! 43 : DIVL  : E15.7 :   Interval (Large)               : -999.
!     CHEAD(44) = ''       !! 44 : STYP  : I10   :   Scaling type                   : 1
      CHEAD(45) = ''       !! 45 : COPTN : A16   :   Option (character)             :
      CHEAD(46) = ''       !! 46 : IOPTN : I16   :   Option (integer)               : 0
      CHEAD(47) = ''       !! 47 : ROPTN : E16.7 :   Option (real)                  : 0
      CHEAD(48) = ''       !! 48 : DATE1 : I16   :   Time (auxiliary)               :
      CHEAD(49) = ''       !! 49 : DATE2 : A16   :   Time unit (auxiliary)          :
!     CHEAD(50) = ''       !! 50 : MEMbf1 : A16   :   Memo bf1                        :
      CHEAD(51) = ''       !! 51 : MEMO2 : A16   :   Memo O2                        :
      CHEAD(52) = ''       !! 52 : MEMO3 : A16   :   Memo O3                        :
      CHEAD(53) = ''       !! 53 : MEMO4 : A16   :   Memo O4                        :
      CHEAD(54) = ''       !! 54 : MEMO5 : A16   :   Memo O5                        :
      CHEAD(55) = ''       !! 55 : MEMO6 : A16   :   Memo O6                        :
      CHEAD(56) = ''       !! 56 : MEMO7 : A16   :   Memo O7                        :
      CHEAD(57) = ''       !! 57 : MEMO8 : A16   :   Memo O8                        :
      CHEAD(58) = ''       !! 58 : MEMO9 : A16   :   Memo O9                        :
      CHEAD(59) = ''       !! 59 : MEM1O : A16   :   Memo 1O                        :
      CHEAD(60) = ''       !! 60 : CDATE : A16   :   Data creation date             :
      CHEAD(61) = ''       !! 61 : CSIGN : A16   :   Data creator                   :
      CHEAD(62) = ''       !! 62 : MDATE : A16   :   Data change date               :
      CHEAD(63) = ''       !! 63 : MSIGN : A16   :   Data modified by               :
      CHEAD(64) = ''       !! 64 : SIZE  : I10   :   Array size                     : 0
      write(CHEAD(39), '(E15.7)') -999.0
      write(CHEAD(40), '(E15.7)') -999.0
      write(CHEAD(41), '(E15.7)') -999.0
      write(CHEAD(42), '(E15.7)') -999.0
      write(CHEAD(43), '(E15.7)') -999.0
      write(CHEAD(44), '(I10)') 1
      write(CHEAD(29), '(A16)') 'TX'
      write(CHEAD(32), '(A16)') 'TY'
      write(CHEAD(35), '(A16)') 'TZ'
      write(CHEAD(30), '(I10)') 1
      write(CHEAD(31), '(I10)') nx
      write(CHEAD(33), '(I10)') 1
      write(CHEAD(34), '(I10)') ny
      write(CHEAD(36), '(I10)') 1
      write(CHEAD(37), '(I10)') nz
      write(CHEAD(64), '(I10)') nx*ny*nz
      CHEAD(27) = '09990101 000000 '
      CHEAD(50) = '  09990101000000'
      CHEAD(27) = YEAR//' 000000 '
      CHEAD(50) = '  '//YEAR//'000000'
      CHEAD(38) = 'UR4'

      OPEN(13, FILE=CFILE, FORM='UNFORMATTED', ACCESS='APPend')
      write(13) CHEAD
      write(13) (((Real(DAT3(I,J,K)),I=1,nx),J=1,ny),K=1,nz)
      CLOSE(13)

      end subroutine write_gtoolfile

      subroutine set_parallel_range( N1, N2, NPR, IRA, IST, IEN )

      INTEGER N1, N2
      INTEGER NPR, IRA
      INTEGER IST, IEN
      INTEGER IWORK1, IWORK2

      IWORK1 = ( N2 - N1 + 1 ) / NPR
      IWORK2 = mod( N2 - N1 + 1, NPR )
      IST = IRA*IWORK1 + N1 + min( IRA, IWORK2 )
      IEN = IST + IWORK1 - 1
      if ( IWORK2 > IRA ) then
        IEN = IEN + 1
      end if

      end subroutine set_parallel_range
