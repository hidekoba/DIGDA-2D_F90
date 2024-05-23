!######################################################################

! --- information -----------------------------------------------------
!
!  Three-dimensional sediment model
!
!  HISTORY
!     '14.08.01  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_sdmain
      use module_zcfile
      use module_digda_grid
      implicit none

! ==================================================
!   Time
! ==================================================

      INTEGER   TSTR, TEND, TENDYR
      PARAMETER (TSTR   =      1)
      PARAMETER (TENDYR =      1)  !!      10 [yr]
!     PARAMETER (TENDYR =     10)  !!      10 [yr]
!     PARAMETER (TENDYR =     20)  !!      20 [yr]
!     PARAMETER (TENDYR =     30)  !!      30 [yr]
!     PARAMETER (TENDYR =     40)  !!      40 [yr]
!     PARAMETER (TENDYR =     50)  !!      50 [yr]
!     PARAMETER (TENDYR =    100)  !!     100 [yr]
!     PARAMETER (TENDYR =    200)  !!     200 [yr]
!     PARAMETER (TENDYR =    300)  !!     300 [yr]
!     PARAMETER (TENDYR =    400)  !!     400 [yr]
!     PARAMETER (TENDYR =    500)  !!     500 [yr]
!     PARAMETER (TENDYR =   1000)  !!   1,000 [yr]
!     PARAMETER (TENDYR =   2000)  !!   2,000 [yr]
!     PARAMETER (TENDYR =   3000)  !!   3,000 [yr]
!     PARAMETER (TENDYR =   4000)  !!   4,000 [yr]
!     PARAMETER (TENDYR =   5000)  !!   5,000 [yr]
!     PARAMETER (TENDYR =  10000)  !!  10,000 [yr]
!     PARAMETER (TENDYR =  20000)  !!  20,000 [yr]
!     PARAMETER (TENDYR =  30000)  !!  30,000 [yr]
!     PARAMETER (TENDYR =  40000)  !!  40,000 [yr]
!     PARAMETER (TENDYR =  50000)  !!  50,000 [yr]
!     PARAMETER (TENDYR = 100000)  !! 100,000 [yr]
!     PARAMETER (TENDYR = 200000)  !! 200,000 [yr]
!     PARAMETER (TENDYR = 300000)  !! 300,000 [yr]
!     PARAMETER (TENDYR = 400000)  !! 400,000 [yr]
!     PARAMETER (TENDYR = 500000)  !! 500,000 [yr]
      PARAMETER (TEND    = TENDYR / TSYR )
!     PARAMETER (TEND    =      1 )

      INTEGER   TOUT, TOUTYR
      PARAMETER (TOUTYR  =     120000)
      PARAMETER (TOUT    = TOUTYR / TSYR )

! ==================================================
!   Sediments
! ==================================================

      REAL*8 SDTX(NXYDIM, SDZDIM, SEDTDM)
      SAVE   SDTX
      REAL*8 SDTXBUF(NXDIM, NYDIM, SDZDIM, SEDTDM)

! --------------------------------------------------
!   Particles
!      Clay, POM, CaCO3, Opal, agetracer
! --------------------------------------------------

!-- Concentration :
!      clay [g/m^3], pom, caco3, opal [mMol/m^3], agetracer [yr]
      REAL*8 SCLA(NXY, SDZ),   SNCLA(NXY, SDZ)
      REAL*8 SPOM(NXY, SDZ),   SNPOM(NXY, SDZ)
      REAL*8 SCAR(NXY, SDZ),   SNCAR(NXY, SDZ)
      REAL*8 SOPA(NXY, SDZ),   SNOPA(NXY, SDZ)
      REAL*8 SAGE(NXY, SDZ),   SNAGE(NXY, SDZ)
      REAL*8 S13POM(NXY, SDZ), SN13POM(NXY, SDZ)
      REAL*8 S14POM(NXY, SDZ), SN14POM(NXY, SDZ)
      REAL*8 S13CAR(NXY, SDZ), SN13CAR(NXY, SDZ)
      REAL*8 S14CAR(NXY, SDZ), SN14CAR(NXY, SDZ)
      REAL*8 SCLAG(NXYDIM, SDZDIM)
      REAL*8 SPOMG(NXYDIM, SDZDIM)
      REAL*8 SCARG(NXYDIM, SDZDIM)
      REAL*8 SOPAG(NXYDIM, SDZDIM)
      REAL*8 SAGEG(NXYDIM, SDZDIM)
      REAL*8 S13POMG(NXYDIM, SDZDIM)
      REAL*8 S14POMG(NXYDIM, SDZDIM)
      REAL*8 S13CARG(NXYDIM, SDZDIM)
      REAL*8 S14CARG(NXYDIM, SDZDIM)

!-- Weight percent [%] :
      REAL*8 PCLA(NXY, SDZ)
      REAL*8 PPOM(NXY, SDZ)
      REAL*8 PCAR(NXY, SDZ)
      REAL*8 POPA(NXY, SDZ)
      REAL*8 WT(NXY, SDZ)
      REAL*8 PCLAG(NXYDIM, SDZDIM)
      REAL*8 PPOMG(NXYDIM, SDZDIM)
      REAL*8 PCARG(NXYDIM, SDZDIM)
      REAL*8 POPAG(NXYDIM, SDZDIM)
      REAL*8 WTG(NXYDIM, SDZDIM)

!-- Dissolution rate [1/s], Dissolution term [mMol/m^3/s] :
      REAL*8 RCLA(NXY, SDZ),   TRCLA(NXY, SDZ)
      REAL*8 RPOM(NXY, SDZ),   TRPOM(NXY, SDZ)
      REAL*8 RCAR(NXY, SDZ),   TRCAR(NXY, SDZ)
      REAL*8 ROPA(NXY, SDZ),   TROPA(NXY, SDZ)
      REAL*8 R13POM(NXY, SDZ), TR13POM(NXY, SDZ)
      REAL*8 R14POM(NXY, SDZ), TR14POM(NXY, SDZ)
      REAL*8 R13CAR(NXY, SDZ), TR13CAR(NXY, SDZ)
      REAL*8 R14CAR(NXY, SDZ), TR14CAR(NXY, SDZ)

! --------------------------------------------------
!   Dissolved material
!      Phosphate, Dissolved Inorganic Carbon,
!      Alkalinity, Dissolved Oxygen, Silicate
! --------------------------------------------------

!-- Concentration [mMol/m^3] :
      REAL*8 CNIT(NXY, SDZ),   CNNIT(NXY, SDZ)
      REAL*8 CDIC(NXY, SDZ),   CNDIC(NXY, SDZ)
      REAL*8 CALK(NXY, SDZ),   CNALK(NXY, SDZ)
      REAL*8 CDO2(NXY, SDZ),   CNDO2(NXY, SDZ)
      REAL*8 CSIL(NXY, SDZ),   CNSIL(NXY, SDZ)
      REAL*8 C13DIC(NXY, SDZ), CN13DIC(NXY, SDZ)
      REAL*8 C14DIC(NXY, SDZ), CN14DIC(NXY, SDZ)
      REAL*8 CDC13(NXY, SDZ),  CDC14(NXY, SDZ)
      REAL*8 CCO32(NXY, SDZ)
      REAL*8 SATCO32(NXY)
      REAL*8 SATSIL(NXY,SDZ)
      REAL*8 CNITG(NXYDIM, SDZDIM)
      REAL*8 CDICG(NXYDIM, SDZDIM)
      REAL*8 CALKG(NXYDIM, SDZDIM)
      REAL*8 CDO2G(NXYDIM, SDZDIM)
      REAL*8 CSILG(NXYDIM, SDZDIM)
      REAL*8 C13DICG(NXYDIM, SDZDIM)
      REAL*8 C14DICG(NXYDIM, SDZDIM)
      REAL*8 CDC13G(NXYDIM, SDZDIM)
      REAL*8 CDC14G(NXYDIM, SDZDIM)
      REAL*8 CCO32G(NXYDIM, SDZDIM)
      REAL*8 SATCO32G(NXYDIM)
      REAL*8 SATSILG(NXYDIM,SDZDIM)

!-- Concentration (bottom water) :
!      temp [C], salinity [psu], the others [mMol/m^3]
      REAL*8 BWTEM(NXY)
      REAL*8 BWSAL(NXY)
      REAL*8 BWNIT(NXY)
      REAL*8 BWDIC(NXY)
      REAL*8 BWALK(NXY)
      REAL*8 BWDO2(NXY)
      REAL*8 BWSIL(NXY)
      REAL*8 BW13DIC(NXY)
      REAL*8 BW14DIC(NXY)
      REAL*8 BWTEMG(NXYDIM)
      REAL*8 BWSALG(NXYDIM)
      REAL*8 BWNITG(NXYDIM)
      REAL*8 BWDICG(NXYDIM)
      REAL*8 BWALKG(NXYDIM)
      REAL*8 BWDO2G(NXYDIM)
      REAL*8 BWSILG(NXYDIM)
      REAL*8 BW13DICG(NXYDIM)
      REAL*8 BW14DICG(NXYDIM)

!-- Concentration (total) :
!      total phosphorous [mMol/m^3], total silicon [mMol/m^3]
      REAL*8 PT(NXY)
      REAL*8 SIT(NXY)

!-- Dissolution term [mMol/m^3/s] :
      REAL*8 TRNIT(NXY, SDZ)
      REAL*8 TRDIC(NXY, SDZ)
      REAL*8 TRALK(NXY, SDZ)
      REAL*8 TRDO2(NXY, SDZ)
      REAL*8 TRSIL(NXY, SDZ)
      REAL*8 TR13DIC(NXY, SDZ)
      REAL*8 TR14DIC(NXY, SDZ)

!-- Sea-bottom Flux :
!      clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
      REAL*8 FICLA(NXY),   FOCLA(NXY)
      REAL*8 FIPOM(NXY),   FOPOM(NXY)
      REAL*8 FICAR(NXY),   FOCAR(NXY)
      REAL*8 FIOPA(NXY),   FOOPA(NXY)
      REAL*8 FIAGE(NXY)
      REAL*8 FI13POM(NXY), FO13POM(NXY)
      REAL*8 FI14POM(NXY), FO14POM(NXY)
      REAL*8 FI13CAR(NXY), FO13CAR(NXY)
      REAL*8 FI14CAR(NXY), FO14CAR(NXY)
      REAL*8 FICLAG(NXYDIM), FOCLAG(NXYDIM)
      REAL*8 FIPOMG(NXYDIM), FOPOMG(NXYDIM)
      REAL*8 FICARG(NXYDIM), FOCARG(NXYDIM)
      REAL*8 FIOPAG(NXYDIM), FOOPAG(NXYDIM)
      REAL*8 FIAGEG(NXYDIM), FOAGEG(NXYDIM)
      REAL*8 FI13POMG(NXYDIM), FO13POMG(NXYDIM)
      REAL*8 FI14POMG(NXYDIM), FO14POMG(NXYDIM)
      REAL*8 FI13CARG(NXYDIM), FO13CARG(NXYDIM)
      REAL*8 FI14CARG(NXYDIM), FO14CARG(NXYDIM)

!-- Sea-bottom Flux [mMol/m^3*cm/s] :
      REAL*8 FINIT(NXY),   FONIT(NXY)
      REAL*8 FIDIC(NXY),   FODIC(NXY)
      REAL*8 FIALK(NXY),   FOALK(NXY)
      REAL*8 FIDO2(NXY),   FODO2(NXY)
      REAL*8 FISIL(NXY),   FOSIL(NXY)
      REAL*8 FI13DIC(NXY), FO13DIC(NXY)
      REAL*8 FI14DIC(NXY), FO14DIC(NXY)

!-- Sediment-bottom Flux :
!      clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
      REAL*8 BFCLA(NXY)
      REAL*8 BFPOM(NXY)
      REAL*8 BFCAR(NXY)
      REAL*8 BFOPA(NXY)
      REAL*8 BF13POM(NXY)
      REAL*8 BF14POM(NXY)
      REAL*8 BF13CAR(NXY)
      REAL*8 BF14CAR(NXY)
      REAL*8 BFCLAG(NXYDIM)
      REAL*8 BFPOMG(NXYDIM)
      REAL*8 BFCARG(NXYDIM)
      REAL*8 BFOPAG(NXYDIM)
      REAL*8 BF13POMG(NXYDIM)
      REAL*8 BF14POMG(NXYDIM)
      REAL*8 BF13CARG(NXYDIM)
      REAL*8 BF14CARG(NXYDIM)

!-- Density, Molar mass :
!      Density [g/cm^3], Molar mass [g/Mol],
!      Molar mass / Density [cm^3/Mol]
      REAL*8 RHOCLA, RHOPOM, RHOCAR, RHOOPA
      REAL*8 MOLCLA, MOLPOM, MOLCAR, MOLOPA
      REAL*8 MCLA,   MPOM,   MCAR,   MOPA
      SAVE   RHOCLA, RHOPOM, RHOCAR, RHOOPA
      SAVE   MOLCLA, MOLPOM, MOLCAR, MOLOPA
      SAVE   MCLA,   MPOM,   MCAR,   MOPA
      DATA   RHOCLA, RHOPOM, RHOCAR, RHOOPA &
     &    / 2.60D0, 2.60D0, 2.60D0, 2.60D0 /
          !! 2.6, 2.6, 2.6, 2.6 [g/cm^3]
      DATA   MOLCLA, MOLPOM, MOLCAR, MOLOPA &
     &    / 430.0D0, 12.0D0, 100.0D0,  60.0D0 /
          !! 430, 12, 100, 60 [g/Mol]
      REAL*8 RHO13POM, RHO14POM, RHO13CAR, RHO14CAR
      REAL*8 MOL13POM, MOL14POM, MOL13CAR, MOL14CAR
      REAL*8 M13POM,   M14POM,   M13CAR,   M14CAR
      SAVE   RHO13POM, RHO14POM, RHO13CAR, RHO14CAR
      SAVE   MOL13POM, MOL14POM, MOL13CAR, MOL14CAR
      SAVE   M13POM,   M14POM,   M13CAR,   M14CAR
      DATA   RHO13POM, RHO14POM, RHO13CAR, RHO14CAR &
     &    / 2.60D0, 2.60D0, 2.60D0, 2.60D0 /
          !! 2.6, 2.6, 2.6, 2.6 [g/cm^3]
      DATA   MOL13POM, MOL14POM, MOL13CAR, MOL14CAR &
     &    /  13.0D0, 14.0D0, 101.0D0, 102.0D0 /
          !! 13, 14, 101, 102 [g/Mol]

!-- Redfield ratio :
      REAL*8 FP, PP, NP, CP, OP
          !! Fe:P, P:P, N:P, C:P, O:P
      SAVE   FP, PP, NP, CP, OP
      DATA   FP, PP, NP, CP, OP &
     &    / 4.70D-4, 1.0D0, 16.0D0, 117.0D0, 170.0D0 /
          !! 0.00047, 1, 16, 117, 170

!-- Advection, Diffusion :
!      Sedimentation rate
      REAL*8 W(NXY, SDZ+1)
      REAL*8 WG(NXYDIM, SDZDIM+1)
!      Porosity
      REAL*8 PHI0, Z0
      SAVE   PHI0, Z0
      DATA   PHI0, Z0 &
     &    / 0.9D0, 1.823D0 /
          !! 0.9, 1.823 [m]
      REAL*8 PHI(SDZ), PHII(SDZ+1)
      SAVE   PHI, PHII
!      Bioturbation diffusivity
      REAL*8 DB0, ZB
      SAVE   DB0, ZB
      DATA   DB0, ZB &
     &    / 3.0D-5, 1.0D-1 /
          !! 3.0*10-5 [m^2/yr], 0.10 [m]
      REAL*8 DBI(SDZ+1)
      SAVE   DBI
!      Pw Effective molecular diffusivity
      REAL*8 DSW
      REAL*8 DSWNIT, DSWDIC, DSWALK, DSWDO2, DSWSIL
      REAL*8 DSW13DIC, DSW14DIC
      REAL*8 NOPA
      SAVE   DSW
      SAVE   DSWNIT, DSWDIC, DSWALK, DSWDO2, DSWSIL
      SAVE   DSW13DIC, DSW14DIC
      SAVE   NOPA
      DATA   DSW &
     &    / 4.59D-10 /
          !! 4.59*10-10 [m^2/s]
      DATA   DSWNIT, DSWDIC, DSWALK, DSWDO2, DSWSIL &
     &    / 5.0D-10, 6.4D-10, 6.4D-10, 1.2D-9, 5.0D-10 /
          !! 5.0*10-10, 6.4*10-10, 6.4*10-10, 1.2*10-9, 5.0*10-10  [m^2/s]
      DATA   DSW13DIC, DSW14DIC &
     &    / 6.4D-10, 6.4D-10 /
          !! 6.4*10-10, 6.4*10-10 [m^2/s]
      DATA   NOPA &
     &    / 2.50D0 /
      REAL*8 DEINIT(SDZ+1)
      REAL*8 DEIDIC(SDZ+1)
      REAL*8 DEIALK(SDZ+1)
      REAL*8 DEIDO2(SDZ+1)
      REAL*8 DEISIL(SDZ+1)
      REAL*8 DEI13DIC(SDZ+1)
      REAL*8 DEI14DIC(SDZ+1)
      SAVE   DEINIT, DEIDIC, DEIALK, DEIDO2, DEISIL
      SAVE   DEI13DIC, DEI14DIC

!-- 13C, 14C :
      REAL*8 PDB              !! standard 13C/12C
      REAL*8 DSTD             !! standard 14C/12C
      REAL*8 EOPH, EOSH, EPS  !! fractionaton factor
      REAL*8 FEOPH, FEOSH, R13, R14
      REAL*8 THALF, ADCAY     !! half decay time [yr], decay constant
      REAL*8 DLTC13, DLTC14
      SAVE   PDB
      SAVE   DSTD
      SAVE   EOPH, EOSH
      SAVE   FEOPH, FEOSH
      SAVE   THALF, ADCAY
      DATA   PDB &
     &    / 0.0112372D0 /
      DATA   DSTD &
     &    / 1.1892D-12 /
      DATA   EOPH, EOSH &
     &    / -22.0, 0.0 /
      DATA   EPS &
     &    / 1.0D-10 /
      DATA   THALF &
     &    / 5370.D0 /

!-- MASK :
      REAL*8 ZW(NZ)
      SAVE ZW
      REAL*8 BDEP(NXY), BDEPG(NXYDIM) !! ocean depth [m]
      SAVE BDEP, BDEPG

!-- Parallel computing :
      INTEGER, ALLOCATABLE :: NRCVC1(:),  NDISP1(:)
      INTEGER, ALLOCATABLE :: NRCVC2(:),  NDISP2(:)
      INTEGER, ALLOCATABLE :: NRCVC3(:),  NDISP3(:)
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
      LOGICAL SDFIRST
      SAVE SDFIRST
      DATA SDFIRST / .TRUE. /
      LOGICAL SDOFIRST
      SAVE SDOFIRST
      DATA SDOFIRST / .TRUE.  /

      LOGICAL DEIARC
      SAVE DEIARC
      NAMELIST /NMSDDEAR/ DEIARC
      DATA DEIARC / .TRUE. /  !! Archer's DEI formulation ?

      LOGICAL  CALC13, CALC14
      SAVE     CALC13, CALC14
      NAMELIST /NMSDISO/ CALC13, CALC14
      DATA CALC13, CALC14 / .TRUE., .TRUE. /  !! carbon isotopes ?

!-- etc :
      INTEGER I, J, IJ, K, L, IC
      CHARACTER SEDTNM(SEDTDM)*16
      SAVE SEDTNM
      INTEGER TT

      LOGICAL  WOTSTP
      SAVE     WOTSTP
      DATA     WOTSTP / .TRUE. /  !! writeout time step ?
      INTEGER     NFSOTS
      PARAMETER ( NFSOTS = 999 )

      CHARACTER CFILE*128, CTITLE*128, CUNIT*128

end module globals_sdmain

! =====================================================================

module module_sdmain
      use globals_sdmain
      use module_digda_clay
      use module_digda_pom
      use module_digda_caco3
      use module_digda_sio2l
      use module_digda_step_solid
      use module_digda_step_liquid
      use module_digda_chemeq
      implicit none
contains

      SUBROUTINE SD

!-- SEDIMENT TRACERS : setup --

IF (PARALLEL) THEN
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

      IF ( NPROCS /= NCPU ) THEN
         IF ( MYRANK == 0 ) THEN
            PRINT *, '=== ERROR ==='
            CALL MPI_FINALIZE(IERR)
            STOP
         END IF
      END IF

      DO IRANK = 0, NPROCS-1
         CALL PARA_RANGE( &
     &                    1, NXYDIM, NPROCS, IRANK, &
     &                    JJSTR(IRANK), JJEND(IRANK) &
     &                  )
      END DO
END IF

      OPEN(NFPARM, file="PARAMET")
      OPEN(NFSTDO, file="STDOUT")

      WRITE(NFSTDO, *) '*** SDMAINSET ***'
      WRITE(NFSTDO, *) '-----------------------------'
      WRITE(NFSTDO, *) '-- Setup of Sediment Model --'
      WRITE(NFSTDO, *) '-----------------------------'
      REWIND(NFPARM)
      READ (NFPARM, NMSDDEAR, END=99)
99 continue
      WRITE(NFSTDO, NMSDDEAR)
      REWIND(NFPARM)
      READ (NFPARM, NMSDISO, END=98)
98 continue
      WRITE(NFSTDO, NMSDISO)

!     ### for Sediment Model ###
      SEDTNM(1)  = 'SCLA'
      SEDTNM(2)  = 'SPOM'
      SEDTNM(3)  = 'SCAR'
      SEDTNM(4)  = 'SOPA'
      SEDTNM(5)  = 'SAGE'
      SEDTNM(6)  = 'CNIT'
      SEDTNM(7)  = 'CDIC'
      SEDTNM(8)  = 'CALK'
      SEDTNM(9)  = 'CSIL'
      SEDTNM(10) = 'CDO2'
      SEDTNM(11) = 'S13POM'
      SEDTNM(12) = 'S13CAR'
      SEDTNM(13) = 'C13DIC'
      SEDTNM(14) = 'S14POM'
      SEDTNM(15) = 'S14CAR'
      SEDTNM(16) = 'C14DIC'

!CCCCCCALL SDFILESET
      CALL SDCLAYSET
      CALL SDPOMSET
      CALL SDCACO3SET
      CALL SDOPALSET
      CALL SDCHEMEQSET
      CALL SDSOLIDSET
      CALL SDLIQUIDSET

!-- SEDIMENT TRACERS : setup (end) --

!     -- Ocean depth --
      CFILE="bound/bdep"
      CALL RDFILE(  BDEPG, NXYDIM,  CFILE )
IF (PARALLEL) THEN
      DO IJ = 1, NXY
         BDEP(IJ) = MAX( BDEPG(MYRANK*NXY+IJ), 0.0D0 )
      END DO
ELSE
      DO IJ = 1, NXY
         BDEP(IJ) = MAX( BDEPG(IJ), 0.0D0 )
      END DO
END IF

!     -- Convert unit --

      Z0  = Z0  * 1.0D+2                       !! [m] --> [cm]

      DB0 = DB0 * 1.0D+4 / 3.1536D+7           !! [m^2/yr] --> [cm^2/s]
      ZB  = ZB  * 1.0D+2                       !! [m] --> [cm]

      DO K = SKST, SKEN
         LZ(K)   = LZ(K)   * 1.0D+2            !! [m] --> [cm]
         DLZ(K)  = DLZ(K)  * 1.0D+2            !! [m] --> [cm]
      END DO
      DO K = SKST, SKEN+1
         LZI(K)  = LZI(K)  * 1.0D+2            !! [m] --> [cm]
         DLZI(K) = DLZI(K) * 1.0D+2            !! [m] --> [cm]
      END DO

      DSWNIT   = DSWNIT   * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      DSWDIC   = DSWDIC   * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      DSWALK   = DSWALK   * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      DSWDO2   = DSWDO2   * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      DSWSIL   = DSWSIL   * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      IF ( CALC13 ) THEN
      DSW13DIC = DSW13DIC * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      END IF
      IF ( CALC14 ) THEN
      DSW14DIC = DSW14DIC * 1.0D+4             !! [m^2/s] --> [cm^2/s]
      END IF
      DSW    = DSW        * 1.0D+4             !! [m^2/s] --> [cm^2/s]

      MCLA   = 1.0D-6 /  RHOCLA                !! [m^3/g]
      MPOM   = 1.0D-9 / ( RHOPOM / MOLPOM )    !! [m^3/mMol]
      MCAR   = 1.0D-9 / ( RHOCAR / MOLCAR )    !! [m^3/mMol]
      MOPA   = 1.0D-9 / ( RHOOPA / MOLOPA )    !! [m^3/mMol]
      IF ( CALC13 ) THEN
      M13POM = 1.0D-9 / ( RHOPOM / MOL13POM )  !! [m^3/mMol]
      M13CAR = 1.0D-9 / ( RHOCAR / MOL13CAR )  !! [m^3/mMol]
      END IF
      IF ( CALC14 ) THEN
      M14POM = 1.0D-9 / ( RHOPOM / MOL14POM )  !! [m^3/mMol]
      M14CAR = 1.0D-9 / ( RHOCAR / MOL14CAR )  !! [m^3/mMol]
      THALF  = THALF * 3.1536D+7               !! [yr] --> [s]
      ADCAY  = LOG( 2.0D0 ) / THALF / 365.0D0 / 8.64D4
      END IF
      IF ( CALC13 .or. CALC14 ) THEN
      FEOPH  = 1.0D0 + EOPH * 1.0D-3
      FEOSH  = 1.0D0 + EOSH * 1.0D-3
      END IF

      DO IJ = 1, NXY
         BDEP(IJ) = BDEP(IJ) * 1.0D+2          !! [m] --> [cm]
      END DO

!     -- Porosity, Diffusivity --

      DO K = SKST, SKEN
         PHI(K) = PHI0 * EXP( - LZ(K)  / Z0 )
      END DO

      DO K = SKST, SKEN+1
         PHII(K) = PHI0 * EXP( - LZI(K) / Z0 )
         DBI(K)  = ( DB0 / ( 1.0D0 + EXP( 2.0D0 * ( LZI(K) - ZB ) ) ) )
      END DO

      DO K = SKST, SKEN+1
         IF (DEIARC) THEN
            DEINIT(K)   = DSWNIT   * ( PHII(K) ** 3.0D0 )
            DEIDIC(K)   = DSWDIC   * ( PHII(K) ** 3.0D0 )
            DEIALK(K)   = DSWALK   * ( PHII(K) ** 3.0D0 )
            DEIDO2(K)   = DSWDO2   * ( PHII(K) ** 3.0D0 )
            DEISIL(K)   = DSWSIL   * ( PHII(K) ** 3.0D0 )
            IF ( CALC13 ) THEN
            DEI13DIC(K) = DSW13DIC * ( PHII(K) ** 3.0D0 )
            END IF
            IF ( CALC14 ) THEN
            DEI14DIC(K) = DSW14DIC * ( PHII(K) ** 3.0D0 )
            END IF
         ELSE
            DEINIT(K)   = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            DEIDIC(K)   = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            DEIALK(K)   = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            DEIDO2(K)   = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            DEISIL(K)   = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            IF ( CALC13 ) THEN
            DEI13DIC(K) = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            END IF
            IF ( CALC14 ) THEN
            DEI14DIC(K) = ( PHII(K) ** ( NOPA - 1.0D0 ) ) * DSW
            END IF
         END IF
      END DO

!     -- Set boundary conditions --

      CFILE="bound/ficla"
      CALL RDFILE( FICLAG, NXYDIM,  CFILE )
      CFILE="bound/fipom"
      CALL RDFILE( FIPOMG, NXYDIM,  CFILE )
      CFILE="bound/ficar"
      CALL RDFILE( FICARG, NXYDIM,  CFILE )
      CFILE="bound/fiopa"
      CALL RDFILE( FIOPAG, NXYDIM,  CFILE )
IF (PARALLEL) THEN
      DO IJ = 1, NXY
         FICLA(IJ) = MAX( FICLAG(MYRANK*NXY+IJ), 0.0D0 )
         FIPOM(IJ) = MAX( FIPOMG(MYRANK*NXY+IJ), 0.0D0 )
         FICAR(IJ) = MAX( FICARG(MYRANK*NXY+IJ), 0.0D0 )
         FIOPA(IJ) = MAX( FIOPAG(MYRANK*NXY+IJ), 0.0D0 )
         FIAGE(IJ) = 0.0D0
      END DO
ELSE
      DO IJ = 1, NXY
         FICLA(IJ) = MAX( FICLAG(IJ), 0.0D0 )
         FIPOM(IJ) = MAX( FIPOMG(IJ), 0.0D0 )
         FICAR(IJ) = MAX( FICARG(IJ), 0.0D0 )
         FIOPA(IJ) = MAX( FIOPAG(IJ), 0.0D0 )
         FIAGE(IJ) = 0.0D0
      END DO
END IF
      CFILE="bound/btem"
      CALL RDFILE( BWTEMG, NXYDIM,  CFILE )
      CFILE="bound/bsal"
      CALL RDFILE( BWSALG, NXYDIM,  CFILE )
      CFILE="bound/bnit"
      CALL RDFILE( BWNITG, NXYDIM,  CFILE )
      CFILE="bound/bdic"
      CALL RDFILE( BWDICG, NXYDIM,  CFILE )
      CFILE="bound/balk"
      CALL RDFILE( BWALKG, NXYDIM,  CFILE )
      CFILE="bound/bdo2"
      CALL RDFILE( BWDO2G, NXYDIM,  CFILE )
      CFILE="bound/bsil"
      CALL RDFILE( BWSILG, NXYDIM,  CFILE )
IF (PARALLEL) THEN
      DO IJ = 1, NXY
         BWTEM(IJ)   = BWTEMG(MYRANK*NXY+IJ)
         BWSAL(IJ)   = BWSALG(MYRANK*NXY+IJ)
         BWNIT(IJ)   = MAX( BWNITG(MYRANK*NXY+IJ), 0.0D0 )
         BWDIC(IJ)   = MAX( BWDICG(MYRANK*NXY+IJ), 0.0D0 )
         BWALK(IJ)   = MAX( BWALKG(MYRANK*NXY+IJ), 0.0D0 )
         BWDO2(IJ)   = MAX( BWDO2G(MYRANK*NXY+IJ), 0.0D0 )
         BWSIL(IJ)   = MAX( BWSILG(MYRANK*NXY+IJ), 0.0D0 )
      END DO
ELSE
      DO IJ = 1, NXY
         BWTEM(IJ)   = BWTEMG(IJ)
         BWSAL(IJ)   = BWSALG(IJ)
         BWNIT(IJ)   = MAX( BWNITG(IJ), 0.0D0 )
         BWDIC(IJ)   = MAX( BWDICG(IJ), 0.0D0 )
         BWALK(IJ)   = MAX( BWALKG(IJ), 0.0D0 )
         BWDO2(IJ)   = MAX( BWDO2G(IJ), 0.0D0 )
         BWSIL(IJ)   = MAX( BWSILG(IJ), 0.0D0 )
      END DO
END IF

      IF ( CALC13 ) THEN
      DO IJ = IJTSTR, IJTEND
         FI13POM(IJ) = FEOPH * PDB  * FIPOM(IJ)
         FI13CAR(IJ) = FEOSH * PDB  * FICAR(IJ)
      END DO
      CFILE="bound/bc13"
      CALL RDFILE(  BW13DICG, NXYDIM,  CFILE )
IF (PARALLEL) THEN
      DO IJ = 1, NXY
         BW13DIC(IJ) = MAX( BW13DICG(MYRANK*NXY+IJ), 0.0D0 )
      END DO
ELSE
      DO IJ = 1, NXY
         BW13DIC(IJ) = MAX( BW13DICG(IJ), 0.0D0 )
      END DO
END IF
      END IF

      IF ( CALC14 ) THEN
      DO IJ = IJTSTR, IJTEND
         FI14POM(IJ) = FEOPH * DSTD * FIPOM(IJ)
         FI14CAR(IJ) = FEOSH * DSTD * FICAR(IJ)
      END DO
      CFILE="bound/bc14"
      CALL RDFILE(  BW14DICG, NXYDIM,  CFILE )
IF (PARALLEL) THEN
      DO IJ = 1, NXY
         BW14DIC(IJ) = MAX( BW14DICG(MYRANK*NXY+IJ), 0.0D0 )
      END DO
ELSE
      DO IJ = 1, NXY
         BW14DIC(IJ) = MAX( BW14DICG(IJ), 0.0D0 )
      END DO
END IF
      END IF

      DO IJ = 1, NXY
         SIT(IJ) = 0.51D0
         PT(IJ)  = 7.70D0
      END DO

!     DO L = 1, SEDTDM
!        DO K = 1, SDZDIM
!           DO J = 1, NYDIM
!              DO I = 1, NXDIM
!                 IJ = I + (J - 1) * NXDIM
!                 SDTXBUF(I, J, K, L) = SDTX(IJ, K, L)
!              END DO
!           END DO
!        END DO
!     END DO
!     DO L = 1, SEDTDM
!        DO K = SKST, SKEN
!           DO J = JSTR, JEND
!              DO I = 1, ISTR-1
!                 SDTXBUF(I, J, K, L) = SDTXBUF(IEND-ISTR+1+I, J, K, L)
!              END DO
!              DO I = IEND+1, NXDIM
!                 SDTXBUF(I, J, K, L) = SDTXBUF(ISTR-1+I-IEND, J, K, L)
!              END DO
!           END DO
!        END DO
!     END DO
!     DO L = 1, SEDTDM
!        DO K = 1, SDZDIM
!           DO J = 1, NYDIM
!              DO I = 1, NXDIM
!                 IJ = I + (J - 1) * NXDIM
!                 SDTX(IJ, K, L) = SDTXBUF(I, J, K, L)
!              END DO
!           END DO
!        END DO
!     END DO
!     DO L = 1, SEDTDM
!        DO K = 1, SKST-1
!           DO IJ = 1, NXY
!              SDTX(IJ, K, L) = SDTX(IJ, SKST, L)
!           END DO
!        END DO
!     END DO

!     -- Time Integration --

      DO TT = TSTR, TEND

         YY = INT( TSYR * TT )

      IF (SDFIRST) THEN
         DO K = SKST, SKEN
            DO IJ = IJTSTR, IJTEND
               SCLA(IJ, K) = RHOCLA * 1.0D+6
               SPOM(IJ, K) = 0.0D0
               SCAR(IJ, K) = 0.0D0
               SOPA(IJ, K) = 0.0D0
               SAGE(IJ, K) = 0.0D0
               CNIT(IJ, K) = BWNIT(IJ)
               CDIC(IJ, K) = BWDIC(IJ)
               CALK(IJ, K) = BWALK(IJ)
               CDO2(IJ, K) = BWDO2(IJ)
               CSIL(IJ, K) = BWSIL(IJ)
               IF ( CALC13 ) THEN
               S13POM(IJ, K) = 0.0D0
               S13CAR(IJ, K) = 0.0D0
               C13DIC(IJ, K) = BW13DIC(IJ)
               END IF
               IF ( CALC14 ) THEN
               S14POM(IJ, K) = 0.0D0
               S14CAR(IJ, K) = 0.0D0
               C14DIC(IJ, K) = BW14DIC(IJ)
               END IF
               PCLA(IJ, K) = 1.00+2
               PPOM(IJ, K) = 0.0D0
               PCAR(IJ, K) = 0.0D0
               POPA(IJ, K) = 0.0D0
            END DO
         END DO
         SDFIRST = .FALSE.
      ELSE
IF (PARALLEL) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               SCLA(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  1)
               SPOM(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  2)
               SCAR(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  3)
               SOPA(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  4)
               SAGE(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  5)
               CNIT(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  6)
               CDIC(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  7)
               CALK(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  8)
               CSIL(IJ, K) = SDTX(MYRANK*NXY+IJ, K,  9)
               CDO2(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 10)
            END DO
         END DO
         IF ( CALC13 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               S13POM(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 11)
               S13CAR(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 12)
               C13DIC(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 13)
            END DO
         END DO
         END IF
         IF ( CALC14 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               S14POM(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 14)
               S14CAR(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 15)
               C14DIC(IJ, K) = SDTX(MYRANK*NXY+IJ, K, 16)
            END DO
         END DO
         END IF
ELSE
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               SCLA(IJ, K) = SDTX(IJ, K,  1)
               SPOM(IJ, K) = SDTX(IJ, K,  2)
               SCAR(IJ, K) = SDTX(IJ, K,  3)
               SOPA(IJ, K) = SDTX(IJ, K,  4)
               SAGE(IJ, K) = SDTX(IJ, K,  5)
               CNIT(IJ, K) = SDTX(IJ, K,  6)
               CDIC(IJ, K) = SDTX(IJ, K,  7)
               CALK(IJ, K) = SDTX(IJ, K,  8)
               CSIL(IJ, K) = SDTX(IJ, K,  9)
               CDO2(IJ, K) = SDTX(IJ, K, 10)
            END DO
         END DO
         IF ( CALC13 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               S13POM(IJ, K) = SDTX(IJ, K, 11)
               S13CAR(IJ, K) = SDTX(IJ, K, 12)
               C13DIC(IJ, K) = SDTX(IJ, K, 13)
            END DO
         END DO
         END IF
         IF ( CALC14 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXY
               S14POM(IJ, K) = SDTX(IJ, K, 14)
               S14CAR(IJ, K) = SDTX(IJ, K, 15)
               C14DIC(IJ, K) = SDTX(IJ, K, 16)
            END DO
         END DO
         END IF
END IF
      END IF

!     -- Saturated carbonate ion concentration --

      CALL SDCHEMEQ( &
     &    BWTEM, BWSAL, CDIC, CALK, &
     &    SIT, PT, &
     &    CCO32 )

      CALL SATCAR( &
     &    BWTEM, BWSAL, BDEP, &
     &    SATCO32 )

!     -- Dissolution term --

      CALL DISCLA( &
     &    SCLA, &
     &    RCLA )

      CALL DISPOM( &
     &    SPOM, CDO2, CP, OP, PHI, &
     &    RPOM )

      CALL DISCAR( &
     &    SCAR, CCO32, SATCO32, PHI, &
     &    RCAR )

      CALL DISOPA( &
     &    SOPA, CSIL, PCLA, POPA, BWTEM, &
     &    ROPA, SATSIL )

!     -- Dissolution rate --
!        CLA, POM, CAR, OPA

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TRCLA(IJ, K)   = ( 1.0D0 - PHI(K) ) * RCLA(IJ, K) &
     &                     * SCLA(IJ, K)
            TRPOM(IJ, K)   = ( 1.0D0 - PHI(K) ) * RPOM(IJ, K) &
     &                     * SPOM(IJ, K)
            TRCAR(IJ, K)   = ( 1.0D0 - PHI(K) ) * RCAR(IJ, K) &
     &                     * SCAR(IJ, K)
            TROPA(IJ, K)   = ( 1.0D0 - PHI(K) ) * ROPA(IJ, K) &
     &                     * SOPA(IJ, K)
         END DO
      END DO

      IF ( CALC13 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TR13POM(IJ, K) = ( 1.0D0 - PHI(K) ) * RPOM(IJ, K) &
     &                     * S13POM(IJ, K)
            TR13CAR(IJ, K) = ( 1.0D0 - PHI(K) ) * RCAR(IJ, K) &
     &                     * S13CAR(IJ, K)
         END DO
      END DO
      END IF

      IF ( CALC14 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TR14POM(IJ, K) = ( 1.0D0 - PHI(K) ) * RPOM(IJ, K) &
     &                     * S14POM(IJ, K)
            TR14CAR(IJ, K) = ( 1.0D0 - PHI(K) ) * RCAR(IJ, K) &
     &                     * S14CAR(IJ, K)
         END DO
      END DO
      END IF

!     -- Sedimentation rate --

      DO IJ = IJTSTR, IJTEND
         W(IJ, SKST) = ( ( ( FICLA(IJ) * MCLA ) &
     &                   + ( FIPOM(IJ) * MPOM ) &
     &                   + ( FICAR(IJ) * MCAR ) &
     &                   + ( FIOPA(IJ) * MOPA ) ) &
     &               / ( 1.0D0 - PHII(SKST) ) )
         DO K = SKST, SKEN
            W(IJ, K+1) = ( ( 1.0D0 - PHII(K) ) / ( 1.0D0 - PHII(K+1) ) &
     &                   * W(IJ, K) ) &
     &                 - ( DLZ(K) * ( ( TRPOM(IJ, K) * MPOM ) &
     &                              + ( TRCAR(IJ, K) * MCAR ) &
     &                              + ( TROPA(IJ, K) * MOPA ) ) &
     &                   / ( 1.0D0 - PHII(K+1) ) )
         END DO
      END DO

!     -- Particles (t --> t+dt) --
!        Clay, POM, CaCO3, Opal

      CALL SOLID_CLAY( &
     &    FICLA, TRCLA, W, PHI, PHII, DBI, SCLA, RHOCLA, &
     &    SNCLA, BFCLA )

      CALL SOLID( &
     &    FIPOM, TRPOM, W, PHI, PHII, DBI, SPOM, &
     &    SNPOM, BFPOM )

      CALL SOLID( &
     &    FICAR, TRCAR, W, PHI, PHII, DBI, SCAR, &
     &    SNCAR, BFCAR )

      CALL SOLID( &
     &    FIOPA, TROPA, W, PHI, PHII, DBI, SOPA, &
     &    SNOPA, BFOPA )

      CALL AGE( &
     &    FIAGE, W, PHI, PHII, DBI, SAGE, &
     &    SNAGE )

      IF ( CALC13 ) THEN
      CALL SOLID( &
     &    FI13POM, TR13POM, W, PHI, PHII, DBI, S13POM, &
     &    SN13POM, BF13POM )

      CALL SOLID( &
     &    FI13CAR, TR13CAR, W, PHI, PHII, DBI, S13CAR, &
     &    SN13CAR, BF13CAR )
      END IF

      IF ( CALC14 ) THEN
      CALL SOLID( &
     &    FI14POM, TR14POM, W, PHI, PHII, DBI, S14POM, &
     &    SN14POM, BF14POM )

      CALL SOLID( &
     &    FI14CAR, TR14CAR, W, PHI, PHII, DBI, S14CAR, &
     &    SN14CAR, BF14CAR )
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            SN14POM(IJ, K) = SN14POM(IJ, K) &
     &                     - ADCAY * TS * S14POM(IJ, K)
            SN14CAR(IJ, K) = SN14CAR(IJ, K) &
     &                     - ADCAY * TS * S14CAR(IJ, K)
         END DO
      END DO
      END IF

!     -- Dissolution rate --
!        NIT, DIC, ALK, DO2, SIL

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TRNIT(IJ, K)   = ( ( PP / CP ) * TRPOM(IJ, K) )
            TRDIC(IJ, K)   = ( TRPOM(IJ, K) + TRCAR(IJ, K) )
            TRALK(IJ, K)   = ( - ( NP / CP ) * TRPOM(IJ, K) &
     &                         + ( 2.0D0 * TRCAR(IJ, K) ) )
            TRDO2(IJ, K)   = - ( ( OP / CP ) * TRPOM(IJ, K) )
            IF ( CDO2(IJ, K) <= 0.0D0 ) THEN
               TRDO2(IJ, K) = 0.0D0
            END IF
            TRSIL(IJ, K)   = TROPA(IJ, K)
         END DO
      END DO

      IF ( CALC13 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TR13DIC(IJ, K) = ( TR13POM(IJ, K) + TR13CAR(IJ, K) )
         END DO
      END DO
      END IF

      IF ( CALC14 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            TR14DIC(IJ, K) = ( TR14POM(IJ, K) + TR14CAR(IJ, K) )
         END DO
      END DO
      END IF

!     -- Dissolved material (t --> t+dt) --
!        NIT, DIC, ALK, DO2, SIL

      CALL LIQUID( &
     &    BWNIT, TRNIT, PHI, PHII, DEINIT, CNIT, &
     &    CNNIT )

      CALL LIQUID( &
     &    BWDIC, TRDIC, PHI, PHII, DEIDIC, CDIC, &
     &    CNDIC )

      CALL LIQUID( &
     &    BWALK, TRALK, PHI, PHII, DEIALK, CALK, &
     &    CNALK )

      CALL LIQUID_O2( &
     &    BWDO2, TRDO2, PHI, PHII, DEIDO2, CDO2, &
     &    CNDO2 )

      CALL LIQUID( &
     &    BWSIL, TRSIL, PHI, PHII, DEISIL, CSIL, &
     &    CNSIL )

      IF ( CALC13 ) THEN
      CALL LIQUID( &
     &    BW13DIC, TR13DIC, PHI, PHII, DEI13DIC, C13DIC, &
     &    CN13DIC )
      END IF

      IF ( CALC14 ) THEN
      CALL LIQUID( &
     &    BW14DIC, TR14DIC, PHI, PHII, DEI14DIC, C14DIC, &
     &    CN14DIC )
      END IF
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            CN14DIC(IJ, K) = CN14DIC(IJ, K) &
     &                     - ADCAY * TS * C14DIC(IJ, K)
         END DO
      END DO

!     -- Particles weight percent --

      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            WT(IJ, K) = SNCLA(IJ, K) &
     &                + ( SNPOM(IJ, K) * MOLPOM * 1.0D-3 ) &
     &                + ( SNCAR(IJ, K) * MOLCAR * 1.0D-3 ) &
     &                + ( SNOPA(IJ, K) * MOLOPA * 1.0D-3 )
            IF ( WT(IJ, K) == 0.0D0 ) THEN
               PCLA(IJ, K) = 0.0D0
               PPOM(IJ, K) = 0.0D0
               PCAR(IJ, K) = 0.0D0
               POPA(IJ, K) = 0.0D0
            ELSE
               PCLA(IJ, K) = &
     &            ( SNCLA(IJ, K) / WT(IJ, K) ) &
     &            * 1.0D+2
               PPOM(IJ, K) = &
     &            ( SNPOM(IJ, K) * MOLPOM * 1.0D-3 / WT(IJ, K) ) &
     &            * 1.0D+2
               PCAR(IJ, K) = &
     &            ( SNCAR(IJ, K) * MOLCAR * 1.0D-3 / WT(IJ, K) ) &
     &            * 1.0D+2
               POPA(IJ, K) = &
     &            ( SNOPA(IJ, K) * MOLOPA * 1.0D-3 / WT(IJ, K) ) &
     &            * 1.0D+2
            END IF
         END DO
      END DO

!     -- Isotopic fractionation --

      IF ( CALC13 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            IF ( CNDIC(IJ, K) /= 0.0D0 ) THEN
               R13 = CN13DIC(IJ, K) / MAX( CNDIC(IJ, K), EPS )
               CDC13(IJ, K) = &
     &            ( 1.0D0 / 0.0112372D0 &
     &              * R13 / ( 1.0D0 - R13 ) - 1.0D0 ) * 1.0D3
            ELSE
               CDC13(IJ, K) = 0.0D0
            END IF
         END DO
      END DO
      END IF

      IF ( CALC14 ) THEN
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            IF ( CNDIC(IJ, K) /= 0.0D0 ) THEN
               R14 = CN14DIC(IJ, K) / MAX( CNDIC(IJ, K), EPS )
               CDC14(IJ, K) = &
     &            ( 1.0D0 / 1.1892D-12 &
     &              * R14 / ( 1.0D0 - R14 ) - 1.0D0 ) * 1.0D3
            ELSE
               CDC14(IJ, K) = 0.0D0
            END IF
         END DO
      END DO
      DO K = SKST, SKEN
         DO IJ = IJTSTR, IJTEND
            DLTC13 = CDC13(IJ, K) * 1.0D-3  !! [permil]
            DLTC14 = CDC14(IJ, K) * 1.0D-3  !! [permil]
            CDC14(IJ, K) = &
     &         ( DLTC14 - 2.0D0 * ( DLTC13 + 0.025D0 ) &
     &                          * ( 1.0D0 + DLTC14 ) &
     &         ) * 1.0D+3   !! [permil]
         END DO
      END DO
      END IF

!     -- Sea-bottom flux --

      DO IJ = IJTSTR, IJTEND
         FONIT(IJ) = - ( PHII(SKST) * DEINIT(SKST) &
     &      * ( BWNIT(IJ) - CNNIT(IJ, SKST) ) / DLZI(SKST) )
         FODIC(IJ) = - ( PHII(SKST) * DEIDIC(SKST) &
     &      * ( BWDIC(IJ) - CNDIC(IJ, SKST) ) / DLZI(SKST) )
         FOALK(IJ) = - ( PHII(SKST) * DEIALK(SKST) &
     &      * ( BWALK(IJ) - CNALK(IJ, SKST) ) / DLZI(SKST) )
         FODO2(IJ) = - ( PHII(SKST) * DEIDO2(SKST) &
     &      * ( BWDO2(IJ) - CNDO2(IJ, SKST) ) / DLZI(SKST) )
         FOSIL(IJ) = - ( PHII(SKST) * DEISIL(SKST) &
     &      * ( BWSIL(IJ) - CNSIL(IJ, SKST) ) / DLZI(SKST) )
         FINIT(IJ)  = FIPOM(IJ) * ( PP / CP )
         FIDIC(IJ)  = FIPOM(IJ) + FICAR(IJ)
         FIALK(IJ)  = - FIPOM(IJ) * ( NP / CP ) + 2.0D0 * FICAR(IJ)
         FIDO2(IJ)  = FIPOM(IJ) * ( OP / CP )
         FISIL(IJ)  = FIOPA(IJ)
         IF ( CALC13 ) THEN
         FO13DIC(IJ) = - ( PHII(SKST) * DEI13DIC(SKST) &
     &      * ( BW13DIC(IJ) - CN13DIC(IJ, SKST) ) / DLZI(SKST) )
         FI13DIC(IJ)  = FI13POM(IJ) + FI13CAR(IJ)
         END IF
         IF ( CALC14 ) THEN
         FO14DIC(IJ) = - ( PHII(SKST) * DEI14DIC(SKST) &
     &      * ( BW14DIC(IJ) - CN14DIC(IJ, SKST) ) / DLZI(SKST) )
         FI14DIC(IJ)  = FI14POM(IJ) + FI14CAR(IJ)
         END IF
      END DO

      DO IJ = IJTSTR, IJTEND
         FOPOM(IJ) = CP * FONIT(IJ)
         FOCAR(IJ) = 0.5D0 * ( + NP * FONIT(IJ) + FOALK(IJ) )
!         FOPOM(IJ) = ( CP / ( NP + 2.0D0 * CP ) ) &
!     &      * ( 2.0D0 * FODIC(IJ) - FOALK(IJ) )
!         FOCAR(IJ) = ( CP / ( NP + 2.0D0 * CP ) ) &
!     &      * ( ( NP / CP ) * FODIC(IJ) + FOALK(IJ) )
         FOOPA(IJ) = FOSIL(IJ)
      END DO

!     -- OUTPUT --
IF (PARALLEL) THEN
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(MYRANK*NXY+IJ, K,  1) = SNCLA(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  2) = SNPOM(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  3) = SNCAR(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  4) = SNOPA(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  5) = SNAGE(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  6) = CNNIT(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  7) = CNDIC(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  8) = CNALK(IJ, K)
            SDTX(MYRANK*NXY+IJ, K,  9) = CNSIL(IJ, K)
            SDTX(MYRANK*NXY+IJ, K, 10) = CNDO2(IJ, K)
         END DO
      END DO
      IF ( CALC13 ) THEN
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(MYRANK*NXY+IJ, K, 11) = SN13POM(IJ, K)
            SDTX(MYRANK*NXY+IJ, K, 12) = SN13CAR(IJ, K)
            SDTX(MYRANK*NXY+IJ, K, 13) = CN13DIC(IJ, K)
         END DO
      END DO
      END IF
      IF ( CALC14 ) THEN
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(MYRANK*NXY+IJ, K, 14) = SN14POM(IJ, K)
            SDTX(MYRANK*NXY+IJ, K, 15) = SN14CAR(IJ, K)
            SDTX(MYRANK*NXY+IJ, K, 16) = CN14DIC(IJ, K)
         END DO
      END DO
      END IF
ELSE
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(IJ, K,  1) = SNCLA(IJ, K)
            SDTX(IJ, K,  2) = SNPOM(IJ, K)
            SDTX(IJ, K,  3) = SNCAR(IJ, K)
            SDTX(IJ, K,  4) = SNOPA(IJ, K)
            SDTX(IJ, K,  5) = SNAGE(IJ, K)
            SDTX(IJ, K,  6) = CNNIT(IJ, K)
            SDTX(IJ, K,  7) = CNDIC(IJ, K)
            SDTX(IJ, K,  8) = CNALK(IJ, K)
            SDTX(IJ, K,  9) = CNSIL(IJ, K)
            SDTX(IJ, K, 10) = CNDO2(IJ, K)
         END DO
      END DO
      IF ( CALC13 ) THEN
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(IJ, K, 11) = SN13POM(IJ, K)
            SDTX(IJ, K, 12) = SN13CAR(IJ, K)
            SDTX(IJ, K, 13) = CN13DIC(IJ, K)
         END DO
      END DO
      END IF
      IF ( CALC14 ) THEN
      DO K = SKST, SKEN
         DO IJ = 1, NXY
            SDTX(IJ, K, 14) = SN14POM(IJ, K)
            SDTX(IJ, K, 15) = SN14CAR(IJ, K)
            SDTX(IJ, K, 16) = CN14DIC(IJ, K)
         END DO
      END DO
      END IF
END IF

!     -- Time step --
      IF ( WOTSTP ) THEN
         WRITE(NFSOTS, *), 'Time Step :', TT, 'step', YY, 'year'
      END IF

!     -- Conclude --

!CCCCCIF ( TT == TEND ) THEN
      IF ( TT == TEND .OR. MOD(TT, TOUT) == 0 ) THEN

IF (PARALLEL) THEN

      IF (SDOFIRST) THEN
      ALLOCATE(NRCVC1(0:NPROCS-1))
      ALLOCATE(NDISP1(0:NPROCS-1))
      ALLOCATE(NRCVC2(0:NPROCS-1))
      ALLOCATE(NDISP2(0:NPROCS-1))
      ALLOCATE(NRCVC3(0:NPROCS-1))
      ALLOCATE(NDISP3(0:NPROCS-1))

      ALLOCATE(DUMMY01(NXYZDM+NPROCS))
      ALLOCATE(DUMMY02(NXYZDM+NPROCS))
      ALLOCATE(DUMMY03(NXYZDM+NPROCS))
      ALLOCATE(DUMMY04(NXYZDM+NPROCS))
      ALLOCATE(DUMMY05(NXYZDM+NPROCS))
      ALLOCATE(DUMMY11(NXYZDM+NPROCS))
      ALLOCATE(DUMMY12(NXYZDM+NPROCS))
      ALLOCATE(DUMMY13(NXYZDM+NPROCS))
      ALLOCATE(DUMMY14(NXYZDM+NPROCS))
      ALLOCATE(DUMMY15(NXYZDM+NPROCS))
      ALLOCATE(DUMMY16(NXYZDM+NPROCS))
      ALLOCATE(DUMMY21(NXYZDM+NPROCS))
      ALLOCATE(DUMMY22(NXYDIM+NPROCS))
      ALLOCATE(DUMMY23(NXYZDM+NXYDIM+NPROCS))
      ALLOCATE(DUMMY31(NXYDIM+NPROCS))
      ALLOCATE(DUMMY32(NXYDIM+NPROCS))
      ALLOCATE(DUMMY33(NXYDIM+NPROCS))
      ALLOCATE(DUMMY34(NXYDIM+NPROCS))
      ALLOCATE(DUMMY41(NXYDIM+NPROCS))
      ALLOCATE(DUMMY42(NXYDIM+NPROCS))
      ALLOCATE(DUMMY43(NXYDIM+NPROCS))
      ALLOCATE(DUMMY44(NXYDIM+NPROCS))
      ALLOCATE(DUMMY51(NXYDIM+NPROCS))
      ALLOCATE(DUMMY52(NXYDIM+NPROCS))
      ALLOCATE(DUMMY53(NXYDIM+NPROCS))
      ALLOCATE(DUMMY54(NXYDIM+NPROCS))
      IF ( CALC13 ) THEN
      ALLOCATE(DUMMY06(NXYZDM+NPROCS))
      ALLOCATE(DUMMY08(NXYZDM+NPROCS))
      ALLOCATE(DUMMY17(NXYZDM+NPROCS))
      ALLOCATE(DUMMY19(NXYZDM+NPROCS))
      ALLOCATE(DUMMY35(NXYDIM+NPROCS))
      ALLOCATE(DUMMY37(NXYDIM+NPROCS))
      ALLOCATE(DUMMY45(NXYDIM+NPROCS))
      ALLOCATE(DUMMY47(NXYDIM+NPROCS))
      ALLOCATE(DUMMY55(NXYDIM+NPROCS))
      ALLOCATE(DUMMY57(NXYDIM+NPROCS))
      END IF
      IF ( CALC14 ) THEN
      ALLOCATE(DUMMY07(NXYZDM+NPROCS))
      ALLOCATE(DUMMY09(NXYZDM+NPROCS))
      ALLOCATE(DUMMY18(NXYZDM+NPROCS))
      ALLOCATE(DUMMY20(NXYZDM+NPROCS))
      ALLOCATE(DUMMY36(NXYDIM+NPROCS))
      ALLOCATE(DUMMY38(NXYDIM+NPROCS))
      ALLOCATE(DUMMY46(NXYDIM+NPROCS))
      ALLOCATE(DUMMY48(NXYDIM+NPROCS))
      ALLOCATE(DUMMY56(NXYDIM+NPROCS))
      ALLOCATE(DUMMY58(NXYDIM+NPROCS))
      END IF

      ISIZE1 = NXY
      DO IRANK = 0, NPROCS-1
         NRCVC1(IRANK) = NXY
         NDISP1(IRANK) = NXY * IRANK
      END DO
      ISIZE2 = NXYZ
      DO IRANK = 0, NPROCS-1
         NRCVC2(IRANK) = NXYZ
         NDISP2(IRANK) = NXYZ * IRANK
      END DO
      ISIZE3 = ( NXYZ + NXY )
      DO IRANK = 0, NPROCS-1
         NRCVC3(IRANK) = ( NXYZ + NXY )
         NDISP3(IRANK) = ( NXYZ + NXY ) * IRANK
      END DO

      SDOFIRST = .FALSE.
      END IF

      CALL MPI_GATHERV( SCLA,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY01, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SPOM,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY02, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SCAR,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY03, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SOPA,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY04, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SAGE,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY05, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CNIT,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY11, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CDIC,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY12, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CALK,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY13, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CDO2,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY14, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CSIL,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY15, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CCO32,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY16, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SATSIL,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY21, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( SATCO32, ISIZE1,         MPI_REAL8, &
     &                  DUMMY22, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( W,       ISIZE3,         MPI_REAL8, &
     &                  DUMMY23, NRCVC3, NDISP3, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FICLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY31, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FIPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY32, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FICAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY33, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FIOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY34, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FOCLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY41, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FOPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY42, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FOCAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY43, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FOOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY44, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BFCLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY51, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BFPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY52, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BFCAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY53, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BFOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY54, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )

      IF ( CALC13 ) THEN
      CALL MPI_GATHERV( S13POM,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY06, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( S13CAR,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY08, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( C13DIC,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY17, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CDC13,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY19, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FI13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY35, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FI13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY37, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FO13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY45, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FO13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY47, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BF13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY55, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BF13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY57, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      END IF

      IF ( CALC14 ) THEN
      CALL MPI_GATHERV( S14POM,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY07, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( S14CAR,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY09, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( C14DIC,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY18, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( CDC14,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY20, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FI14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY36, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FO14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY46, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( FO14CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY48, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BF14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY56, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      CALL MPI_GATHERV( BF14CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY58, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      END IF

      IF ( MYRANK .EQ. IROOT ) THEN
         DO IRANK = 0, NPROCS-1
            IJKG = NDISP2(IRANK)
            DO K = SKST, SKEN
               DO IJ = JJSTR(IRANK), JJEND(IRANK)
                  IJKG = IJKG + 1
                  SCLAG(IJ, K)   = DUMMY01(IJKG)
                  SPOMG(IJ, K)   = DUMMY02(IJKG)
                  SCARG(IJ, K)   = DUMMY03(IJKG)
                  SOPAG(IJ, K)   = DUMMY04(IJKG)
                  SAGEG(IJ, K)   = DUMMY05(IJKG)
                  CNITG(IJ, K)   = DUMMY11(IJKG)
                  CDICG(IJ, K)   = DUMMY12(IJKG)
                  CALKG(IJ, K)   = DUMMY13(IJKG)
                  CDO2G(IJ, K)   = DUMMY14(IJKG)
                  CSILG(IJ, K)   = DUMMY15(IJKG)
                  CCO32G(IJ, K)  = DUMMY16(IJKG)
               IF ( CALC13 ) THEN
                  S13POMG(IJ, K) = DUMMY06(IJKG)
                  S13CARG(IJ, K) = DUMMY08(IJKG)
                  C13DICG(IJ, K) = DUMMY17(IJKG)
                  CDC13G(IJ, K)  = DUMMY19(IJKG)
               END IF
               IF ( CALC14 ) THEN
                  S14POMG(IJ, K) = DUMMY07(IJKG)
                  S14CARG(IJ, K) = DUMMY09(IJKG)
                  C14DICG(IJ, K) = DUMMY18(IJKG)
                  CDC14G(IJ, K)  = DUMMY20(IJKG)
               END IF
                  SATSILG(IJ, K) = DUMMY21(IJKG)
               END DO
            END DO
            IJKG = NDISP3(IRANK)
            DO K = SKST, SKEN+1
               DO IJ = JJSTR(IRANK), JJEND(IRANK)
                  IJKG = IJKG + 1
                  WG(IJ, K)      = DUMMY23(IJKG)
               END DO
            END DO
            DO IJ = JJSTR(IRANK), JJEND(IRANK)
               SATCO32G(IJ) = DUMMY22(IJ)
               FICLAG(IJ)   = DUMMY31(IJ)
               FIPOMG(IJ)   = DUMMY32(IJ)
               FICARG(IJ)   = DUMMY33(IJ)
               FIOPAG(IJ)   = DUMMY34(IJ)
               FOCLAG(IJ)   = DUMMY41(IJ)
               FOPOMG(IJ)   = DUMMY42(IJ)
               FOCARG(IJ)   = DUMMY43(IJ)
               FOOPAG(IJ)   = DUMMY44(IJ)
               BFCLAG(IJ)   = DUMMY51(IJ)
               BFPOMG(IJ)   = DUMMY52(IJ)
               BFCARG(IJ)   = DUMMY53(IJ)
               BFOPAG(IJ)   = DUMMY54(IJ)
               IF ( CALC13 ) THEN
               FI13POMG(IJ) = DUMMY35(IJ)
               FI13CARG(IJ) = DUMMY37(IJ)
               FO13POMG(IJ) = DUMMY45(IJ)
               FO13CARG(IJ) = DUMMY47(IJ)
               BF13POMG(IJ) = DUMMY55(IJ)
               BF13CARG(IJ) = DUMMY57(IJ)
               END IF
               IF ( CALC14 ) THEN
               FI14POMG(IJ) = DUMMY36(IJ)
               FI14CARG(IJ) = DUMMY38(IJ)
               FO14POMG(IJ) = DUMMY46(IJ)
               FO14CARG(IJ) = DUMMY48(IJ)
               BF14POMG(IJ) = DUMMY56(IJ)
               BF14CARG(IJ) = DUMMY58(IJ)
               END IF
            END DO
         END DO

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
               WTG(IJ, K)   = SCLAG(IJ, K) &
     &                      + ( SPOMG(IJ, K) * MOLPOM * 1.0D-3 ) &
     &                      + ( SCARG(IJ, K) * MOLCAR * 1.0D-3 ) &
     &                      + ( SOPAG(IJ, K) * MOLOPA * 1.0D-3 )
               PCLAG(IJ, K) = SCLAG(IJ, K) &
     &                      / WTG(IJ, K) * 1.0D+2
               PPOMG(IJ, K) = SPOMG(IJ, K) * MOLPOM * 1.0D-3 &
     &                      / WTG(IJ, K) * 1.0D+2
               PCARG(IJ, K) = SCARG(IJ, K) * MOLCAR * 1.0D-3 &
     &                      / WTG(IJ, K) * 1.0D+2
               POPAG(IJ, K) = SOPAG(IJ, K) * MOLOPA * 1.0D-3 &
     &                      / WTG(IJ, K) * 1.0D+2
            END DO
         END DO

      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

ELSE

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
               SCLAG(IJ, K) = SCLA(IJ, K)
               SPOMG(IJ, K) = SPOM(IJ, K)
               SCARG(IJ, K) = SCAR(IJ, K)
               SOPAG(IJ, K) = SOPA(IJ, K)
               SAGEG(IJ, K) = SAGE(IJ, K)
            END DO
         END DO

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
               PCLAG(IJ, K) = PCLA(IJ, K)
               PPOMG(IJ, K) = PPOM(IJ, K)
               PCARG(IJ, K) = PCAR(IJ, K)
               POPAG(IJ, K) = POPA(IJ, K)
            END DO
         END DO

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
                CNITG(IJ, K)   = CNIT(IJ, K)
                CDICG(IJ, K)   = CDIC(IJ, K)
                CALKG(IJ, K)   = CALK(IJ, K)
                CDO2G(IJ, K)   = CDO2(IJ, K)
                CSILG(IJ, K)   = CSIL(IJ, K)
                CCO32G(IJ, K)  = CCO32(IJ, K)
            END DO
         END DO

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
                SATSILG(IJ, K) = SATSIL(IJ, K)
            END DO
         END DO

         DO IJ = 1, NXYDIM
            SATCO32G(IJ) = SATCO32(IJ)
         END DO

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
                WG(IJ, K) = W(IJ, K)
            END DO
         END DO

         IF ( CALC13 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
                S13POMG(IJ, K) = S13POM(IJ, K)
                S13CARG(IJ, K) = S13CAR(IJ, K)
                C13DICG(IJ, K) = C13DIC(IJ, K)
                CDC13G(IJ, K)  = CDC13(IJ, K)
            END DO
         END DO
         END IF

         IF ( CALC14 ) THEN
         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
                S14POMG(IJ, K) = S14POM(IJ, K)
                S14CARG(IJ, K) = S14CAR(IJ, K)
                C14DICG(IJ, K) = C14DIC(IJ, K)
                CDC14G(IJ, K)  = CDC14(IJ, K)
            END DO
         END DO
         END IF

         DO IJ = 1, NXYDIM
            FICLAG(IJ)     = FICLA(IJ)
            FIPOMG(IJ)     = FIPOM(IJ)
            FICARG(IJ)     = FICAR(IJ)
            FIOPAG(IJ)     = FIOPA(IJ)
            FOCLAG(IJ)     = FOCLA(IJ)
            FOPOMG(IJ)     = FOPOM(IJ)
            FOCARG(IJ)     = FOCAR(IJ)
            FOOPAG(IJ)     = FOOPA(IJ)
            BFCLAG(IJ)     = BFCLA(IJ)
            BFPOMG(IJ)     = BFPOM(IJ)
            BFCARG(IJ)     = BFCAR(IJ)
            BFOPAG(IJ)     = BFOPA(IJ)
         END DO

         IF ( CALC13 ) THEN
         DO IJ = 1, NXYDIM
            FI13POMG(IJ)   = FI13POM(IJ)
            FI13CARG(IJ)   = FI13CAR(IJ)
            FO13POMG(IJ)   = FO13POM(IJ)
            FO13CARG(IJ)   = FO13CAR(IJ)
            BF13POMG(IJ)   = BF13POM(IJ)
            BF13CARG(IJ)   = BF13CAR(IJ)
         END DO
         END IF

         IF ( CALC14 ) THEN
         DO IJ = 1, NXYDIM
            FI14POMG(IJ)   = FI14POM(IJ)
            FI14CARG(IJ)   = FI14CAR(IJ)
            FO14POMG(IJ)   = FO14POM(IJ)
            FO14CARG(IJ)   = FO14CAR(IJ)
            BF14POMG(IJ)   = BF14POM(IJ)
            BF14CARG(IJ)   = BF14CAR(IJ)
         END DO
         END IF

END IF

IF (PARALLEL) THEN
      IF ( MYRANK .EQ. IROOT ) THEN
         CFILE  = 'scla'
         CTITLE = 'Clay'
         CUNIT  = 'g/m^3'
         CALL FILEOUT(      YY,   SCLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'spom'
         CTITLE = 'OrgM'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   SPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'scar'
         CTITLE = 'CaCO3'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   SCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'sopa'
         CTITLE = 'SiO2'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   SOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'sage'
         CTITLE = 'Age tracer'
         CUNIT  = 'ND'
         CALL FILEOUT(      YY,   SAGEG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cnit'
         CTITLE = 'Phosphate'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   CNITG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdic'
         CTITLE = 'DIC'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   CDICG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'calk'
         CTITLE = 'Alkalinity'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   CALKG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdo2'
         CTITLE = 'Oxygen'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   CDO2G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'csil'
         CTITLE = 'Silicate'
         CUNIT  = 'mmol/m^3'
         CALL FILEOUT(      YY,   CSILG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'ficla'
         CTITLE = 'Clay Bottom Flux'
         CUNIT  = 'g/m^3*cm/s'
         CALL FILEOUT(      YY,  FICLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fipom'
         CTITLE = 'OrgM Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FIPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'ficar'
         CTITLE = 'CaCO3 Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FICARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fiopa'
         CTITLE = 'SiO2 Bottom Flux'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FIOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fopom'
         CTITLE = 'OrgM Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FOPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'focar'
         CTITLE = 'CaCO3 Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FOCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'foopa'
         CTITLE = 'SiO2 Diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  FOOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfcla'
         CTITLE = 'Clay Burial'
         CUNIT  = 'g/m^3*cm/s'
         CALL FILEOUT(      YY,  BFCLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'bfpom'
         CTITLE = 'OrgM Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  BFPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfcar'
         CTITLE = 'CaCO3 Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  BFCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfopa'
         CTITLE = 'SiO2 Burial'
         CUNIT  = 'mMol/m^3*cm/s'
         CALL FILEOUT(      YY,  BFOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'w'
         CTITLE = 'Sediment Rate'
         CUNIT  = 'cm/s'
         CALL FILEOUT(      YY,      WG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM, SDZDIM+1 )
         CFILE  = 'cdc13'
         CTITLE = 'delta carbon-13'
         CUNIT  = 'permil'
         CALL FILEOUT(      YY,  CDC13G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdc14'
         CTITLE = 'Delta carbon-14'
         CUNIT  = 'permil'
         CALL FILEOUT(      YY,  CDC14G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
END IF

      END IF ! -- End of Output

!     CALL SDFILEOUT(
!    I     SEDTNM,
!    I     SDTX )

      END DO ! -- End of Time integration

      end subroutine SD

end module module_sdmain

! =====================================================================
      SUBROUTINE RDFILE( FLX, NSIZE, CFILE)

      IMPLICIT NONE

      INTEGER NSIZE
      REAL*8  FLX(NSIZE)
      CHARACTER CFILE*128
      CHARACTER CHEAD(64)*16

      OPEN(11, FILE=CFILE, FORM='UNFORMATTED')
      READ(11) CHEAD
      READ(11) FLX
      CLOSE(11)

      END SUBROUTINE RDFILE

! =====================================================================
      SUBROUTINE FILEOUT( &
     &           YY,    DAT3,  CFILE, CTITLE,  CUNIT, &
     &           NX,    NY,    NZ     )

      IMPLICIT NONE

      INTEGER I, J, K
      INTEGER NX, NY, NZ
      REAL*8  DAT3(NX,NY,NZ)
      CHARACTER CFILE*128, CTITLE*128, CUNIT*128
      CHARACTER CHEAD(64)*16
      CHARACTER YEAR*8
      INTEGER YY

      DATA CHEAD / 64*'' /
      WRITE(YEAR, 900) YY
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
!     CHEAD(30) = ''       !! 30 : ASTR1 : I10   :   Start of axis grid number      : 1
!     CHEAD(31) = ''       !! 31 : AEND1 : I10   :   End of axis grid number        :
!     CHEAD(32) = ''       !! 32 : AITM2 : A16   :   Axis grid identification name  :
!     CHEAD(33) = ''       !! 33 : ASTR2 : I10   :   Start of axis grid number      : 1
!     CHEAD(34) = ''       !! 34 : AEND2 : I10   :   End of axis grid number        :
!     CHEAD(35) = ''       !! 35 : AITM3 : A16   :   Axis grid identification name  :
!     CHEAD(36) = ''       !! 36 : ASTR3 : I10   :   Start of axis grid number      : 1
!     CHEAD(37) = ''       !! 37 : AEND3 : I10   :   End of axis grid number        :
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
!     CHEAD(50) = ''       !! 50 : MEMO1 : A16   :   Memo O1                        :
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
      WRITE(CHEAD(39), '(E15.7)') -999.0
      WRITE(CHEAD(40), '(E15.7)') -999.0
      WRITE(CHEAD(41), '(E15.7)') -999.0
      WRITE(CHEAD(42), '(E15.7)') -999.0
      WRITE(CHEAD(43), '(E15.7)') -999.0
      WRITE(CHEAD(44), '(I10)') 1
      WRITE(CHEAD(29), '(A16)') 'TX'
      WRITE(CHEAD(32), '(A16)') 'TY'
      WRITE(CHEAD(35), '(A16)') 'TZ'
      WRITE(CHEAD(30), '(I10)') 1
      WRITE(CHEAD(31), '(I10)') NX
      WRITE(CHEAD(33), '(I10)') 1
      WRITE(CHEAD(34), '(I10)') NY
      WRITE(CHEAD(36), '(I10)') 1
      WRITE(CHEAD(37), '(I10)') NZ
      WRITE(CHEAD(64), '(I10)') NX*NY*NZ
      CHEAD(27) = '09990101 000000 '
      CHEAD(50) = '  09990101000000'
      CHEAD(27) = YEAR//' 000000 '
      CHEAD(50) = '  '//YEAR//'000000'
      CHEAD(38) = 'UR4'

      OPEN(13, FILE=CFILE, FORM='UNFORMATTED', ACCESS='APPEND')
      WRITE(13) CHEAD
      WRITE(13) (((REAL(DAT3(I,J,K)),I=1,NX),J=1,NY),K=1,NZ)
      CLOSE(13)

      END SUBROUTINE FILEOUT

      SUBROUTINE PARA_RANGE( N1, N2, NPR, IRA, IST, IEN )

      INTEGER N1, N2
      INTEGER NPR, IRA
      INTEGER IST, IEN
      INTEGER IWORK1, IWORK2

      IWORK1 = ( N2 - N1 + 1 ) / NPR
      IWORK2 = MOD( N2 - N1 + 1, NPR )
      IST = IRA*IWORK1 + N1 + MIN( IRA, IWORK2 )
      IEN = IST + IWORK1 - 1
      IF ( IWORK2 > IRA ) THEN
        IEN = IEN + 1
      END IF

      END SUBROUTINE PARA_RANGE
