!######################################################################

! --- information -----------------------------------------------------
!
!  Three-dimensional sediment model
!
!  HISTORY
!     '14.08.01  H.Kobayashi :
!
! ---------------------------------------------------------------------

module globals_digda_main
  use module_zcfile
  use module_digda_grid
  use module_digda_parameter
  use module_digda_timestep
  use module_digda_mpi
  implicit none

! ! Porosity
! double precision :: PHI(SDZ)
! double precision :: PHII(SDZ+1)

! ! Bioturbation diffusivity
! double precision :: DBI(SDZ+1)

! ! Porewater diffusivity
! double precision :: DEINIT(SDZ+1)
! double precision :: DEIDIC(SDZ+1)
! double precision :: DEIALK(SDZ+1)
! double precision :: DEIDO2(SDZ+1)
! double precision :: DEISIL(SDZ+1)
! double precision :: DEI13DIC(SDZ+1)
! double precision :: DEI14DIC(SDZ+1)

  ! Carbon isotopes
  double precision :: FEOPH, FEOSH, R13, R14, ADCAY, THALF, DLTC13, DLTC14

  ! Option --
  logical :: SDFIRST = .true.
  logical :: SDOFIRST = .true.

  logical :: DEIARC = .true.  !! Archer's DEI formulation ?
  namelist /NMSDDEAR/ DEIARC

  logical :: CALC13 = .true. !! calculate carbon-13 ?
  logical :: CALC14 = .true. !! calculate carbon-14 ?
  namelist /NMSDISO/ CALC13, CALC14

  ! etc
  character(16) :: SEDTNM(SEDTDM)
  character(128) :: CFILE, CTITLE, CUNIT

contains

end module globals_digda_main

! =====================================================================

module module_digda_main
  use netcdf_utilities
  use gtool_utilities
  use globals_digda_main
  use module_digda_clay
  use module_digda_pom
  use module_digda_carbonate
  use module_digda_sio2l
  use module_digda_step_solid
  use module_digda_step_liquid
  use module_digda_chemeq
  implicit none

contains

  subroutine digda_main( &
 &  bdepg, &
 &  bfrcg, &
 &  sdtx )

    ! --------------------------------------------------
    !   Forcing of Particles and Dissolved material
    !      Clay, POM, CaCO3, Opal, agetracer
    !      Phosphate, Dissolved Inorganic Carbon,
    !      Alkalinity, Dissolved Oxygen, Silicate
    ! --------------------------------------------------

    double precision, intent(in) :: BDEPG(NXYDIM)
    double precision, intent(in) :: BFRCG(NXYDIM, BENFDM)

    double precision, intent(inout) :: SDTX(NXYDIM, SDZDIM, SEDTDM)
    double precision :: SDTXBUF(NXDIM, NYDIM, SDZDIM, SEDTDM)

    ! --------------------------------------------------
    !   Particles
    !     Clay, POM, CaCO3, Opal, agetracer
    ! --------------------------------------------------

    ! Concentration
    !   clay [g/m^3], pom, caco3, opal [mMol/m^3], agetracer [yr]
    double precision :: SCLA(NXY, SDZ),   SNCLA(NXY, SDZ)
    double precision :: SPOM(NXY, SDZ),   SNPOM(NXY, SDZ)
    double precision :: SCAR(NXY, SDZ),   SNCAR(NXY, SDZ)
    double precision :: SOPA(NXY, SDZ),   SNOPA(NXY, SDZ)
    double precision :: SAGE(NXY, SDZ),   SNAGE(NXY, SDZ)
    double precision :: S13POM(NXY, SDZ), SN13POM(NXY, SDZ)
    double precision :: S14POM(NXY, SDZ), SN14POM(NXY, SDZ)
    double precision :: S13CAR(NXY, SDZ), SN13CAR(NXY, SDZ)
    double precision :: S14CAR(NXY, SDZ), SN14CAR(NXY, SDZ)
    double precision :: SCLAG(NXYDIM, SDZDIM)
    double precision :: SPOMG(NXYDIM, SDZDIM)
    double precision :: SCARG(NXYDIM, SDZDIM)
    double precision :: SOPAG(NXYDIM, SDZDIM)
    double precision :: SAGEG(NXYDIM, SDZDIM)
    double precision :: S13POMG(NXYDIM, SDZDIM)
    double precision :: S14POMG(NXYDIM, SDZDIM)
    double precision :: S13CARG(NXYDIM, SDZDIM)
    double precision :: S14CARG(NXYDIM, SDZDIM)

    ! Weight percent [%]
    double precision :: PCLA(NXY, SDZ)
    double precision :: PPOM(NXY, SDZ)
    double precision :: PCAR(NXY, SDZ)
    double precision :: POPA(NXY, SDZ)
    double precision :: WT(NXY, SDZ)
    double precision :: PCLAG(NXYDIM, SDZDIM)
    double precision :: PPOMG(NXYDIM, SDZDIM)
    double precision :: PCARG(NXYDIM, SDZDIM)
    double precision :: POPAG(NXYDIM, SDZDIM)
    double precision :: WTG(NXYDIM, SDZDIM)

    ! Dissolution rate [1/s], Dissolution term [mMol/m^3/s]
    double precision :: RCLA(NXY, SDZ),   TRCLA(NXY, SDZ)
    double precision :: RPOM(NXY, SDZ),   TRPOM(NXY, SDZ)
    double precision :: RCAR(NXY, SDZ),   TRCAR(NXY, SDZ)
    double precision :: ROPA(NXY, SDZ),   TROPA(NXY, SDZ)
    double precision :: R13POM(NXY, SDZ), TR13POM(NXY, SDZ)
    double precision :: R14POM(NXY, SDZ), TR14POM(NXY, SDZ)
    double precision :: R13CAR(NXY, SDZ), TR13CAR(NXY, SDZ)
    double precision :: R14CAR(NXY, SDZ), TR14CAR(NXY, SDZ)

    ! --------------------------------------------------
    !   Dissolved material
    !     Phosphate, Dissolved Inorganic Carbon,
    !     Alkalinity, Dissolved Oxygen, Silicate
    ! --------------------------------------------------

    ! Concentration [mMol/m^3]
    double precision :: CNIT(NXY, SDZ),   CNNIT(NXY, SDZ)
    double precision :: CDIC(NXY, SDZ),   CNDIC(NXY, SDZ)
    double precision :: CALK(NXY, SDZ),   CNALK(NXY, SDZ)
    double precision :: CDO2(NXY, SDZ),   CNDO2(NXY, SDZ)
    double precision :: CSIL(NXY, SDZ),   CNSIL(NXY, SDZ)
    double precision :: C13DIC(NXY, SDZ), CN13DIC(NXY, SDZ)
    double precision :: C14DIC(NXY, SDZ), CN14DIC(NXY, SDZ)
    double precision :: CDC13(NXY, SDZ),  CDC14(NXY, SDZ)
    double precision :: CCO32(NXY, SDZ)
    double precision :: SATCO32(NXY)
    double precision :: SATSIL(NXY,SDZ)
    double precision :: CNITG(NXYDIM, SDZDIM)
    double precision :: CDICG(NXYDIM, SDZDIM)
    double precision :: CALKG(NXYDIM, SDZDIM)
    double precision :: CDO2G(NXYDIM, SDZDIM)
    double precision :: CSILG(NXYDIM, SDZDIM)
    double precision :: C13DICG(NXYDIM, SDZDIM)
    double precision :: C14DICG(NXYDIM, SDZDIM)
    double precision :: CDC13G(NXYDIM, SDZDIM)
    double precision :: CDC14G(NXYDIM, SDZDIM)
    double precision :: CCO32G(NXYDIM, SDZDIM)
    double precision :: SATCO32G(NXYDIM)
    double precision :: SATSILG(NXYDIM,SDZDIM)

    ! Concentration (bottom water)
    !   temperature [C], salinity [psu], the others [mMol/m^3]
    double precision :: BWTEM(NXY)
    double precision :: BWSAL(NXY)
    double precision :: BWNIT(NXY)
    double precision :: BWDIC(NXY)
    double precision :: BWALK(NXY)
    double precision :: BWDO2(NXY)
    double precision :: BWSIL(NXY)
    double precision :: BW13DIC(NXY)
    double precision :: BW14DIC(NXY)
    double precision :: BWTEMG(NXYDIM)
    double precision :: BWSALG(NXYDIM)
    double precision :: BWNITG(NXYDIM)
    double precision :: BWDICG(NXYDIM)
    double precision :: BWALKG(NXYDIM)
    double precision :: BWDO2G(NXYDIM)
    double precision :: BWSILG(NXYDIM)
    double precision :: BW13DICG(NXYDIM)
    double precision :: BW14DICG(NXYDIM)

    ! Concentration (total)
    !   total phosphorous [mMol/m^3], total silicon [mMol/m^3]
    double precision :: PT(NXY)
    double precision :: SIT(NXY)

    ! Dissolution term [mMol/m^3/s]
    double precision :: TRNIT(NXY, SDZ)
    double precision :: TRDIC(NXY, SDZ)
    double precision :: TRALK(NXY, SDZ)
    double precision :: TRDO2(NXY, SDZ)
    double precision :: TRSIL(NXY, SDZ)
    double precision :: TR13DIC(NXY, SDZ)
    double precision :: TR14DIC(NXY, SDZ)

    ! --------------------------------------------------
    !   Fluxes
    ! --------------------------------------------------

    ! Sedimentation rate
    double precision :: W(NXY, SDZ+1)
    double precision :: WG(NXYDIM, SDZDIM+1)

    ! Sea-bottom Flux
    !   clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
    double precision :: FICLA(NXY),   FOCLA(NXY)
    double precision :: FIPOM(NXY),   FOPOM(NXY)
    double precision :: FICAR(NXY),   FOCAR(NXY)
    double precision :: FIOPA(NXY),   FOOPA(NXY)
    double precision :: FIAGE(NXY)
    double precision :: FI13POM(NXY), FO13POM(NXY)
    double precision :: FI14POM(NXY), FO14POM(NXY)
    double precision :: FI13CAR(NXY), FO13CAR(NXY)
    double precision :: FI14CAR(NXY), FO14CAR(NXY)

    double precision :: FICLAG(NXYDIM), FOCLAG(NXYDIM)
    double precision :: FIPOMG(NXYDIM), FOPOMG(NXYDIM)
    double precision :: FICARG(NXYDIM), FOCARG(NXYDIM)
    double precision :: FIOPAG(NXYDIM), FOOPAG(NXYDIM)
    double precision :: FIAGEG(NXYDIM), FOAGEG(NXYDIM)
    double precision :: FI13POMG(NXYDIM), FO13POMG(NXYDIM)
    double precision :: FI14POMG(NXYDIM), FO14POMG(NXYDIM)
    double precision :: FI13CARG(NXYDIM), FO13CARG(NXYDIM)
    double precision :: FI14CARG(NXYDIM), FO14CARG(NXYDIM)

    ! Sea-bottom Flux [mMol/m^3*cm/s]
    double precision :: FINIT(NXY),   FONIT(NXY)
    double precision :: FIDIC(NXY),   FODIC(NXY)
    double precision :: FIALK(NXY),   FOALK(NXY)
    double precision :: FIDO2(NXY),   FODO2(NXY)
    double precision :: FISIL(NXY),   FOSIL(NXY)
    double precision :: FI13DIC(NXY), FO13DIC(NXY)
    double precision :: FI14DIC(NXY), FO14DIC(NXY)

    ! Sediment-bottom Flux
    !   clay [g/m^3*cm/s], the others [mMol/m^3*cm/s]
    double precision :: BFCLA(NXY)
    double precision :: BFPOM(NXY)
    double precision :: BFCAR(NXY)
    double precision :: BFOPA(NXY)
    double precision :: BF13POM(NXY)
    double precision :: BF14POM(NXY)
    double precision :: BF13CAR(NXY)
    double precision :: BF14CAR(NXY)
    double precision :: BFCLAG(NXYDIM)
    double precision :: BFPOMG(NXYDIM)
    double precision :: BFCARG(NXYDIM)
    double precision :: BFOPAG(NXYDIM)
    double precision :: BF13POMG(NXYDIM)
    double precision :: BF14POMG(NXYDIM)
    double precision :: BF13CARG(NXYDIM)
    double precision :: BF14CARG(NXYDIM)

    ! --------------------------------------------------
    !   Others
    ! --------------------------------------------------

    ! MASK
    ! ocean depth [m]
    double precision :: BDEP(NXY)

    ! Parallel computing
    integer, allocatable :: NRCVC1(:),  NDISP1(:)
    integer, allocatable :: NRCVC2(:),  NDISP2(:)
    integer, allocatable :: NRCVC3(:),  NDISP3(:)
    double precision, allocatable :: DUMMY01(:), DUMMY02(:), DUMMY03(:)
    double precision, allocatable :: DUMMY04(:), DUMMY05(:), DUMMY06(:)
    double precision, allocatable :: DUMMY07(:), DUMMY08(:), DUMMY09(:)
    double precision, allocatable :: DUMMY11(:), DUMMY12(:), DUMMY13(:)
    double precision, allocatable :: DUMMY14(:), DUMMY15(:), DUMMY16(:)
    double precision, allocatable :: DUMMY17(:), DUMMY18(:), DUMMY19(:)
    double precision, allocatable :: DUMMY20(:)
    double precision, allocatable :: DUMMY21(:), DUMMY22(:), DUMMY23(:)
    double precision, allocatable :: DUMMY31(:), DUMMY32(:), DUMMY33(:)
    double precision, allocatable :: DUMMY34(:), DUMMY35(:), DUMMY36(:)
    double precision, allocatable :: DUMMY37(:), DUMMY38(:)
    double precision, allocatable :: DUMMY41(:), DUMMY42(:), DUMMY43(:)
    double precision, allocatable :: DUMMY44(:), DUMMY45(:), DUMMY46(:)
    double precision, allocatable :: DUMMY47(:), DUMMY48(:)
    double precision, allocatable :: DUMMY51(:), DUMMY52(:), DUMMY53(:)
    double precision, allocatable :: DUMMY54(:), DUMMY55(:), DUMMY56(:)
    double precision, allocatable :: DUMMY57(:), DUMMY58(:)

    ! etc
    character(len=100) :: FILENAME, VARNAME, VARUNIT
    double precision :: LON_VALUES(NXDIM)
    double precision :: LAT_VALUES(NYDIM)
    double precision :: DEPTH_VALUES(SDZDIM)
    double precision :: VAR2D_VALUES(NXDIM, NYDIM, 1)
    double precision :: VAR3D_VALUES(NXDIM, NYDIM, SDZDIM, 1)

    integer :: I, J, IJ, K, L, IC
    integer :: TT, YY

!-- SEDIMENT TRACERS : setup --

IF (use_mpi) THEN
    call initialize_digda_mpi
    DO IRANK = 0, NPROCS-1
      call set_parallel_range( &
   &       1, NXYDIM, NPROCS, IRANK, &
   &       JJSTR(IRANK), JJEND(IRANK) &
   &       )
    END DO
END IF

      OPEN(NFPARM, FILE="PARAMET")
      OPEN(NFSTDO, FILE="STDOUT")
      OPEN(NFSOTS, FILE="TSTEP")

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

      call set_digda_clay
      call set_digda_pom
      call set_digda_carbonate
      call set_digda_silicate
      call set_digda_chemeq
      call set_digda_step_solid
      call set_digda_step_liquid

!-- SEDIMENT TRACERS : setup (end) --

!     -- Ocean depth --
IF (use_mpi) THEN
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
      THALF  = THALFYR * 3.1536D+7             !! [yr] --> [s]
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

      DO I = 1, NXYDIM
        FICLAG(I) = BFRCG(I, 1)
        FIPOMG(I) = BFRCG(I, 2)
        FICARG(I) = BFRCG(I, 3)
        FIOPAG(I) = BFRCG(I, 4)
        FIAGEG(I) = BFRCG(I, 5)
        BWTEMG(I) = BFRCG(I, 6)
        BWSALG(I) = BFRCG(I, 7)
        BWNITG(I) = BFRCG(I, 8)
        BWDICG(I) = BFRCG(I, 9)
        BWALKG(I) = BFRCG(I, 10)
        BWSILG(I) = BFRCG(I, 11)
        BWDO2G(I) = BFRCG(I, 12)
        BW13DICG(I) = BFRCG(I, 13)
        BW14DICG(I) = BFRCG(I, 14)
      END DO

IF (use_mpi) THEN
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
IF (use_mpi) THEN
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
IF (use_mpi) THEN
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
IF (use_mpi) THEN
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
         SDFIRST = .false.
      ELSE
IF (use_mpi) THEN
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

      call calc_digda_chemeq( &
     &    BWTEM, BWSAL, CDIC, CALK, &
     &    SIT, PT, &
     &    CCO32 )

      call calc_carbonate_saturation( &
     &    BWTEM, BWSAL, BDEP, &
     &    SATCO32 )

!     -- Dissolution term --

      call calc_clay_dissolution( &
     &    SCLA, &
     &    RCLA )

      call calc_pom_dissolution( &
     &    SPOM, CDO2, CP, OP, PHI, &
     &    RPOM )

      call calc_carbonate_dissolution( &
     &    SCAR, CCO32, SATCO32, PHI, &
     &    RCAR )

      call calc_silicate_dissolution( &
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

      call calc_step_solid_clay( &
     &    FICLA, TRCLA, W, PHI, PHII, DBI, SCLA, RHOCLA, &
     &    SNCLA, BFCLA )

      call calc_step_solid( &
     &    FIPOM, TRPOM, W, PHI, PHII, DBI, SPOM, &
     &    SNPOM, BFPOM )

      call calc_step_solid( &
     &    FICAR, TRCAR, W, PHI, PHII, DBI, SCAR, &
     &    SNCAR, BFCAR )

      call calc_step_solid( &
     &    FIOPA, TROPA, W, PHI, PHII, DBI, SOPA, &
     &    SNOPA, BFOPA )

      call calc_step_solid_age( &
     &    FIAGE, W, PHI, PHII, DBI, SAGE, &
     &    SNAGE )

      IF ( CALC13 ) THEN
      call calc_step_solid( &
     &    FI13POM, TR13POM, W, PHI, PHII, DBI, S13POM, &
     &    SN13POM, BF13POM )

      call calc_step_solid( &
     &    FI13CAR, TR13CAR, W, PHI, PHII, DBI, S13CAR, &
     &    SN13CAR, BF13CAR )
      END IF

      IF ( CALC14 ) THEN
      call calc_step_solid( &
     &    FI14POM, TR14POM, W, PHI, PHII, DBI, S14POM, &
     &    SN14POM, BF14POM )

      call calc_step_solid( &
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

      call calc_step_liquid( &
     &    BWNIT, TRNIT, PHI, PHII, DEINIT, CNIT, &
     &    CNNIT )

      call calc_step_liquid( &
     &    BWDIC, TRDIC, PHI, PHII, DEIDIC, CDIC, &
     &    CNDIC )

      call calc_step_liquid( &
     &    BWALK, TRALK, PHI, PHII, DEIALK, CALK, &
     &    CNALK )

      call calc_step_liquid_o2( &
     &    BWDO2, TRDO2, PHI, PHII, DEIDO2, CDO2, &
     &    CNDO2 )

      call calc_step_liquid( &
     &    BWSIL, TRSIL, PHI, PHII, DEISIL, CSIL, &
     &    CNSIL )

      IF ( CALC13 ) THEN
      call calc_step_liquid( &
     &    BW13DIC, TR13DIC, PHI, PHII, DEI13DIC, C13DIC, &
     &    CN13DIC )
      END IF

      IF ( CALC14 ) THEN
      call calc_step_liquid( &
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
IF (use_mpi) THEN
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
!     IF ( WOHIST ) THEN
        IF ( TT == TEND .OR. MOD(TT, TOUT) == 0 ) THEN
!     ELSE
!       IF ( TT == TEND ) THEN
!     END IF

IF (use_mpi) THEN

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

      SDOFIRST = .false.
      END IF

      call MPI_GATHERV( SCLA,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY01, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SPOM,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY02, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SCAR,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY03, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SOPA,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY04, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SAGE,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY05, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CNIT,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY11, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CDIC,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY12, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CALK,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY13, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CDO2,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY14, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CSIL,    ISIZE2,         MPI_REAL8, &
     &                  DUMMY15, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CCO32,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY16, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SATSIL,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY21, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( SATCO32, ISIZE1,         MPI_REAL8, &
     &                  DUMMY22, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( W,       ISIZE3,         MPI_REAL8, &
     &                  DUMMY23, NRCVC3, NDISP3, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FICLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY31, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FIPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY32, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FICAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY33, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FIOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY34, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FOCLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY41, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FOPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY42, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FOCAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY43, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FOOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY44, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BFCLA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY51, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BFPOM,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY52, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BFCAR,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY53, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BFOPA,   ISIZE1,         MPI_REAL8, &
     &                  DUMMY54, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )

      IF ( CALC13 ) THEN
      call MPI_GATHERV( S13POM,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY06, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( S13CAR,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY08, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( C13DIC,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY17, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CDC13,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY19, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FI13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY35, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FI13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY37, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FO13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY45, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FO13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY47, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BF13POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY55, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BF13CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY57, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      END IF

      IF ( CALC14 ) THEN
      call MPI_GATHERV( S14POM,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY07, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( S14CAR,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY09, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( C14DIC,  ISIZE2,         MPI_REAL8, &
     &                  DUMMY18, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( CDC14,   ISIZE2,         MPI_REAL8, &
     &                  DUMMY20, NRCVC2, NDISP2, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FI14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY36, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FO14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY46, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( FO14CAR, ISIZE1,         MPI_REAL8, &
     &                  DUMMY48, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BF14POM, ISIZE1,         MPI_REAL8, &
     &                  DUMMY56, NRCVC1, NDISP1, MPI_REAL8, &
     &                  IROOT,   MPI_COMM_WORLD, IERR     )
      call MPI_GATHERV( BF14CAR, ISIZE1,         MPI_REAL8, &
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
              !FICLAG(IJ)   = DUMMY31(IJ)
              !FIPOMG(IJ)   = DUMMY32(IJ)
              !FICARG(IJ)   = DUMMY33(IJ)
              !FIOPAG(IJ)   = DUMMY34(IJ)
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
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)

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
          ! FICLAG(IJ)     = FICLA(IJ)
          ! FIPOMG(IJ)     = FIPOM(IJ)
          ! FICARG(IJ)     = FICAR(IJ)
          ! FIOPAG(IJ)     = FIOPA(IJ)
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

IF (use_mpi) THEN
      IF ( MYRANK .EQ. IROOT ) THEN

         DO I = 1, NXDIM
         LON_VALUES(I) = ((I - 1) * 360.0 / NXDIM) + (180.0 / NXDIM)
         END DO

         DO I = 1, NYDIM
         LAT_VALUES(I) = -90.0 + ((I - 1) * 180.0 / NYDIM) + (90.0 / NYDIM)
         END DO

!        FILENAME = "scla.nc"
!        VARNAME = "Clay"
!        VARUNIT = "g m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = SCLAG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "spom.nc"
!        VARNAME = "OrgM"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = SPOMG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "scar.nc"
!        VARNAME = "CaCO3"
!        VARUNIT = "g m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = SCARG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "sopa.nc"
!        VARNAME = "SiO2"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = SOPAG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "sage.nc"
!        VARNAME = "Age tracer"
!        VARUNIT = "year"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = SAGEG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "cnit.nc"
!        VARNAME = "Phosphate"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CNITG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "cdic.nc"
!        VARNAME = "DIC"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CDICG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "calk.nc"
!        VARNAME = "Alkalinity"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CALKG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "csil.nc"
!        VARNAME = "Silicate"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CSILG(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "cdo2.nc"
!        VARNAME = "Oxygen"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CDO2G(IJ, K)
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

         FILENAME = "ficla.nc"
         VARNAME = "Clay_bottom_flux"
         VARUNIT = "g m-3 * cm s-1"
         DO J = 1, NYDIM
           DO I = 1, NXDIM
             IJ = I + (J - 1) * NXDIM
             VAR2D_VALUES(I, J, 1) = FICLAG(IJ)
           END DO
         END DO
         call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

         FILENAME = "fipom.nc"
         VARNAME = "OrgM_bottom_flux"
         VARUNIT = "mMol m-3 * cm s-1"
         DO J = 1, NYDIM
           DO I = 1, NXDIM
             IJ = I + (J - 1) * NXDIM
             VAR2D_VALUES(I, J, 1) = FIPOMG(IJ)
           END DO
         END DO
         call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

         FILENAME = "ficar.nc"
         VARNAME = "CaCO3_bottom_flux"
         VARUNIT = "mMol m-3 * cm s-1"
         DO J = 1, NYDIM
           DO I = 1, NXDIM
             IJ = I + (J - 1) * NXDIM
             VAR2D_VALUES(I, J, 1) = FICARG(IJ)
           END DO
         END DO
         call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

         FILENAME = "fiopa.nc"
         VARNAME = "SiO2_bottom_flux"
         VARUNIT = "mMol m-3 * cm s-1"
         DO J = 1, NYDIM
           DO I = 1, NXDIM
             IJ = I + (J - 1) * NXDIM
             VAR2D_VALUES(I, J, 1) = FIOPAG(IJ)
           END DO
         END DO
         call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "focla.nc"
!        VARNAME = "Clay_diffusion"
!        VARUNIT = "g m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = FOCLAG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "fopom.nc"
!        VARNAME = "OrgM_diffusion"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = FOPOMG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "focar.nc"
!        VARNAME = "CaCO3_diffusion"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = FOCARG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "foopa.nc"
!        VARNAME = "SiO2_diffusion"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = FOOPAG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "bfcla.nc"
!        VARNAME = "Clay_burial"
!        VARUNIT = "g m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = BFCLAG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "bfpom.nc"
!        VARNAME = "OrgM_burial"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = BFPOMG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "bfcar.nc"
!        VARNAME = "CaCO3_burial"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = BFCARG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "bfopa.nc"
!        VARNAME = "SiO2_burial"
!        VARUNIT = "mMol m-3 * cm s-1"
!        DO J = 1, NYDIM
!          DO I = 1, NXDIM
!            IJ = I + (J - 1) * NXDIM
!            VAR2D_VALUES(I, J, 1) = BFOPAG(IJ)
!          END DO
!        END DO
!        call WRITE_NETCDF_2D(FILENAME, VARNAME, VARUNIT, VAR2D_VALUES, LON_VALUES, LAT_VALUES, YY)

!        FILENAME = "w.nc"
!        VARNAME = "Sediment_rate"
!        VARUNIT = "mmol m-3"
!        DO K = 1, SDZDIM+1
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = WG(IJ, K)
!            END DO
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ+1, YY)

!        FILENAME = "cdc13.nc"
!        VARNAME = "delta_carbon-13"
!        VARUNIT = "permil"
!        DO K = 1, SDZDIM
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CDC13G(IJ, K)
!            END DO
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

!        FILENAME = "cdc14.nc"
!        VARNAME = "Delta_carbon-14"
!        VARUNIT = "permil"
!        DO K = 1, SDZDIM+1
!          DO J = 1, NYDIM
!            DO I = 1, NXDIM
!              IJ = I + (J - 1) * NXDIM
!              VAR3D_VALUES(I, J, K, 1) = CDC14G(IJ, K)
!            END DO
!          END DO
!        END DO
!        call WRITE_NETCDF_3D(FILENAME, VARNAME, VARUNIT, VAR3D_VALUES, LON_VALUES, LAT_VALUES, LZ, YY)

         CFILE  = 'scla'
         CTITLE = 'Clay'
         CUNIT  = 'g/m^3'
         call write_gtoolfile(      YY,   SCLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'spom'
         CTITLE = 'OrgM'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   SPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'scar'
         CTITLE = 'CaCO3'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   SCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'sopa'
         CTITLE = 'SiO2'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   SOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'sage'
         CTITLE = 'Age tracer'
         CUNIT  = 'year'
         call write_gtoolfile(      YY,   SAGEG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cnit'
         CTITLE = 'Phosphate'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   CNITG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdic'
         CTITLE = 'DIC'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   CDICG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'calk'
         CTITLE = 'Alkalinity'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   CALKG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'csil'
         CTITLE = 'Silicate'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   CSILG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdo2'
         CTITLE = 'Oxygen'
         CUNIT  = 'mmol/m^3'
         call write_gtoolfile(      YY,   CDO2G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'ficla'
         CTITLE = 'Clay bottom flux'
         CUNIT  = 'g/m^3*cm/s'
         call write_gtoolfile(      YY,  FICLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fipom'
         CTITLE = 'OrgM bottom flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FIPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'ficar'
         CTITLE = 'CaCO3 bottom flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FICARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fiopa'
         CTITLE = 'SiO2 bottom flux'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FIOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'fopom'
         CTITLE = 'OrgM diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FOPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'focar'
         CTITLE = 'CaCO3 diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FOCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'foopa'
         CTITLE = 'SiO2 diffusion'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  FOOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfcla'
         CTITLE = 'Clay burial'
         CUNIT  = 'g/m^3*cm/s'
         call write_gtoolfile(      YY,  BFCLAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'bfpom'
         CTITLE = 'OrgM burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  BFPOMG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfcar'
         CTITLE = 'CaCO3 burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  BFCARG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'bfopa'
         CTITLE = 'SiO2 burial'
         CUNIT  = 'mMol/m^3*cm/s'
         call write_gtoolfile(      YY,  BFOPAG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,       1  )
         CFILE  = 'w'
         CTITLE = 'Sediment rate'
         CUNIT  = 'cm/s'
         call write_gtoolfile(      YY,      WG,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM, SDZDIM+1 )
         CFILE  = 'cdc13'
         CTITLE = 'delta carbon-13'
         CUNIT  = 'permil'
         call write_gtoolfile(      YY,  CDC13G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )
         CFILE  = 'cdc14'
         CTITLE = 'Delta carbon-14'
         CUNIT  = 'permil'
         call write_gtoolfile(      YY,  CDC14G,   CFILE,  CTITLE,   CUNIT, &
     &                   NXDIM,   NYDIM,  SDZDIM  )

         DO K = SKST, SKEN
            DO IJ = 1, NXYDIM
               SDTX(IJ, K, 1) = SCLAG(IJ, K)
               SDTX(IJ, K, 2) = SPOMG(IJ, K)
               SDTX(IJ, K, 3) = SCARG(IJ, K)
               SDTX(IJ, K, 4) = SOPAG(IJ, K)
               SDTX(IJ, K, 5) = SAGEG(IJ, K)
               SDTX(IJ, K, 6) = CNITG(IJ, K)
               SDTX(IJ, K, 7) = CDICG(IJ, K)
               SDTX(IJ, K, 8) = CALKG(IJ, K)
               SDTX(IJ, K, 9) = CSILG(IJ, K)
               SDTX(IJ, K, 10) = CDO2G(IJ, K)
            IF ( CALC13 ) THEN
               SDTX(IJ, K, 11) = S13POMG(IJ, K)
               SDTX(IJ, K, 12) = S13CARG(IJ, K)
               SDTX(IJ, K, 13) = C13DICG(IJ, K)
            END IF
            IF ( CALC14 ) THEN
               SDTX(IJ, K, 14) = S14POMG(IJ, K)
               SDTX(IJ, K, 15) = S14CARG(IJ, K)
               SDTX(IJ, K, 16) = C14DICG(IJ, K)
            END IF
            END DO
         END DO

      END IF

      call MPI_BCAST(SDTX, NXYDIM*SDZDIM*SEDTDM, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)

END IF

      END IF ! -- End of Output

      END DO ! -- End of Time integration

      end subroutine digda_main

end module module_digda_main
