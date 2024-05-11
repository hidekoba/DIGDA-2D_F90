!######################################################################
module sddim_module
      implicit none

#include "parallel_config.h"
#include "grid_config.h"
#include "time_step_config.h"

#ifdef PARALLEL
      INTEGER   IERR, NPROCS, MYRANK
      INTEGER   IREQ, ISTATUS
      INTEGER   NCPU
      INTEGER   IRANK, ISIZE1, ISIZE2, ISIZE3, IIST, IIEN, IROOT, IJKG
      INTEGER   ISEND1, ISEND2, ISEND3, ISEND4, ISEND5
      INTEGER   IRECV1, IRECV2, IRECV3, IRECV4, IRECV5
      INTEGER   IREQ1,  IREQ2,  IREQ3,  IREQ4,  IREQ5
      INTEGER   IREQ11, IREQ12, IREQ13, IREQ14, IREQ15
      INTEGER   IREQ21, IREQ22, IREQ23, IREQ24, IREQ25
      PARAMETER (IROOT  =      0)
!     PARAMETER (NCPU   =      4)
!     PARAMETER (NCPU   =      8)
      PARAMETER (NCPU   =     16)
!     PARAMETER (NCPU   =     18)
      INTEGER   JJSTR(0:NCPU-1), JJEND(0:NCPU-1)
#else
      INTEGER   NCPU
      PARAMETER (NCPU   =      1)
#endif

      INTEGER   SEDTDM
!     PARAMETER (SEDTDM =     10)
      PARAMETER (SEDTDM =     16)
      INTEGER   SDZDIM, SZTDIM
      INTEGER   SKST, SKEN
#ifdef ZGRID_CONSTANT
      PARAMETER (SDZDIM =     20)
      PARAMETER (SKST   =      1, SKEN   =     20)
      PARAMETER (SZTDIM = SDZDIM*SEDTDM)
#endif
#ifdef ZGRID_ARCHER
      PARAMETER (SDZDIM =     16)
      PARAMETER (SKST   =      1, SKEN   =     16)
      PARAMETER (SZTDIM = SDZDIM*SEDTDM)
#endif
#ifdef ZGRID_HEINZE
      PARAMETER (SDZDIM =     10)
      PARAMETER (SKST   =      1, SKEN   =     10)
      PARAMETER (SZTDIM = SDZDIM*SEDTDM)
#endif
      REAL*8 LZ(SKST:SKEN),     LZI(SKST:SKEN+1)
      REAL*8 DLZI(SKST:SKEN+1), DLZ(SKST:SKEN)
#ifdef ZGRID_CONSTANT
!     Constant
      DATA   LZ   / 2.500D-3, 7.500D-3, 1.250D-2, 1.750D-2, 2.250D-2, &
     &              2.750D-2, 3.250D-2, 3.750D-2, 4.250D-2, 4.750D-2, &
     &              5.250D-2, 5.750D-2, 6.250D-2, 6.750D-2, 7.250D-2, &
     &              7.750D-2, 8.250D-2, 8.750D-2, 9.250D-2, 9.750D-2 /
      DATA   LZI  / 0.000D0,  5.000D-3, 1.000D-2, 1.500D-2, 2.000D-2, &
     &              2.500D-2, 3.000D-2, 3.500D-2, 4.000D-2, 4.500D-2, &
     &              5.000D-2, 5.500D-2, 6.000D-2, 6.500D-2, 7.000D-2, &
     &              7.500D-2, 8.000D-2, 8.500D-2, 9.000D-2, 9.500D-2, &
     &              1.000D-1                                         /
      DATA   DLZ  / 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3 /
      DATA   DLZI / 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              5.000D-3                                         /
#endif
#ifdef ZGRID_ARCHER
!     Archer et al. (2002)
      DATA   LZ   / 2.500D-4, 1.750D-3, 5.000D-3, 1.100D-2, 2.000D-2, &
     &              3.250D-2, 5.000D-2, 7.250D-2, 1.025D-1, 1.400D-1, &
     &              1.850D-1, 2.300D-1, 2.750D-1, 4.000D-1, 6.250D-1, &
     &              8.750D-1 /
      DATA   LZI  / 0.000D0,  5.000D-4, 3.000D-3, 7.000D-3, 1.500D-2, &
     &              2.500D-2, 4.000D-2, 6.000D-2, 8.500D-2, 1.200D-1, &
     &              1.600D-1, 2.100D-1, 2.500D-1, 3.000D-1, 5.000D-1, &
     &              7.500D-1, 1.000D0 /
      DATA   DLZ  / 5.000D-4, 2.500D-3, 4.000D-3, 8.000D-3, 1.000D-2, &
     &              1.500D-2, 2.000D-2, 2.500D-2, 3.500D-2, 4.000D-2, &
     &              5.000D-2, 4.000D-2, 5.000D-2, 2.000D-1, 2.500D-1, &
     &              2.500D-1 /
      DATA   DLZI / 5.000D-4, 1.500D-3, 3.250D-3, 6.000D-3, 9.000D-3, &
     &              1.250D-2, 1.750D-2, 2.250D-2, 3.000D-2, 3.750D-2, &
     &              4.500D-2, 4.500D-2, 4.500D-2, 1.250D-1, 2.250D-1, &
     &              2.500D-1, 2.500D-1 /
#endif &
#ifdef ZGRID_HEINZE &
!     Heinze et al. (1999) &
      DATA   LZ   / 1.500D-3, 4.500D-3, 8.500D-3, 1.350D-2, 1.850D-2, &
     &              2.600D-2, 3.600D-2, 4.600D-2, 6.325D-2, 8.775D-2 /
      DATA   LZI  / 0.000D0,  3.000D-3, 6.000D-3, 1.100D-2, 1.600D-2, &
     &              2.100D-2, 3.100D-2, 4.100D-2, 5.100D-2, 7.550D-2, &
     &              1.000D-1                                         /
      DATA   DLZ  / 3.000D-3, 3.000D-3, 5.000D-3, 5.000D-3, 5.000D-3, &
     &              1.000D-2, 1.000D-2, 1.000D-2, 2.450D-2, 2.450D-2 /
      DATA   DLZI / 3.000D-3, 3.000D-3, 4.000D-3, 5.000D-3, 5.000D-3, &
     &              7.500D-3, 1.000D-2, 1.000D-2, 1.725D-2, 2.450D-2, &
     &              2.450D-2                                         /
#endif

!     -- zcdim.F --
      INTEGER   NXDIM, NYDIM, NZDIM
      INTEGER   NXYDIM
#ifdef ONE_DIM
      PARAMETER (NXDIM  =      1)
      PARAMETER (NYDIM  =      1)
      PARAMETER (NZDIM  =     44)
      PARAMETER (NXYDIM  = NXDIM*NYDIM)
#endif
#ifdef Grid_128x64
      PARAMETER (NXDIM  =    128)
      PARAMETER (NYDIM  =     64)
      PARAMETER (NZDIM  =     44)
      PARAMETER (NXYDIM  = NXDIM*NYDIM)
#endif
#ifdef Grid_256x192
      PARAMETER (NXDIM  =    256)
      PARAMETER (NYDIM  =    192)
      PARAMETER (NZDIM  =     44)
      PARAMETER (NXYDIM  = NXDIM*NYDIM)
#endif
#ifdef Grid_360x180
      PARAMETER (NXDIM  =    360)
      PARAMETER (NYDIM  =    180)
      PARAMETER (NZDIM  =     44)
      PARAMETER (NXYDIM  = NXDIM*NYDIM)
#endif
#ifdef Grid_2160x1080
      PARAMETER (NXDIM  =   2160)
      PARAMETER (NYDIM  =   1080)
      PARAMETER (NZDIM  =     44)
      PARAMETER (NXYDIM  = NXDIM*NYDIM)
#endif

      INTEGER   NX, NY, NZ, SDZ
      INTEGER   NXY, NXYZ, NXYZDM
      PARAMETER (NX     = NXDIM)
      PARAMETER (NY     = NYDIM/NCPU)
      PARAMETER (NZ     = NZDIM)
      PARAMETER (SDZ    = SDZDIM)
      PARAMETER (NXY    = NX*NY)
      PARAMETER (NXYZ   = NXY*SDZ)
      PARAMETER (NXYZDM = NXYDIM*SDZ)

      INTEGER ISTR, IEND
      PARAMETER(ISTR   =      1, IEND   =     NX)
      INTEGER JSTR, JEND
      PARAMETER(JSTR   =      1, JEND   =     NY)
      INTEGER KSTR, KEND
      PARAMETER(KSTR    =     1, KEND   =     NZ)
      INTEGER   IJTSTR, IJTEND
      PARAMETER (IJTSTR =     1, IJTEND =    NXY)

      REAL*8    TSYR, TS
#ifdef Time_Step_1y
      PARAMETER (TSYR  = 365.0D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_1_2y
      PARAMETER (TSYR  = 182.5D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_50d
      PARAMETER (TSYR  =  50.0D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_30d
      PARAMETER (TSYR  =  1.0D0/12.0D0)   !! [yr]
!     PARAMETER (TSYR  =  30.0D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_10d
      PARAMETER (TSYR  =  10.0D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_5d
      PARAMETER (TSYR  =  5.00D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_2d
      PARAMETER (TSYR  =  2.00D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_1d
      PARAMETER (TSYR  =  1.00D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_1_10d
      PARAMETER (TSYR  =  0.10D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_1_25d
      PARAMETER (TSYR  =  0.04D0/365D0)   !! [yr]
#endif
#ifdef Time_Step_1_50d
      PARAMETER (TSYR  =  0.02D0/365D0)   !! [yr]
#endif
      PARAMETER (TS = TSYR * 3.1536D+7)   !! [s]
      INTEGER   YY, YEND

      REAL*8 DEPTH(NZDIM)
      DATA   DEPTH /    5.0D0,   10.0D0,   15.0D0,   20.0D0,   25.0D0, &
     &                 30.0D0,   37.5D0,   45.0D0,   55.0D0,   65.0D0, &
     &                 75.0D0,   85.0D0,  100.0D0,  120.0D0,  145.0D0, &
     &                170.0D0,  200.0D0,  250.0D0,  325.0D0,  400.0D0, &
     &                500.0D0,  600.0D0,  750.0D0,  900.0D0, 1050.0D0, &
     &               1200.0D0, 1400.0D0, 1650.0D0, 1900.0D0, 2150.0D0, &
     &               2400.0D0, 2650.0D0, 2900.0D0, 3150.0D0, 3400.0D0, &
     &               3650.0D0, 3900.0D0, 4150.0D0, 4400.0D0, 4650.0D0, &
     &               4900.0D0, 5150.0D0, 5400.0D0, 5500.0D0           /

!     -- zcmask.F
!     REAL*8    AMSKT(NXYDIM, NZDIM)
      INTEGER   NBOT(NXYDIM)

!     -- zcgrid.F
      REAL*8    DZ(NXYDIM, NZDIM)
!     REAL*8    TS

end module sddim_module
