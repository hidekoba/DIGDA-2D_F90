!######################################################################
module zcfile_module
      implicit none

      INTEGER   NCF
      PARAMETER(NCF = 256)
      INTEGER   NFMIN, NFMAX
      PARAMETER(NFMIN = 21, NFMAX = 89)

      INTEGER   NFSTDI, NFSTDO, NFPARM
      PARAMETER(NFSTDI=5, NFSTDO = 6, NFPARM=7)

end module zcfile_module
