      SUBROUTINE MKNSM()
!     PUPROSE: CREATE THE SYMMETRY INDEX VECTOR
!
      use gugx, only: SGS
      use stdalloc, only: mma_allocate
      IMPLICIT None
!
! to get some dimensions
#include "rasdim.fh"
! NSM from rasscf,fh
#include "rasscf.fh"
! NSYM from general.fh
#include "general.fh"
! NGAS and NGSSH from gas.fh
#include "gas.fh"
!
      Integer IGAS, ISYM, LEV, NLEV, NSTA

      NLEV=0
      DO IGAS=1,NGAS
        DO ISYM=1,NSYM
          NSTA=NLEV+1
          NLEV=NLEV+NGSSH(IGAS,ISYM)
          DO LEV=NSTA,NLEV
            NSM(LEV)=ISYM
          END DO
        END DO
      END DO

      If (SGS%nSym/=0) Then
         SGS%nLev=nLev
         Call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')
         SGS%ISM(1:nLev)=NSM(1:nLev)
      End If

      END SUBROUTINE MKNSM

