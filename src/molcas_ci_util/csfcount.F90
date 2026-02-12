!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CSFCOUNT(CIS,NSYM,NUW)

use Symmetry_Info, only: Mul
use gugx, only: CIStruct
use Definitions, only: iwp

implicit none
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp), intent(in) :: NSYM
integer(kind=iwp), intent(out) :: NUW
integer(kind=iwp) :: ISYDWN, ISYM, ISYTOT, ISYUP, MV, N

! CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
! SEPARATED FOR EACH MIDVERTEX AND SYMMETRY

NUW = 0
do MV=1,CIS%nMidV
  do ISYM=1,NSYM
    CIS%IOW(1,ISYM,MV) = NUW*CIS%nIpWlk
    NUW = NUW+CIS%NOW(1,ISYM,MV)
  end do
end do
CIS%nWalk = NUW
do MV=1,CIS%nMidV
  do ISYM=1,NSYM
    CIS%IOW(2,ISYM,MV) = CIS%nWalk*CIS%nIpWlk
    CIS%nWalk = CIS%nWalk+CIS%NOW(2,ISYM,MV)
  end do
end do

! CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
! SEPARATED BY MIDVERTICES AND SYMMETRY.
! FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.

CIS%NCSF(:) = 0
do ISYTOT=1,NSYM
  do MV=1,CIS%nMidV
    do ISYUP=1,NSYM
      ISYDWN = Mul(ISYTOT,ISYUP)
      N = CIS%NOW(1,ISYUP,MV)*CIS%NOW(2,ISYDWN,MV)
      CIS%NOCSF(ISYUP,MV,ISYTOT) = N
      CIS%IOCSF(ISYUP,MV,ISYTOT) = CIS%NCSF(ISYTOT)
      CIS%NCSF(ISYTOT) = CIS%NCSF(ISYTOT)+N
#     include "compiler_features.h"
#     ifdef _BUGGY_INTEL_LLVM_
      ! dummy statement to work around compiler bug, will never be executed
      if (ISYUP > 99) CIS%NCSF(ISYTOT) = -1
#     endif
    end do
  end do
end do

end subroutine CSFCOUNT
