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

subroutine CSFTRA(KEY,CI,AREF)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "mrci.fh"
character(len=4) :: KEY
real(kind=wp) :: CI(NCONF), AREF(NREF,NREF)
integer(kind=iwp) :: I, J
real(kind=wp) :: RSUM, TMP(MXREF) !IFG

if (NREF == 1) return
if (KEY == ' CSF') then
  do I=1,NREF
    RSUM = Zero
    do J=1,NREF
      RSUM = RSUM+AREF(I,J)*CI(IREFX(J))
    end do
    TMP(I) = RSUM
  end do
else
  do I=1,NREF
    RSUM = Zero
    do J=1,NREF
      RSUM = RSUM+AREF(J,I)*CI(IREFX(J))
    end do
    TMP(I) = RSUM
  end do
end if
do I=1,NREF
  CI(IREFX(I)) = TMP(I)
end do

return

end subroutine CSFTRA
