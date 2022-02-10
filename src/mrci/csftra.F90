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

implicit real*8(A-H,O-Z)
character*4 KEY
dimension CI(NCONF), AREF(NREF,NREF)
#include "SysDef.fh"
#include "mrci.fh"
dimension TMP(MXREF)

if (NREF == 1) return
if (KEY == ' CSF') then
  do I=1,NREF
    SUM = 0.0d00
    do J=1,NREF
      SUM = SUM+AREF(I,J)*CI(IREFX(J))
    end do
    TMP(I) = SUM
  end do
else
  do I=1,NREF
    SUM = 0.0d00
    do J=1,NREF
      SUM = SUM+AREF(J,I)*CI(IREFX(J))
    end do
    TMP(I) = SUM
  end do
end if
do I=1,NREF
  CI(IREFX(I)) = TMP(I)
end do

return

end subroutine CSFTRA
