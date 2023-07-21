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

subroutine CHO_REOINI()
!
! Purpose: initializations for vector reordering.

use ChoReo

implicit real*8(a-h,o-z)
#include "cholesky.fh"
#include "choorb.fh"
integer :: I, J, MULD2H
! Statement function
MULD2H(I,J) = ieor(I-1,J-1)+1

call IZERO(NNBST,NSYM)
do ISYMA=1,NSYM
  do ISYMB=1,ISYMA-1
    NABPK(ISYMA,ISYMB) = NBAS(ISYMA)*NBAS(ISYMB)
    NABPK(ISYMB,ISYMA) = NABPK(ISYMA,ISYMB)
    ISYM = MULD2H(ISYMB,ISYMA)
    NNBST(ISYM) = NNBST(ISYM)+NABPK(ISYMA,ISYMB)
  end do
  NABPK(ISYMA,ISYMA) = NBAS(ISYMA)*(NBAS(ISYMA)+1)/2
  NNBST(1) = NNBST(1)+NABPK(ISYMA,ISYMA)
end do

call CHO_OPFVEC(1,0)

end subroutine CHO_REOINI
