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

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Cholesky, only: NABPK, NBAS, NNBST, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISYM, ISYMA, ISYMB

NNBST(1:NSYM) = 0
do ISYMA=1,NSYM
  do ISYMB=1,ISYMA-1
    NABPK(ISYMA,ISYMB) = NBAS(ISYMA)*NBAS(ISYMB)
    NABPK(ISYMB,ISYMA) = NABPK(ISYMA,ISYMB)
    ISYM = MUL(ISYMB,ISYMA)
    NNBST(ISYM) = NNBST(ISYM)+NABPK(ISYMA,ISYMB)
  end do
  NABPK(ISYMA,ISYMA) = nTri_Elem(NBAS(ISYMA))
  NNBST(1) = NNBST(1)+NABPK(ISYMA,ISYMA)
end do

call CHO_OPFVEC(1,0)

end subroutine CHO_REOINI
