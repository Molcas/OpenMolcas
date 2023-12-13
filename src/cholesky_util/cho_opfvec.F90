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

subroutine CHO_OPFVEC(ISYM,IOPT)
!
! Purpose: open/close files for full storage vectors, sym. ISYM.

use Symmetry_Info, only: Mul
use Cholesky, only: LUFV, nSym, REONAM
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ISYM, IOPT
integer(kind=iwp) :: ISYMA, ISYMB, LUNIT
character(len=6) :: FNAME
character(len=*), parameter :: SECNAM = 'CHO_OPFVEC'

if (IOPT == 0) then
  LUFV(1:NSYM,1:NSYM) = -1
else if (IOPT == 1) then
  do ISYMB=1,NSYM
    ISYMA = MUL(ISYMB,ISYM)
    if (ISYMA >= ISYMB) then
      write(FNAME,'(A4,I1,I1)') REONAM,ISYMA,ISYMB
      LUNIT = 7
      call DANAME_MF_WA(LUNIT,FNAME)
      LUFV(ISYMA,ISYMB) = LUNIT
      LUFV(ISYMB,ISYMA) = LUNIT
    end if
  end do
else if (IOPT == 2) then
  do ISYMB=1,NSYM
    ISYMA = MUL(ISYMB,ISYM)
    if (ISYMA >= ISYMB) then
      LUNIT = LUFV(ISYMA,ISYMB)
      call DACLOS(LUNIT)
      LUFV(ISYMA,ISYMB) = -1
      LUFV(ISYMB,ISYMA) = -1
    end if
  end do
else
  call CHO_QUIT('IOPT error in '//SECNAM,104)
end if

end subroutine CHO_OPFVEC
