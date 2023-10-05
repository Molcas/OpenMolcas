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

subroutine CHO_GETSTOR(VCSTOR)
!
! Purpose: get total vector storage (in words).

use Cholesky, only: LuPri, MaxVec, nSym, NumCho
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: VCSTOR(nSym)
integer(kind=iwp) :: ISYM
character(len=*), parameter :: SECNAM = 'CHO_GETSTOR'

do ISYM=1,NSYM
  if (NUMCHO(ISYM) > MAXVEC) then
    write(LUPRI,*) SECNAM,': too many Cholesky vectors in symmetry ',ISYM,': ',NUMCHO(ISYM)
    call CHO_QUIT('Error in '//SECNAM,103)
    VCSTOR(ISYM) = Zero
  else if (NUMCHO(ISYM) < 0) then
    write(LUPRI,*) SECNAM,': negative #Cholesky vectors in symmetry ',ISYM,': ',NUMCHO(ISYM)
    call CHO_QUIT('Error in '//SECNAM,103)
    VCSTOR(ISYM) = Zero
  else
    call CHO_GETSTOR_S(VCSTOR(ISYM),ISYM)
  end if
end do

end subroutine CHO_GETSTOR
