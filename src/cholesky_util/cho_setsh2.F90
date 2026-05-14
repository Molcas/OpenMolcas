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

subroutine CHO_SETSH2(ISHLSO,ISOSHL,NBSTSH,NBAST,NSHELL)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NBAST, ISOSHL(NBAST), NSHELL, NBSTSH(NSHELL)
integer(kind=iwp), intent(out) :: ISHLSO(NBAST)
integer(kind=iwp) :: I, ICOUNT, ISHL

do ISHL=1,NSHELL
  ICOUNT = 0
  I = 0
  do while ((I < NBAST) .and. (ICOUNT < NBSTSH(ISHL)))
    I = I+1
    if (ISOSHL(I) == ISHL) then
      ICOUNT = ICOUNT+1
      ISHLSO(I) = ICOUNT
    end if
  end do
end do

end subroutine CHO_SETSH2
