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

subroutine faibj5(LENBUF,JTURN,IBUF,BUF,AIBJ,ABIJ)

use Definitions, only: wp, iwp

implicit none
#include "mrci.fh"
integer(kind=iwp) :: LENBUF, JTURN, IBUF(NBITM3+2)
real(kind=wp) :: BUF(NBITM3), ABIJ(NVSQ), AIBJ(NVSQ)
integer(kind=iwp) :: i

if (LENBUF > 0) then
  if (JTURN == 1) then
    do i=1,LENBUF
      aibj(IBUF(i)) = buf(i)
    end do
  else
    do i=1,LENBUF
      abij(IBUF(i)) = buf(i)
    end do
  end if
end if

return

end subroutine faibj5
