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

subroutine Cho_PrtMaxMem(Location)
!
! Purpose: print max. available memory block.

use Cholesky, only: LuPri
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Location
integer(kind=iwp) :: l, lMax
real(kind=wp) :: dlMax
character(len=2) :: Unt

l = len(Location)
if (l < 1) then
  write(Lupri,'(/,A)') 'Largest memory block available @<UNKNOWN>:'
else
  write(Lupri,'(/,A,A,A)') 'Largest memory block available @',Location(1:l),':'
end if
call mma_maxDBLE(lMax)
call Cho_Word2Byte(lMax,8,dlMax,Unt)
write(Lupri,'(3X,I10,A,F10.3,A,A)') lMax,' 8-byte words; ',dlMax,' ',Unt

end subroutine Cho_PrtMaxMem
