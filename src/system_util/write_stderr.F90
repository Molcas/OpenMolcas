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

subroutine write_stderr(msg)

use Para_Info, only: MyRank
use Definitions, only: u0

implicit none
character(len=*), intent(in) :: msg

write(u0,'(a,i6,a,1x,a)') '[ process ',myrank,']:',trim(msg)
call xflush(u0)

end subroutine write_stderr
