!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine zlcase(line)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: line
integer(kind=iwp) :: k, idx

do k=1,len(line)
   idx = ichar(line(k:k))
   !if ((97 <= idx).and.(idx <= 122)) line(k:k) = char(idx-32)
   if ((65 <= idx).and.(idx <= 90)) line(k:k) = char(idx+32)
end do

return

end subroutine zlcase
