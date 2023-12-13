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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************

!----------------------------------------------------------------------*
! A function that will return the double factorial. We do not expect   *
! big numbers, so we do it brute-force. Observe that N must be odd, but*
! to skip the if-sentence, we assume that the one who calls this       *
! function has seen to that.                                           *
!----------------------------------------------------------------------*
function iDubFac(N)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iDubFac
integer(kind=iwp), intent(in) :: N
integer(kind=iwp) :: k

iDubFac = 1
do k=3,N,2
  iDubFac = iDubFac*k
end do

return

end function iDubFac
