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

subroutine SECEQ(A,B,C,NAL,IFT,FAC)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAL, IFT
real(kind=wp), intent(in) :: A(NAL,NAL), B(NAL,NAL), FAC
real(kind=wp), intent(out) :: C((NAL*(NAL+1))/2)
integer(kind=iwp) :: IAB, NA, NB

if (IFT == 0) then
  IAB = 0
  do NA=1,NAL
    do NB=1,NA-1
      IAB = IAB+1
      C(IAB) = B(NB,NA)+A(NA,NB)
    end do
    IAB = IAB+1
    C(IAB) = FAC*A(NA,NA)
  end do
else
  IAB = 0
  do NA=1,NAL
    do NB=1,NA-1
      IAB = IAB+1
      C(IAB) = B(NB,NA)-A(NA,NB)
    end do
    IAB = IAB+1
    C(IAB) = Zero
  end do
end if

return

end subroutine SECEQ
