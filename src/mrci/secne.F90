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

subroutine SECNE(A,B,C,NAL,NBL,IFT)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAL, NBL, IFT
real(kind=wp), intent(in) :: A(NAL,NBL), B(NBL,NAL)
real(kind=wp), intent(out) :: C(NBL,NAL)
integer(kind=iwp) :: NA, NB

if (IFT == 0) then
  do NA=1,NAL
    do NB=1,NBL
      C(NB,NA) = B(NB,NA)+A(NA,NB)
    end do
  end do
else
  do NA=1,NAL
    do NB=1,NBL
      C(NB,NA) = B(NB,NA)-A(NA,NB)
    end do
  end do
end if

return

end subroutine SECNE
