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

subroutine SECNE(A,B,C,NAL,NBL,NSIJ,IFT)

implicit real*8(A-H,O-Z)
dimension A(NAL,NBL), B(NBL,NAL), C(NBL,NAL)

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
! Avoid unused argument warnings
if (.false.) call Unused_integer(NSIJ)

end subroutine SECNE
