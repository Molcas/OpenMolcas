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

subroutine PckInt_mck(abab,nZeta,nab,ab,Zeta)

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 abab(nZeta,nab,nab), ab(nZeta), Zeta(nZeta)

! Integrals

do iZeta=1,nZeta
  xMax = 0.0d0
  do iab=1,nab
    xTest = abs(abab(iZeta,iab,iab))
    if (xTest > xMax) then
      xMax = xTest
    end if
  end do
  ab(iZeta) = sqrt(xMax)
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Zeta)

end subroutine PckInt_mck
