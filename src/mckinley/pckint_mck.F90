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

!#define _DEBUGPRINT_
subroutine PckInt_mck(abab,nZeta,nab,ab)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, nab
real(kind=wp), intent(in) :: abab(nZeta,nab,nab)
real(kind=wp), intent(out) :: ab(nZeta)
integer(kind=iwp) :: iab, iZeta
real(kind=wp) :: xMax, xTest

! Integrals

do iZeta=1,nZeta
  xMax = Zero
  do iab=1,nab
    xTest = abs(abab(iZeta,iab,iab))
    if (xTest > xMax) xMax = xTest
  end do
  ab(iZeta) = sqrt(xMax)
end do

#ifdef _DEBUGPRINT_
call recprt('abab',' ',abab,nZeta,nab*2)
call recprt('ab',' ',ab,nZeta,1)
#endif

return

end subroutine PckInt_mck
