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

subroutine NonSym(nStab,jStab,A,Tx)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nStab, jStab(0:nStab-1)
real(kind=wp), intent(in) :: A(3)
real(kind=wp), intent(inout) :: Tx(3)
integer(kind=iwp) :: iStab

do iStab=0,nStab-1
  if ((A(1) /= Zero) .and. btest(jStab(iStab),0)) cycle
  if ((A(2) /= Zero) .and. btest(jStab(iStab),1)) cycle
  if ((A(3) /= Zero) .and. btest(jStab(iStab),2)) cycle
  if (btest(jStab(iStab),0)) Tx(1) = Zero
  if (btest(jStab(iStab),1)) Tx(2) = Zero
  if (btest(jStab(iStab),2)) Tx(3) = Zero
end do

return

end subroutine NonSym
