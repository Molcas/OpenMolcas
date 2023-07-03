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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 A(3), Tx(3)
integer jStab(0:nStab-1)

do iStab=0,nStab-1
  if ((A(1) /= Zero) .and. (iand(jStab(iStab),1)) /= 0) Go To 10
  if ((A(2) /= Zero) .and. (iand(jStab(iStab),2)) /= 0) Go To 10
  if ((A(3) /= Zero) .and. (iand(jStab(iStab),4)) /= 0) Go To 10
  if (iand(jStab(iStab),1) /= 0) Tx(1) = Zero
  if (iand(jStab(iStab),2) /= 0) Tx(2) = Zero
  if (iand(jStab(iStab),4) /= 0) Tx(3) = Zero
10 continue
end do

return

end subroutine NonSym
