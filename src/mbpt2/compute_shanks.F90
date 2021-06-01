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

subroutine Compute_Shanks(E1,E2,EOrb,lthEOrb,nBas,nFro,nOcc,nSym,E0,Shanks1_E)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lthEOrb, nSym, nBas(nSym), nFro(nSym), nOcc(nSym)
real(kind=wp), intent(in) :: E1, E2, EOrb(lthEOrb)
real(kind=wp), intent(out) :: E0, Shanks1_E
integer(kind=iwp) :: ioff, iorb, iSym, jorb, nOrb
real(kind=wp) :: PotNuc

E0 = Zero
ioff = 0
do iSym=1,nSym
  nOrb = nFro(iSym)+nOcc(iSym)
  do iorb=1,nOrb
    jorb = ioff+iorb
    E0 = E0+EOrb(jorb)
  end do
  ioff = ioff+nBas(iSym)
end do
E0 = Two*E0

call Peek_dScalar('PotNuc',PotNuc)
E0 = E0+PotNuc

! Shanks formula

Shanks1_E = (E2*E0-E1**2)/(E2-Two*E1+E0)

return

end subroutine Compute_Shanks
