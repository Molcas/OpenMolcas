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

subroutine QMPosition(EHam,Cordst,Coord,Forcek,dLJrep,Ract,iQ_Atoms)

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: EHam
integer(kind=iwp), intent(in) :: iQ_Atoms
real(kind=wp), intent(in) :: Cordst(3,iQ_Atoms), Coord(3), Forcek, dLJrep, Ract
integer(kind=iwp) :: iAt
real(kind=wp) :: dDepart, Diff, R

! First the harmonic potential that keeps QM close to centre.

dDepart = (Cordst(1,1)-Coord(1))**2+(Cordst(2,1)-Coord(2))**2+(Cordst(3,1)-Coord(3))**2
EHam = Forcek*Half*dDepart

! Second the repulsion with boundary that keeps QM away from boundary.

do iAt=1,iQ_Atoms
  R = Cordst(1,iAt)**2+Cordst(2,iAt)**2+Cordst(3,iAt)**2
  R = sqrt(R)
  Diff = Ract-R
  EHam = EHam+(dLJRep/Diff)**12
end do

return

end subroutine QMPosition
