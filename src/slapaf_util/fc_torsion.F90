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

function FC_Torsion(jSet,Coor,iTabAtoms,nMax,mAtoms)

use Slapaf_Info, only: Covalent_Bond
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FC_Torsion
integer(kind=iwp), intent(in) :: jSet(4), nMax, mAtoms, iTabAtoms(2,0:nMax,mAtoms)
real(kind=wp), intent(in) :: Coor(3,4)
integer(kind=iwp) :: i, iAtom, iAtom_, iHit, ii, jAtom, jAtom_
real(kind=wp) :: Gij(3)
integer(kind=iwp), parameter :: iCase(2,3) = reshape([1,2,2,3,3,4],[2,3])

FC_Torsion = Zero
do ii=1,3
  iHit = 0
  iAtom = jSet(iCase(1,ii))
  jAtom = jSet(iCase(2,ii))

  ! Loop over the neighbor atoms of iAtom to see if jAtom is a neighbor. If
  ! it is check that the bond is classified as a covalent bond.

  do i=1,iTabAtoms(1,0,iAtom)
    if ((jAtom == iTabAtoms(1,i,iAtom)) .and. (iTabAtoms(2,i,iAtom) == Covalent_Bond)) then
      iHit = 1
      iAtom_ = iCase(1,ii)
      jAtom_ = iCase(2,ii)
      Gij(ii) = One/sqrt((Coor(1,iAtom_)-Coor(1,jAtom_))**2+(Coor(2,iAtom_)-Coor(2,jAtom_))**2+(Coor(3,iAtom_)-Coor(3,jAtom_))**2)
    end if
  end do
  if (iHit == 0) return
end do
FC_Torsion = Gij(1)*Gij(2)*Gij(3)

return

end function FC_Torsion
