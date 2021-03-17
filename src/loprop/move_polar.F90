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

subroutine Move_Polar(Polar,EC,nAtoms,nij,iANr,Bond_Threshold)

! Distributes the contributions from the bonds that doesn't fullfill the requiremen
! Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom)) to the
! two atoms involved in the bond.
!
! Multipole moments are moved in the Move_Prop subroutine!

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nij, iANr(nAtoms)
real(kind=wp), intent(inout) :: Polar(6,nij)
real(kind=wp), intent(in) :: EC(3,nij), Bond_Threshold
integer(kind=iwp) :: iAtom, ii, ij, jAtom, jj
logical(kind=iwp) :: Bond_Ok
logical(kind=iwp), external :: Check_Bond

do iAtom=2,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=1,iAtom-1
    jj = jAtom*(jAtom+1)/2
    Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),iANr(jAtom),Bond_Threshold)
    ij = iAtom*(iAtom-1)/2+jAtom
    if (.not. Bond_Ok) then

      ! Move half of the bond polarizabilities to each atom

      Polar(:,ii) = Polar(:,ii)+Half*Polar(:,ij)
      Polar(:,jj) = Polar(:,jj)+Half*Polar(:,ij)
      Polar(:,ij) = Zero
    end if
  end do
end do

return

end subroutine Move_Polar
