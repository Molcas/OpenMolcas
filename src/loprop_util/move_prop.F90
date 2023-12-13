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

subroutine Move_Prop(rMP,EC,lMax,nElem,nAtoms,nPert,nij,iANr,Bond_Threshold)
! Distributes the contributions from the bonds that doesn't fulfill the requirement
! Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
! two atoms involved in the bond.
!
! Polarizabilities are moved in the Move_Polar subroutine!

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax, nElem, nAtoms, nPert, nij, iANr(nAtoms)
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1)
real(kind=wp), intent(in) :: EC(3,nij), Bond_Threshold
integer(kind=iwp) :: iAtom, iElem, ii, ij, iPert, jAtom, jj
logical(kind=iwp) :: Bond_OK, Check_Bond

do iAtom=2,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=1,iAtom-1
    jj = jAtom*(jAtom+1)/2
    Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),iANr(jAtom),Bond_Threshold)
    if (.not. Bond_OK) then
      ij = iAtom*(iAtom-1)/2+jAtom
      do iPert=0,nPert-1

        ! First move half of the bond properties to iAtom

        do iElem=0,nElem-1
          rMP(ij,iElem,iPert) = rMP(ij,iElem,iPert)*Half
        end do

        call ReExpand(rMP(1,0,iPert),nij,nElem,EC(1,ij),EC(1,ii),ij,lMax)

        do iElem=0,nElem-1
          rMP(ii,iElem,iPert) = rMP(ii,iElem,iPert)+rMP(ij,iElem,iPert)
        end do

        ! Then move the other half of the bond properties to jAtom

        call ReExpand(rMP(1,0,iPert),nij,nElem,EC(1,ii),EC(1,jj),ij,lMax)

        do iElem=0,nElem-1
          rMP(jj,iElem,iPert) = rMP(jj,iElem,iPert)+rMP(ij,iElem,iPert)
        end do

        ! Set local properties to zero

        rMP(ij,:,iPert) = Zero
      end do
    end if
  end do
end do

return

end subroutine Move_Prop
