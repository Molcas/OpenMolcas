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

subroutine Move_Xhole(rX,EC,nAtoms,nij,iANr,Bond_Threshold)
! Distribute the XHole-integrals if the bond lengths does not fulfill:
! Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
! two atoms involved in the bond.
! Ripped from move_prop.

implicit real*8(A-H,O-Z)
#include "real.fh"
dimension rX(nij), EC(3,nij)
dimension iAnr(nAtoms)
logical Bond_OK, Check_Bond

do iAtom=2,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=1,iAtom-1
    jj = jAtom*(jAtom+1)/2
    Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),iANr(jAtom),Bond_Threshold)
    if (.not. Bond_OK) then
      ij = iAtom*(iAtom-1)/2+jAtom

      ! First move half of the bond properties to iAtom

      rX(ij) = rX(ij)*Half
      rX(ii) = rX(ii)+rX(ij)

      ! Then move the other half of the bond properties to jAtom

      rX(jj) = rX(jj)+rX(ij)

      ! Set local properties to zero

      rX(ij) = Zero
    end if
  end do
end do

return

end subroutine Move_Xhole
