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

subroutine lattXPol(Grid,nGrid,nGrid_Eff,PolEff,DipEff,XF,nXF,nOrd_XF,nPolComp)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrid, nXF, nOrd_XF, nPolComp
integer(kind=iwp), intent(inout) :: nGrid_Eff
real(kind=wp), intent(inout) :: Grid(3,nGrid), PolEff(nPolComp,nGrid), DipEff(nGrid)
real(kind=wp), intent(in) :: XF(nXF)
integer(kind=iwp) :: Inc, iOrdOp, iXF, j
! Statement function for Cartesian index
integer(kind=iwp) :: ixyz, nElem
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

! Calculate number of entries per XFIELD point
Inc = 3
do iOrdOp=0,nOrd_XF
  Inc = Inc+nElem(iOrdOp)
end do
Inc = Inc+6  !iXpolType always > 0

! Insert XFIELD polarisabilities into Grid
do iXF=1,nXF
  nGrid_Eff = nGrid_Eff+1
  do j=1,nPolComp
    PolEff(j,nGrid_Eff) = XF(iXF*Inc-6+j)
  end do
  DipEff(nGrid_Eff) = Zero
  do j=1,3
    Grid(j,nGrid_Eff) = XF((iXF-1)*Inc+j)
  end do
end do

return

end subroutine lattXPol
