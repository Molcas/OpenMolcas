!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Rmat,Debug)
! Thomas Bondo Pedersen, December 2005.
!
! Purpose: compute the gradient of the Pipek-Mezey functional.

use Constants, only: Zero, Four
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: GradNorm, Rmat(nOrb2Loc,nOrb2Loc)
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: i, iAtom, j
real(kind=wp) :: Fun, Rjj

RMat(:,:) = Zero
do iAtom=1,nAtoms
  do j=1,nOrb2Loc
    Rjj = PA(j,j,iAtom)
    do i=1,nOrb2Loc
      Rmat(i,j) = Rmat(i,j)+PA(i,j,iAtom)*Rjj
    end do
  end do
end do

GradNorm = Zero
do i=1,nOrb2Loc-1
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm+(Rmat(i,j)-Rmat(j,i))**2
  end do
end do
GradNorm = Four*sqrt(GradNorm)

if (Debug) then
  Fun = Zero
  do i=1,nOrb2Loc
    Fun = Fun+Rmat(i,i)
  end do
  write(u6,*) 'GetGrad_PM: functional = Tr(R) = ',Fun
end if

end subroutine GetGrad_PM
