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

subroutine Hess_Vec(nAtoms,Hess,EVec,nDim)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nDim
real(kind=wp), intent(inout) :: Hess(nTri_Elem(3*nAtoms))
real(kind=wp), intent(out) :: EVec(nDim,nDim)
integer(kind=iwp) :: iElem, iQ, nQ
real(kind=wp) :: rZ
real(kind=wp), parameter :: ThrD = 1.0e-10_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the eigenvalues and the eigenvectors of the Hessian

!#define _DEBUGPRINT_

! Set up a unit matrix

nQ = nDim
call unitmat(EVec,nQ)

! Compute eigenvalues and eigenvectors

call NIDiag_new(Hess,EVec,nQ,nQ)
call JacOrd(Hess,EVec,nQ,nQ)

do iQ=1,nQ
  ! Fix standard direction.
  rZ = Zero
  do iElem=1,nQ
    if (abs(EVec(iElem,iQ)) > abs(rZ)+ThrD) rZ = EVec(iElem,iQ)
  end do
  if (rZ < Zero) EVec(:,iQ) = -EVec(:,iQ)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' Eigenvectors','(12f6.2)',EVec,nDim,nDim)
call TriPrt(' Eigenvalues','(12ES9.2)',Hess,nDim)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Hess_Vec
