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

subroutine Hess_Vec(nAtoms,Hess,EVec,mAtoms,nDim)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAtoms, mAtoms, nDim
real(kind=wp) :: Hess(3*nAtoms*(3*nAtoms+1)/2), EVec((3*mAtoms)**2)
integer(kind=iwp) :: iElem, ij, iQ, nQ
real(kind=wp) :: rZ
real(kind=wp), parameter :: ThrD = 1.0e-10_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the eigenvalues and the eigenvectors of the Hessian

!define _DEBUGPRINT_

! Set up a unit matrix

nQ = nDim
call dcopy_(nQ*nQ,[Zero],0,EVec,1)
call dcopy_(nQ,[One],0,EVec,nQ+1)

! Compute eigenvalues and eigenvectors

call NIDiag_new(Hess,EVec,nQ,nQ)
call JacOrd(Hess,EVec,nQ,nQ)

do iQ=1,nQ
  ! Fix standard direction.
  rZ = Zero
  do iElem=1,nQ
    ij = (iQ-1)*nQ+iElem
    if (abs(EVec(ij)) > abs(rZ)+ThrD) rZ = EVec(ij)
  end do
  ij = (iQ-1)*nQ+1
  if (rZ < Zero) call DScal_(nQ,-One,EVec(ij),1)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' Eigenvectors','(12f6.2)',EVec,nDim,nDim)
call TriPrt(' Eigenvalues','(12E8.2)',Hess,nDim)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Hess_Vec
