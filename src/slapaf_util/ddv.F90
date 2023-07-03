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

subroutine ddV(Cart,mTtAtm,Hess,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)

use Symmetry_Info, only: VarR, VarT

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
real*8 Cart(3,mTtAtm+nHidden), Hess((3*mTtAtm)*(3*mTtAtm+1)/2)
integer iANr(mTtAtm+nHidden), iTabBonds(3,nBonds), iTabAtoms(2,0:nMax,mTtAtm+nHidden)
real*8, allocatable :: HBig(:)
logical VRSave, VTSave

! Temporary big hessian
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (nHidden > 0) then
  nTot = mTtAtm+nHidden
  call mma_allocate(HBig,(3*nTot)*(3*nTot+1)/2,Label='HBig')

  ! Temporary turn on the translational/rotational invariance

  VRSave = VarR
  VTSave = VarT
  VarR = .false.
  VarT = .false.
  call ddV_(Cart,nTot,HBig,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
  VarR = VRSave
  VarT = VTSave
  call dCopy_((3*mTtAtm)*(3*mTtAtm+1)/2,HBig,1,Hess,1)
# ifdef _DEBUGPRINT_
  write(6,*) 'DDV: Improved Hessian'
  call RecPrt('Coord (with hidden atoms):',' ',Cart,3,nTot)
  call TriPrt('Hessian (hidden atoms):',' ',HBig,3*nTot)
  call TriPrt('Hessian (normal):',' ',Hess,3*mTtAtm)
# endif
  call mma_deallocate(HBig)
else
  call ddV_(Cart,mTtAtm,Hess,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
end if

end subroutine ddV
