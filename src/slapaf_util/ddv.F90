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

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: VarR, VarT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mTtAtm, nHidden, iANr(mTtAtm+nHidden), nBonds, iTabBonds(3,nBonds), nMax, &
                                 iTabAtoms(2,0:nMax,mTtAtm+nHidden)
real(kind=wp), intent(in) :: Cart(3,mTtAtm+nHidden)
real(kind=wp), intent(out) :: Hess(nTri_Elem(3*mTtAtm))
integer(kind=iwp) :: nTot
logical(kind=iwp) :: VRSave, VTSave
real(kind=wp), allocatable :: HBig(:)

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
  call mma_allocate(HBig,nTri_Elem(3*nTot),Label='HBig')

  ! Temporary turn on the translational/rotational invariance

  VRSave = VarR
  VTSave = VarT
  VarR = .false.
  VarT = .false.
  call ddV_(Cart,nTot,HBig,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
  VarR = VRSave
  VarT = VTSave
  Hess(:) = HBig(1:nTri_Elem(3*mTtAtm))
# ifdef _DEBUGPRINT_
  write(u6,*) 'DDV: Improved Hessian'
  call RecPrt('Coord (with hidden atoms):',' ',Cart,3,nTot)
  call TriPrt('Hessian (hidden atoms):',' ',HBig,3*nTot)
  call TriPrt('Hessian (normal):',' ',Hess,3*mTtAtm)
# endif
  call mma_deallocate(HBig)
else
  call ddV_(Cart,mTtAtm,Hess,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
end if

end subroutine ddV
