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

subroutine Force(nFix,GrdX,nAtom,nInter,BMx,Iter,Grad,Lbl,Degen)

#ifdef _DEBUGPRINT_
use Slapaf_Info, only: AtomLbl
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nFix, nAtom, nInter, Iter
real(kind=wp), intent(inout) :: GrdX(3*nAtom), Grad(nInter,Iter)
real(kind=wp), intent(in) :: BMx(3*nAtom,3*nAtom), Degen(3*nAtom)
character(len=8), intent(in) :: Lbl(nInter)
real(kind=wp) :: Dummy(1)
real(kind=wp), allocatable :: Frc(:)

#ifdef _DEBUGPRINT_
call RecPrt('In Force:BMx ',' ',BMx,3*nAtom,nInter)
call RecPrt('In Force:Degen ',' ',Degen,1,3*nAtom)
call RecPrt('In Force:GrdX',' ',GrdX,3,nAtom)
#endif

! Frozen internal coordinates are introduced as constraints to the
! energy functional with the Lagrange multipliers method. The new
! energy functional will have zero gradient with respect to the
! frozen parameters.

call mma_allocate(Frc,3*nAtom,Label='Frc')

! Compute the norm of the cartesian force vector.
!
! |dE/dx|=Sqrt(dE/dx|u|dE/dx)

Frc(:) = Degen(:)*GrdX(:)
!                                                                      *
!***********************************************************************
!                                                                      *
! Solve the equation dq/dx dE/dq = dE/dx
!
! B dE/dq = dE/dx

call Eq_Solver('N',3*nAtom,nInter,1,BMx,.true.,Dummy,Frc,Grad(:,Iter))
#ifdef _DEBUGPRINT_
call RecPrt(' Internal Forces in au before FIXIC ',' ',Grad(:,Iter),nInter,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Remove gradient components in the constraints directions.

if (nFix /= 0) call Fixic(nFix,Grad(:,Iter),nInter,BMx,nAtom*3,GrdX,Lbl,Degen)

! Write cartesian symmetry distinct forces which will be relaxed.

#ifdef _DEBUGPRINT_
call PrList('Cartesian forces which will be relaxed hartree/bohr',AtomLbl,nAtom,GrdX,3,nAtom)
#endif

call mma_deallocate(Frc)

return

end subroutine Force
