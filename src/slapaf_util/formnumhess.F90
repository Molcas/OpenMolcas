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

subroutine FormNumHess(nIter,nInter,Delta,nAtom,iNeg,DipM)

use Slapaf_Info, only: Cubic, Curvilinear, Degen, mTROld, dqInt, qInt, Shift, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIter, nInter, nAtom
real(kind=wp), intent(in) :: Delta, DipM(3,nIter)
integer(kind=iwp), intent(out) :: iNeg
#include "print.fh"
integer(kind=iwp) :: i, ii, ij, iPrint, iRout, mTR, nDim, nKtB
real(kind=wp) :: rDum(1)
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: dDipM(:), Degen2(:), FEq(:,:,:), H(:), HB(:), Hx(:), IRInt(:), KtB(:)

!                                                                      *
!***********************************************************************
!                                                                      *
mTR = mTROld
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 182
iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(dDipM,3*(nInter+mTR),Label='dDipM')
dDipM(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the Hessian matrix via a 2-point numerical differentiation.

call mma_allocate(H,nInter**2,Label='H')
call mma_allocate(FEq,merge(nInter,0,Cubic),nInter,nInter,Label='FEq')

call NmHess(Shift,nInter,dqInt,nIter,H,Delta,qInt,FEq,Cubic,DipM,dDipM)

write(u6,*)
write(u6,*) ' Numerical differentiation is finished!'
if (iPrint >= 98) call RecPrt(' Numerical force constant matrix',' ',H,nInter,nInter)

call Add_Info('Numerical Hessian',H,nInter**2,2)
call Put_dArray('Hss_Q',H,nInter**2)
call Put_dArray('Hss_upd',rDum,0)

! That is in internal coordinates, now transform it to Cartesians
! d^2E/dx^2 = dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ

call Qpg_dArray('KtB',Found,nKtB)
if (.not. Found) call Abend()
nDim = nKtB/nInter
call mma_allocate(KtB,nDim*nInter,Label='KtB')
call mma_allocate(HB,nDim*nInter,Label='HB')
call mma_allocate(Hx,nDim**2,Label='Hx')
call mma_allocate(Degen2,nDim,Label='Degen2')
call Get_dArray('KtB',KtB,nKtB)
call DGeMM_('N','T',nInter,nDim,nInter,One,H,nInter,KtB,nDim,Zero,HB,nInter)
call DGeMM_('T','T',nDim,nDim,nInter,One,HB,nInter,KtB,nDim,Zero,Hx,nDim)
i = 0
do ii=1,nAtom
  do ij=1,3
    if (Smmtrc(ij,ii)) then
      i = i+1
      Degen2(i) = Degen(ij,ii)
    end if
  end do
end do

! Compute and add the d^2Q/dx^2 dE/dQ part
if (Curvilinear) call dBuu(Degen2,nInter,nDim,dqInt(1,1),Hx,.true.)

call Put_dArray('Hss_X',Hx,nDim**2)
call mma_deallocate(KtB)
call mma_deallocate(HB)
call mma_deallocate(Hx)
call mma_deallocate(Degen2)
call mma_deallocate(H)

if (Cubic) then
  call RecPrt(' Numerical cubic force constant matrix',' ',FEq,nInter**2,nInter)
  call Add_Info('Numerical anharm. cons.',FEq,nInter**3,2)
end if
call mma_deallocate(FEq)
!                                                                      *
!***********************************************************************
!                                                                      *
! Do an harmonic frequency analysis

call mma_allocate(IRInt,nInter+mTR,Label='IRInt')

call HrmFrq(nAtom,nInter,iNeg,dDipM,mTR,DipM,IRInt)

call Add_Info('Numerical IR Intensities',IRInt,nInter,2)
call mma_deallocate(IRInt)
write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(dDipM)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine FormNumHess
