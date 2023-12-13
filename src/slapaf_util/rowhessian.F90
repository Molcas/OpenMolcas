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
! Copyright (C) Giovanni Ghigo                                         *
!***********************************************************************

subroutine RowHessian(nIter,nInter,Delta)
!***********************************************************************
!                                                                      *
! Object: Numerical estimation of single rows and columns of Hessian   *
! Called from: RlxCtl when lRowH=.true. & iter == NmIter               *
! Author: Giovanni Ghigo, University of Torino, Italy                  *
!                                                                      *
!***********************************************************************

use Slapaf_Info, only: dqInt, mRowH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIter, nInter
real(kind=wp), intent(in) :: Delta
integer(kind=iwp) :: iInter, iRowH, jInter
real(kind=wp) :: dElement, rDum(1)
real(kind=wp), allocatable :: H(:,:)

if (.not. allocated(mRowH)) then
  write(u6,*) 'RowHessian: .NOT.Allocated(mRowH)'
  call Abend()
end if

call mma_allocate(H,nInter,nInter,Label='H')
call Get_dArray('Hss_Q',H,nInter**2)
call Put_dArray('Hss_upd',rDum,0)

#ifdef _DEBUGPRINT_
write(u6,*) 'RowHessian:'
call RecPrt('Initial Hessian',' ',H,nInter,nInter)
call RecPrt('Gradient  dqInt:','(10F9.6)',dqInt,nInter,nIter)
#endif

! Evaluate the Hessian

do iRowH=1,size(mRowH)
  iInter = mRowH(iRowH)
  if (iInter > nIter) then
    write(u6,*) 'RowHessian: iIter>nIter'
    call Abend()
  end if
  H(iInter,:) = (dqInt(:,1)-dqInt(:,iRowH+1))/Delta
  H(:,iInter) = H(iInter,:)
end do

! Symmetrize

do iInter=1,nInter
  do jInter=1,nInter
    dElement = (H(iInter,jInter)+H(jInter,iInter))*Half
    H(iInter,jInter) = dElement
    H(jInter,iInter) = dElement
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Final Hessian',' ',H,nInter,nInter)
#endif
call Put_dArray('Hss_Q',H,nInter**2)
call mma_deallocate(H)

return

end subroutine RowHessian
