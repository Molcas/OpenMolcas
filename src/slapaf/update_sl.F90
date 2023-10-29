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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine Update_sl(Step_Trunc,nWndw,kIter)
!***********************************************************************
!                                                                      *
!     Object: to update coordinates                                    *
!                                                                      *
!    OutPut:                                                           *
!      Step_Trunc     : character label to denote truncation type      *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh                                             *
!             2000                                                     *
!***********************************************************************

use Slapaf_Info, only: Beta, Beta_disp, iter, NmIter, qInt, Shift
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
character, intent(inout) :: Step_Trunc
integer(kind=iwp), intent(in) :: nWndw, kIter
integer(kind=iwp) :: iOpt_RS, iter_, nQQ
real(kind=wp) :: Dummy(1), qBeta, qBeta_Disp
logical Kriging_Hessian, Hide
real(kind=wp), allocatable :: t_qInt(:,:), t_Shift(:,:), tmp(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nQQ = size(qInt,1)

#ifdef _DEBUGPRINT_
call RecPrt('Update_sl: qInt',' ',qInt,nQQ,Iter)
call RecPrt('Update_sl: Shift',' ',Shift,nQQ,Iter-1)
#endif

iOpt_RS = 0
qBeta = Beta
qBeta_Disp = Beta_Disp
!                                                                      *
!***********************************************************************
!                                                                      *
call Mk_Hss_Q()
Kriging_Hessian = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Select between numerical evaluation of the Hessian or a molcular
! geometry optimization.

Hide = Step_Trunc == '#'
Step_Trunc = ' '
if ((iter == NmIter) .and. (NmIter /= 1)) then
!                                                                      *
!***********************************************************************
!                                                                      *
! On the first iteration after a numerical evaluation of the
! Hessian we like the step to be relative to the initial structure.

# ifdef _DEBUGPRINT_
  write(u6,*) 'UpDate_SL: first iteration'
# endif
  iter_ = 1
  call mma_Allocate(t_Shift,nQQ,size(Shift,2),Label='t_Shift')
  t_Shift(:,:) = Shift(:,:)
  Shift(:,:) = Zero

  call mma_Allocate(t_qInt,nQQ,size(qInt,2),Label='t_qInt')
  t_qInt(:,:) = qInt(:,:)
  qInt(:,:) = Zero
  qInt(:,1) = t_qInt(:,1)

  call Update_inner(iter_,Beta,Beta_Disp,Step_Trunc,nWndw,kIter,Kriging_Hessian,qBeta,iOpt_RS,.true.,iter_,qBeta_Disp,Hide)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Move new coordinates to the correct position and compute the
  ! corresponding shift.

  call mma_allocate(Tmp,size(qInt,1),Label='tmp')
  tmp(:) = qInt(:,2)
  qInt(:,:) = t_qInt(:,:)
  qInt(:,iter+1) = tmp(:)
  Shift(:,:) = t_Shift(:,:)
  Shift(:,iter) = tmp(:)-qInt(:,iter)

  call mma_deallocate(tmp)
  call mma_deallocate(t_qInt)
  call mma_deallocate(t_Shift)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Conventional optimization.

  call Update_inner(iter,Beta,Beta_Disp,Step_Trunc,nWndw,kIter,Kriging_Hessian,qBeta,iOpt_RS,.true.,iter,qBeta_Disp,Hide)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the shift in internal coordinates basis.

#ifdef _DEBUGPRINT_
call RecPrt('Shifts in internal coordinate basis / au or rad',' ',Shift,nQQ,Iter)
call RecPrt('qInt in internal coordinate basis / au or rad',' ',qInt,nQQ,Iter+1)
#endif

! Remove unneeded fields from the runfile
Dummy(1) = -Zero
call Put_dArray('BMxOld',Dummy,0)
call Put_dArray('TROld',Dummy,0)

return

end subroutine Update_sl
