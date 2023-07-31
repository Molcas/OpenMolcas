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
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

subroutine NewH(nInter,nIter,dq_orig,g,H,iOptH,mIter)
!***********************************************************************
!                                                                      *
! Object: Driver for inverse Hessian update.                           *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '95                                              *
!***********************************************************************

use Slapaf_Info, only: HUpMet
use NewH_Mod, only: UpdMask
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter, nIter, iOptH, mIter
real(kind=wp), intent(in) :: dq_orig(nInter,nIter), g(nInter,mIter+1)
real(kind=wp), intent(inout) :: H(nInter,nInter)
integer(kind=iwp) :: i
logical(kind=iwp) :: DoMask
real(kind=wp), allocatable :: dg(:), dq(:,:), gi(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' NewH: lIter=',nIter
call RecPrt(' NewH: dq_orig',' ',dq_orig,nInter,nIter)
call RecPrt(' NewH: g',' ',g,nInter,nIter)
call RecPrt(' NewH: H(Old)',' ',H,nInter,nInter)
write(u6,*) ' NewH: iOptH=',iOptH
#endif

! Branch out if the first iteration

if (nIter <= 1) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the MM mask
DoMask = .false.
if (allocated(UpdMask)) then
  if (size(UpdMask) == nInter) DoMask = .true.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(dg,nInter,label='dg')
call mma_allocate(gi,nInter,label='gi')
call mma_allocate(dq,nInter,nIter,label='dq')
dq(:,:) = dq_orig(:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the difference between the gradients. Observe the order!
! This since we are storing the forces rather than the gradients,
! big mistake!!!

dg(:) = g(:,nIter-1)-g(:,nIter)
if (DoMask) then
  do i=1,nInter
    if (UpdMask(i) /= 0) then
      dg(i) = Zero
      dq(i,nIter-1) = Zero
    end if
  end do
end if
#ifdef _DEBUGPRINT_
call RecPrt(' NewH: dg',' ',dg,nInter,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the update

if (btest(iOptH,3)) then

  ! No update

  HUpMet = ' None '

else if (btest(iOptH,0)) then

  !Fletcher (or Meyer) update

  HUpMet = '  F   '
  write(u6,*) 'Deleted option in NewH'
  call Abend()

else if (btest(iOptH,1)) then

  ! Broyden-Powel Symmetric Rank-2 update

  HUpMet = '  BP  '
  write(u6,*) 'Deleted option in NewH'
  call Abend()

else if (btest(iOptH,2)) then

  ! Broyden-Fletcher-Goldfarb-Shanno update

  HUpMet = ' BFGS '
  call DFP(H,nInter,gi,dq(1,nIter-1),dg)

else if (btest(iOptH,4)) then

  ! Murtagh-Sargent-Powell update

  HUpMet = ' MSP  '
  call dGeMV_('N',nInter,nInter,-One,H,nInter,dq(1,nIter-1),1,One,dg,1)
# ifdef _DEBUGPRINT_
  call RecPrt(' NewH: gamma',' ',dg,nInter,1)
# endif
  call MSP(H,dg,dq(1,nIter-1),nInter)

else if (btest(iOptH,5)) then

  ! EU update

  HUpMet = '  EU  '

  ! Some precalculations:
  gi(:) = -g(:,nIter-1)
  if (DoMask) then
    do i=1,nInter
      if (UpdMask(i) /= 0) gi(i) = Zero
    end do
  end if

  call EU(dq(1,nIter-1),dg,gi,H,nInter)

else if (btest(iOptH,6)) then

  ! TS_BFGS update

  HUpMet = 'TSBFGS'
  call TS_BFGS(dq(1,nIter-1),dg,H,nInter)

else

  call WarningMessage(2,'Error in NewH')
  write(u6,*) ' Improper value of iOptH:',iOptH
  call Abend()

end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' NewH:  H(New)',' ',H,nInter,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(dg)
call mma_deallocate(gi)
call mma_deallocate(dq)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine NewH
