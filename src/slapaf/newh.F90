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

use Slapaf_Parameters, only: HUpMet
use NewH_Mod, only: UpdMask
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nInter, nIter, iOptH, mIter
real(kind=wp) :: dq_orig(nInter,nIter), g(nInter,mIter+1), H(nInter,nInter)
integer(kind=iwp) :: i
logical(kind=iwp) :: DoMask
real(kind=wp), allocatable :: dg(:), dq(:,:), gi(:)
! Statement function
logical(kind=iwp) :: Test
Test(i) = iand(iOptH,2**(i-1)) == 2**(i-1)

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
write(u6,*) ' NewH: Test(i)==',(Test(i),i=1,8)
#endif

! Branch out if the first iteration

if (nIter <= 1) Go To 999
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
call dcopy_(nInter*nIter,dq_orig,1,dq,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the difference between the gradients. Observe the order!
! This since we are storing the forces rather than the gradients,
! big mistake!!!

do i=1,nInter
  dg(i) = g(i,nIter-1)-g(i,nIter)
  if (DoMask) then
    if (UpdMask(i) /= 0) then
      dg(i) = Zero
      dq(i,nIter-1) = Zero
    end if
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt(' NewH: dg',' ',dg,nInter,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the update

if (Test(4)) then

  ! No update

  HUpMet = ' None '

else if (Test(1)) then

  !Fletcher (or Meyer) update

  HUpMet = '  F   '
  write(u6,*) 'Deleted option in NewH'
  call Abend()

else if (Test(2)) then

  ! Broyden-Powel Symmetric Rank-2 update

  HUpMet = '  BP  '
  write(u6,*) 'Deleted option in NewH'
  call Abend()

else if (Test(3)) then

  ! Broyden-Fletcher-Goldfarb-Shanno update

  HUpMet = ' BFGS '
  call DFP(H,nInter,gi,dq(1,nIter-1),dg)

else if (Test(5)) then

  ! Murtagh-Sargent-Powell update

  HUpMet = ' MSP  '
  call dGeMV_('N',nInter,nInter,-One,H,nInter,dq(1,nIter-1),1,One,dg,1)
# ifdef _DEBUGPRINT_
  call RecPrt(' NewH: gamma',' ',dg,nInter,1)
# endif
  call MSP(H,dg,dq(1,nIter-1),nInter)

else if (Test(7)) then

  ! EU update

  HUpMet = '  EU  '

  ! Some precalculations:
  do i=1,nInter
    gi(i) = -g(i,nIter-1)
    if (DoMask) then
      if (UpdMask(i) /= 0) gi(i) = Zero
    end if
  end do

  call EU(dq(1,nIter-1),dg,gi,H,nInter)

else if (Test(8)) then

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
999 continue
return

end subroutine NewH
