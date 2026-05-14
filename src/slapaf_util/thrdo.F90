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

subroutine ThrdO(nInter,g,A,e,Fail)
!***********************************************************************
!                                                                      *
! Object: to find the error vector given the gradient, and the Hessian *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: g(nInter), A(nInter,nInter)
real(kind=wp), intent(out) :: e(nInter,2)
logical(kind=iwp), intent(out) :: Fail
integer(kind=iwp) :: i, i0, i1, iRC, iStep, iter, iterMx, itmp
real(kind=wp) :: diff, Test
real(kind=wp), parameter :: Thrd = 1.0e-6_wp

i0 = 1
i1 = 2
iter = 0
iStep = 0
Fail = .true.

! Compile the error vector, starting in position 1
! Observe that g is the force and NOT the gradient

e(:,i0) = g(:)
call DPOTRS('U',nInter,1,A,nInter,e(:,i0),nInter,iRC)
if (iRC /= 0) then
  write(u6,*) 'ThrdO(DPOTRS): iRC=',iRC
  call Abend()
end if

call RecPrt(' ThrdO: e(0)',' ',e(:,i0),nInter,1)

iterMx = 40

do
  iStep = iStep+1

  ! Newton-Raphson scheme

  e(:,i1) = g(:)
  call DPOTRS('U',nInter,1,A,nInter,e(:,i1),nInter,iRC)
  if (iRC /= 0) then
    write(u6,*) 'ThrdO(DPOTRS): iRC=',iRC
    call Abend()
  end if

  !call RecPrt(' ThrdO: e',' ',e(:,i1),nInter,1)
  iter = iter+1

  ! Check if the error vectors are self consistent.

  Test = Zero
  do i=1,nInter
    diff = abs(e(i,i0)-e(i,i1))
    if (diff > Test) Test = diff
  end do
  !write(u6,*) iter,diff

  if (iter > iterMx) then
    call WarningMessage(1,'ThrdO: Exceeded max iterations')
    return
  end if

  if (Test < Thrd) then
    ! Copy converged vectors to slot 1
    if (i1 /= 1) e(:,1) = e(:,i1)
    if (iStep == 10) then
      call RecPrt(' ThrdO: e(Final)',' ',e,nInter,1)
      Fail = .false.
      exit
    else
      iter = 0
    end if
  else
    ! If not change vector slot
    itmp = i0
    i0 = i1
    i1 = itmp
    iStep = iStep-1
  end if
end do
if (Fail) then
  call WarningMessage(2,'Error in ThrdO')
  call Abend()
end if

end subroutine ThrdO
