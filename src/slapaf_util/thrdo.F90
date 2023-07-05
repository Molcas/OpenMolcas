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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 g(nInter), A(nInter,nInter), e(nInter,2)
logical Fail

i0 = 1
i1 = 2
iter = 0
iStep = 0
Fail = .true.

! Compile the error vector, starting in position 1
! Observe that g is the force and NOT the gradient

call dcopy_(nInter,g,1,e(1,i0),1)
call DPOTRS('U',nInter,1,A,nInter,e(1,i0),nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'ThrdO(DPOTRS): iRC=',iRC
  call Abend()
end if

call RecPrt(' ThrdO: e(0)',' ',e(1,i0),nInter,1)

Thrd = 1.D-6
iterMx = 40

do
  iStep = iStep+1

  ! Newton-Raphson scheme

  call dcopy_(nInter,g,1,e(1,i1),1)
  call DPOTRS('U',nInter,1,A,nInter,e(1,i1),nInter,iRC)
  if (iRC /= 0) then
    write(6,*) 'ThrdO(DPOTRS): iRC=',iRC
    call Abend()
  end if

  !call RecPrt(' ThrdO: e',' ',e(1,i1),nInter,1)
  iter = iter+1

  ! Check if the error vectors are self consistent.

  Test = Zero
  do i=1,nInter
    diff = abs(e(i,i0)-e(i,i1))
    if (diff > Test) Test = diff
  end do
  !write(6,*) iter,diff

  if (iter > iterMx) then
    call WarningMessage(1,'ThrdO: Exceeded max iterations')
    return
  end if

  if (Test < Thrd) then
    ! Copy converged vectors to slot 1
    if (i1 /= 1) call dcopy_(nInter,e(1,i1),1,e(1,1),1)
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
