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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine NwShft()
!***********************************************************************
!                                                                      *
! Object: to numerically evaluate the molecular Hessian.               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             May '92                                                  *
!***********************************************************************

use Slapaf_Info, only: Delta, iter, qInt, Shift
use Constants, only: Zero, Two
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nInter
integer(kind=iwp) :: iCount, jCount, jInter, kCount, kInter, lInter

nInter = size(Shift,1)

# ifdef _DEBUGPRINT_
call RecPrt('NwShft  qInt',' ',qInt,nInter,iter)
call RecPrt('NwShft Shift',' ',Shift,nInter,iter-1)
# endif

! Compute the new shift

!write(u6,*) ' iter=',iter
if (iter < 2*nInter+1) then

  ! Shifts for the numerical Hessian

  jInter = (iter+1)/2
  Shift(:,iter) = Zero
  if (mod(iter,2) == 0) then
    Shift(jInter,iter) = -Two*Delta
  else
    ! Undo previous displacement
    if (jInter > 1) Shift(jInter-1,iter) = Delta
    Shift(jInter,iter) = Delta
  end if

else

  ! Shifts for the numerical cubic force constants

  iCount = (iter-2*nInter+3)/4
  !write(u6,*) ' iCount=',iCount
  jCount = 0
  lInter = 0
  outer: do kInter=1,nInter
    do lInter=1,kInter-1
      jCount = jCount+1
      if (jCount == iCount) exit outer
    end do
  end do outer
  if (lInter == 0) then
    call WarningMessage(2,'lInter == 0')
    call Abend()
  end if
  !write(u6,*) 'kInter, lInter=',kInter,lInter
  kCount = iter-2*nInter
  Shift(:,iter) = Zero
  ! Undo last change for numerical Hessian
  if (iCount == 1) Shift(nInter,iter) = Delta
  if (mod(kCount,4) == 1) then
    ! Undo change due to previous pair
    if (lInter /= 1) then
      Shift(kInter,iter) = Delta
      Shift(lInter-1,iter) = Delta
    else if ((lInter == 1) .and. (kInter /= 2)) then
      Shift(kInter-1,iter) = Delta
      Shift(kInter-2,iter) = Delta
    end if
    ! +d,+d
    !write(u6,*) ' +d,+d'
    Shift(kInter,iter) = Shift(kInter,iter)+Delta
    Shift(lInter,iter) = Shift(lInter,iter)+Delta
  else if (mod(kCount,4) == 2) then
    ! -d,+d
    !write(u6,*) ' -d,+d'
    Shift(kInter,iter) = -Two*Delta
    Shift(lInter,iter) = Zero
  else if (mod(kCount,4) == 3) then
    ! +d,-d
    !write(u6,*) ' +d,-d'
    Shift(kInter,iter) = Two*Delta
    Shift(lInter,iter) = -Two*Delta
  else if (mod(kCount,4) == 0) then
    ! -d,-d
    !write(u6,*) ' -d,-d'
    Shift(kInter,iter) = -Two*Delta
    Shift(lInter,iter) = Zero
  end if
end if

! Compute the new parameter set.
qInt(:,iter+1) = qInt(:,iter)+Shift(:,iter)

# ifdef _DEBUGPRINT_
call RecPrt('  qInt',' ',qInt,nInter,iter+1)
call RecPrt(' Shift',' ',Shift,nInter,iter)
# endif

end subroutine NwShft
