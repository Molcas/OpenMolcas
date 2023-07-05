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

use Slapaf_Info, only: Shift, qInt
use Slapaf_parameters, only: iter, Delta

integer nInter

nInter = size(Shift,1)
call NwShft_Internal(Shift,nInter,Iter,Delta,qInt)

contains

subroutine NwShft_Internal(dq,nInter,nIter,Delta,q)
  implicit real*8(A-H,O-Z)
  real*8 dq(nInter,nIter), q(nInter,nIter+1)
# include "real.fh"

# ifdef _DEBUGPRINT_
  call RecPrt('NwShft  q',' ',q,nInter,nIter)
  call RecPrt('NwShft dq',' ',dq,nInter,nIter-1)
# endif

  ! Compute the new shift

  !write(6,*) ' nIter=',nIter
  if (nIter < 2*nInter+1) then

    ! Shifts for the numerical Hessian

    jInter = (nIter+1)/2
    call dcopy_(nInter,[Zero],0,dq(1,nIter),1)
    if (mod(nIter,2) == 0) then
      dq(jInter,nIter) = -Two*Delta
    else
      ! Undo previous displacement
      if (jInter > 1) dq(jInter-1,nIter) = Delta
      dq(jInter,nIter) = Delta
    end if

  else

    ! Shifts for the numerical cubic force constants

    iCount = (nIter-2*nInter+3)/4
    !write(6,*) ' iCount=',iCount
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
    !write(6,*) 'kInter, lInter=',kInter,lInter
    kCount = nIter-2*nInter
    call dcopy_(nInter,[Zero],0,dq(1,nIter),1)
    ! Undo last change for numerical Hessian
    if (iCount == 1) dq(nInter,nIter) = Delta
    if (mod(kCount,4) == 1) then
      ! Undo change due to previous pair
      if (lInter /= 1) then
        dq(kInter,nIter) = Delta
        dq(lInter-1,nIter) = Delta
      else if ((lInter == 1) .and. (kInter /= 2)) then
        dq(kInter-1,nIter) = Delta
        dq(kInter-2,nIter) = Delta
      end if
      ! +d,+d
      !write(6,*) ' +d,+d'
      dq(kInter,nIter) = dq(kInter,nIter)+Delta
      dq(lInter,nIter) = dq(lInter,nIter)+Delta
    else if (mod(kCount,4) == 2) then
      ! -d,+d
      !write(6,*) ' -d,+d'
      dq(kInter,nIter) = -Two*Delta
      dq(lInter,nIter) = Zero
    else if (mod(kCount,4) == 3) then
      ! +d,-d
      !write(6,*) ' +d,-d'
      dq(kInter,nIter) = Two*Delta
      dq(lInter,nIter) = -Two*Delta
    else if (mod(kCount,4) == 0) then
      ! -d,-d
      !write(6,*) ' -d,-d'
      dq(kInter,nIter) = -Two*Delta
      dq(lInter,nIter) = Zero
    end if
  end if

  ! Compute the new parameter set.
  call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
  call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)

# ifdef _DEBUGPRINT_
  call RecPrt('  q',' ',q,nInter,nIter+1)
  call RecPrt(' dq',' ',dq,nInter,nIter)
# endif
end subroutine NwShft_Internal

end subroutine NwShft
