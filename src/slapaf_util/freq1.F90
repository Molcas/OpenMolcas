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

subroutine Freq1()
!***********************************************************************
!                                                                      *
! Object: Displacements for Numerical estimation of single rows and    *
!         columns of Hessian                                           *
!                                                                      *
! Called from: RlxCtl when Allocated(mRowH)=.True.                     *
!                                                                      *
! Author: Giovanni Ghigo, University of Torino, Italy                  *
!***********************************************************************

use Slapaf_Info, only: qInt
use Slapaf_parameters, only: Delta, iter

nInter = size(qInt,1)
call Freq1_Internal(iter,nInter,Delta/2.5d0,qInt)

contains

subroutine Freq1_Internal(nIter,nInter,Delta,qInt)

  use Slapaf_Info, only: Shift, mRowH

  implicit real*8(a-h,o-z)
# include "real.fh"
# include "print.fh"
  real*8 qInt(nInter,nIter+1)

  iRout = 183
  iPrint = nPrint(iRout)

  if (iPrint >= 99) then
    write(6,*) ' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
    call RecPrt('Initial Shift:','(10F9.6)',Shift,nInter,nIter)
    call RecPrt('Initial qInt:','(10F9.6)',qInt,nInter,nIter+1)
  end if

  ! Compute the new shift

  call dcopy_(nInter,[Zero],0,Shift(1,nIter),1)
  nRowH = 0
  if (allocated(mRowH)) nRowH = size(mRowH)
  if (nIter <= nRowH) then
    kInter = mRowH(nIter)
    Shift(kInter,nIter) = Delta
  end if
  if (nIter > 1) then
    jInter = mRowH(nIter-1)
    Shift(jInter,nIter) = -Delta ! Undo previous displacement
  end if

  ! Compute the new parameter set.

  call dcopy_(nInter,qInt(1,nIter),1,qInt(1,nIter+1),1)
  call DaXpY_(nInter,One,Shift(1,nIter),1,qInt(1,nIter+1),1)

  if (iPrint > 5) then
    write(6,*) ' Accumulate the gradient for yet one parameter set'
    write(6,*)
  end if

  if (iPrint >= 98) then
    write(6,*) ' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
    call RecPrt('Final Shift:','(10F9.6)',Shift,nInter,nIter)
    call RecPrt('Final  q:','(10F9.6)',qInt,nInter,nIter+1)
  end if

end subroutine Freq1_Internal

end subroutine Freq1
