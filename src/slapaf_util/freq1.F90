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

use Slapaf_Info, only: Delta, iter, mRowH, qInt, Shift
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
# include "print.fh"
integer(kind=iwp) :: iPrint, iRout, jInter, kInter, nInter, nRowH
real(kind=wp) :: Delta_

nInter = size(qInt,1)
Delta_ = Delta/2.5_wp

iRout = 183
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' [Freq1] nInter, iter, Delta =',nInter,iter,Delta_
  call RecPrt('Initial Shift:','(10F9.6)',Shift,nInter,iter)
  call RecPrt('Initial qInt:','(10F9.6)',qInt,nInter,iter+1)
end if

! Compute the new shift

Shift(:,iter) = Zero
nRowH = 0
if (allocated(mRowH)) nRowH = size(mRowH)
if (iter <= nRowH) then
  kInter = mRowH(iter)
  Shift(kInter,iter) = Delta_
end if
if (iter > 1) then
  jInter = mRowH(iter-1)
  Shift(jInter,iter) = -Delta_ ! Undo previous displacement
end if

! Compute the new parameter set.

qInt(:,iter+1) = qInt(:,iter)+Shift(:,iter)

if (iPrint > 5) then
  write(u6,*) ' Accumulate the gradient for yet one parameter set'
  write(u6,*)
end if

if (iPrint >= 98) then
  write(u6,*) ' [Freq1] nInter, iter, Delta =',nInter,iter,Delta_
  call RecPrt('Final Shift:','(10F9.6)',Shift,nInter,iter)
  call RecPrt('Final  q:','(10F9.6)',qInt,nInter,iter+1)
end if

end subroutine Freq1
