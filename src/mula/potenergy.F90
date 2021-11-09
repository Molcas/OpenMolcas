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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine PotEnergy(A,nMat,iCre,iAnn,energy,grad,Hess,D3,D4,max_term,W,max_ord,nosc,nOscOld)
!  Purpose:
!    Calculate matrix elements of potential energy terms.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Zero, One

!use TabMod
implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 A(0:mdim1,0:ndim1)
integer nMat(0:ndim1,ndim2)
integer iAnn(0:ndim1,ndim2)
integer iCre(0:ndim1,ndim2)
real*8 rdx(4)
real*8 grad(noscold)
real*8 Hess(noscold,noscold)
real*8 D3(noscold,noscold,noscold)
real*8 D4(noscold,noscold,noscold,noscold)
real*8 W(noscold,nosc)
#include "WrkSpc.fh"

! Zeroth order term.
call dcopy_(max_ord+1,[Energy],0,A,max_Ord+2)
rdx(1) = One
rdx(2) = One
rdx(3) = One
rdx(4) = One
call GetMem('Temp','Allo','Real',ipTemp,nOscOld**4)

! First order terms.
if (max_term > 0) then
  call GetMem('grad_2','Allo','Real',ipgrad_2,nOsc)
  call DGEMM_('T','N',1,nOsc,nOscOld,One,grad,nOscOld,W,nOscOld,Zero,Work(ipgrad_2),1)
  call Mul1(nMat,A,icre,iann,Work(ipgrad_2),max_ord,nosc,rdx)
  call GetMem('grad_2','Free','Real',ipgrad_2,nOsc)
end if

! Second order terms.
if (max_term > 1) then
  call GetMem('Hess_2','Allo','Real',ipHess_2,nOsc**2)
  call DGEMM_('T','N',nOscOld,nOsc,nOscOld,One,Hess,nOscOld,W,nOscOld,Zero,Work(ipTemp),nOscOld)
  call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipTemp),nOscOld,W,nOscOld,Zero,Work(ipHess_2),nOsc)
  call Mul2(nMat,A,icre,iann,Work(ipHess_2),max_ord,nosc,rdx)
  call GetMem('Hess_2','Free','Real',ipHess_2,nOsc**2)
end if

! Third order terms.
if (max_term > 2) then
  call GetMem('D3_2','Allo','Real',ipD3_2,nOsc**3)
  call DGEMM_('T','N',nOscOld**2,nOsc,nOscOld,One,D3,nOscOld,W,nOscOld,Zero,Work(ipTemp),nOscOld**2)
  call DGEMM_('T','N',nOsc*nOscOld,nOsc,nOscOld,One,Work(ipTemp),nOscOld,W,nOscOld,Zero,Work(ipD3_2),nOsc*nOscOld)
  call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,Work(ipD3_2),nOscOld,W,nOscOld,Zero,Work(ipTemp),nOsc**2)
  call dcopy_(nOsc**3,Work(ipTemp),1,Work(ipD3_2),1)
  call Mul3(nMat,A,icre,iann,Work(ipD3_2),max_ord,nosc,rdx)
  call GetMem('D3_2','Free','Real',ipD3_2,nOsc**3)
end if

! Fourth order terms.
if (max_term > 3) then
  call GetMem('D4_2','Allo','Real',ipD4_2,nOsc**4)

  call DGEMM_('T','N',nOscOld**3,nOsc,nOscOld,One,D4,nOscOld,W,nOscOld,Zero,Work(ipTemp),nOscOld**3)
  call DGEMM_('T','N',nOsc*nOscOld**2,nOsc,nOscOld,One,Work(ipTemp),nOscOld,W,nOscOld,Zero,Work(ipD4_2),nOsc*nOscOld**2)
  call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,Work(ipD4_2),nOscOld,W,nOscOld,Zero,Work(ipTemp),nOsc**2*nOscOld)
  call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Work(ipTemp),nOscOld,W,nOscOld,Zero,Work(ipD4_2),nOsc**3)
  call Mul4(nMat,A,icre,iann,Work(ipD4_2),max_ord,nosc,rdx)
  call GetMem('D4_2','Free','Real',ipD4_2,nOsc**4)
end if

call GetMem('Temp','Free','Real',ipTemp,nOscOld**4)

end subroutine PotEnergy
