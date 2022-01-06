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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetU_ER(U,R,n)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: compute U = R*[R^T*R]^(-1/2).
!
! (used by ER orbital localisation - hence the _ER)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(out) :: U(n,n)
real(kind=wp), intent(in) :: R(n,n)
#include "WrkSpc.fh"
integer(kind=iwp) :: ipISqrt, ipRTR, ipScr, ipSqrt, iTask, lISqrt, lRTR, lScr, lSqrt, n2, nn

if (n < 1) return

! Allocations.
! ------------

nn = n*(n+1)/2
n2 = n**2

lRTR = n2
lSqrt = n2
lISqrt = n2
lScr = 2*n2+nn
call GetMem('RTR','Allo','Real',ipRTR,lRTR)
call GetMem('Sqrt','Allo','Real',ipSqrt,lSqrt)
call GetMem('ISqrt','Allo','Real',ipISqrt,lISqrt)
call GetMem('Scr','Allo','Real',ipScr,lScr)

! Compute R^T*R.
! --------------

call DGEMM_('T','N',n,n,n,One,R,n,R,n,Zero,Work(ipRTR),n)

! Compute inverse square root of R^T*R.
! -------------------------------------

iTask = 2 ! compute sqrt as well as inverse sqrt
call SqrtMt(Work(ipRTR),n,iTask,Work(ipSqrt),Work(ipISqrt),Work(ipScr))

! Compute U.
! ----------

call DGEMM_('N','N',n,n,n,One,R,n,Work(ipISqrt),n,Zero,U,n)

! De-allocations.
! ---------------

call GetMem('Scr','Free','Real',ipScr,lScr)
call GetMem('ISqrt','Free','Real',ipISqrt,lISqrt)
call GetMem('Sqrt','Free','Real',ipSqrt,lSqrt)
call GetMem('RTR','Free','Real',ipRTR,lRTR)

end subroutine GetU_ER
