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

subroutine Ttotal(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDim
real(kind=wp) :: T1(nDim*nDim), T2(nDim*nDim), T3(nDim*nDim), T4(nDim,nDim), Ttot(nDim,nDim), Ttot_Inv(nDim,nDim)
integer(kind=iwp) :: ipTemp, ipTemp2
#include "WrkSpc.fh"

call Allocate_Work(ipTemp,nDim**2)
call Allocate_Work(ipTemp2,nDim**2)
call Ttotal_(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim,Work(ipTemp),Work(ipTemp2))
call Free_Work(ipTemp2)
call Free_Work(ipTemp)

return

end subroutine Ttotal
