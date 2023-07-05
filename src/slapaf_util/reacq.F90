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
! Copyright (C) 2009, Roland Lindh                                     *
!***********************************************************************

subroutine ReacQ(V_X,nX,V_Q,nQ)
!***********************************************************************
!                                                                      *
!     Objective: Transform the "reaction vector" from cartesian        *
!                coordinates to Internal.                              *
!                                                                      *
!     Roland Lindh, Dept. of Theor. Chem., Lund University, Sweden     *
!     2009                                                             *
!***********************************************************************

use Slapaf_Info, only: BMx, Degen
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nX, nQ
real(kind=wp) :: V_Q(nQ), V_X(nX)
integer(kind=iwp) :: M, N, NRHS

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('BMx',' ',BMx,nX,nQ)
call RecPrt('V_X',' ',V_X,nX,1)
#endif

M = nX
N = nQ
NRHS = 1
call Eq_Solver('T',M,N,NRHS,BMx,.false.,Degen,V_X,V_Q)

#ifdef _DEBUGPRINT_
call RecPrt('V_Q',' ',V_Q,nQ,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ReacQ
