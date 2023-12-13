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
integer(kind=iwp), intent(in) :: nX, nQ
real(kind=wp), intent(in) :: V_X(nX)
real(kind=wp), intent(out) :: V_Q(nQ)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('BMx',' ',BMx,nX,nQ)
call RecPrt('V_X',' ',V_X,nX,1)
#endif

call Eq_Solver('T',nX,nQ,1,BMx,.false.,Degen,V_X,V_Q)

#ifdef _DEBUGPRINT_
call RecPrt('V_Q',' ',V_Q,nQ,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ReacQ
