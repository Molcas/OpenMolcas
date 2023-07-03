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

subroutine ReacX(V_Q,nQ,V_X,nX)
!***********************************************************************
!                                                                      *
!     Objective: Transform the "reaction vector" from internal         *
!                coordinates to Cartesians.                            *
!                                                                      *
!     Roland Lindh, Dept. of Theor. Chem., Lund University, Sweden     *
!     2009                                                             *
!***********************************************************************

use Slapaf_Info, only: BMx

implicit real*8(a-h,o-z)
real*8 V_Q(nQ), V_X(nX)
!                                                                      *
!***********************************************************************
!                                                                      *
call DGEMM_('N','N',nX,1,nQ,1.0d0,BMx,nX,V_Q,nQ,0.0d0,V_X,nX)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ReacX
