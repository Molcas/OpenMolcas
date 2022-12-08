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

subroutine Funi_Init()

use nq_Info, only: Angular_Pruning, Crowding, Fade, Grid_Type, iOpt_Angular, L_Quad, MBC, Moving_Grid, NQ_Direct, nR, ntotgp, Off, &
                   On, Packing, Quadrature, Rotational_Invariance, T_Y, Threshold, R_Max
use Constants, only: Zero, Three, Six
use Definitions, only: wp, iwp

implicit none

!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize the default values of the parameters.

! Default Grid
Quadrature = 'MHL '
nR = 75
L_Quad = 29
Crowding = Three
Fade = Six
MBC = ' '

ntotgp = 0

! Various default thresholds for the integral evaluation.

T_Y = 1.0e-11_wp
Threshold = 1.0e-25_wp

Angular_Pruning = On
Grid_Type = Moving_Grid
Rotational_Invariance = On
NQ_Direct = Off
!NQ_Direct = On
!Packing = On
Packing = Off

! Bit 0: set Lobatto, not set Gauss and Gauss-Legendre
! Bit 1: set scan the whole atomic grid, not set use subset
! Bit 2: set Lebedev, override bit 0
iOpt_Angular = ibset(0_iwp,2)
!                                                                      *
!***********************************************************************
!                                                                      *
R_Max(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Funi_Init
