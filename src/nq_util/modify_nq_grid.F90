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

subroutine Modify_NQ_grid

use Grid_On_Disk
use nq_Info

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
parameter(L_Quad_Low=23,Threshold_High=1.0D-7,nR_Low=50)

!                                                                      *
!***********************************************************************
!                                                                      *
! Reduce the size and the accuracy of the grid temporarily

L_Quad_Save = L_Quad
Threshold_save = Threshold
nR_Save = nR
ThrC = Crowding

L_Quad = min(L_Quad,L_Quad_Low)
if (Quadrature(1:3) == 'LMG') then
  Threshold = max(Threshold_High,Threshold)
else
  nR = min(nR_Low,nR)
end if
Crowding = max(ThrC-Two,One)

write(6,*)
write(6,*) 'Modify the NQ grid!'
write(6,*)
call Funi_Print()
!                                                                      *
!***********************************************************************
!                                                                      *
! Change the Grid set index.

iGrid_Set = Intermediate
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Modify_NQ_grid
