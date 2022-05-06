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

subroutine Reset_NQ_Grid()

use Grid_On_Disk
use nq_Info

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "itmax.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Reset the size and the accuracy of the grid to the requested values.

L_Quad = L_Quad_Save
if (Quadrature(1:3) == 'LMG') then
  Threshold = Threshold_Save
else
  nR = nR_save
end if

Crowding = ThrC

write(6,*)
write(6,'(6X,A)') 'Reset the NQ grid!'
write(6,*)
call Funi_Print()
write(6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
! Change the Grid set index

iGrid_Set = final
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Reset_NQ_Grid
