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

subroutine Funi_Print()

use nq_Grid, only: nGridMax
use nq_Info
implicit real*8(A-H,O-Z)
logical Check
logical, external :: Reduce_Prt
! Statement function
Check(i,j) = iand(i,2**(j-1)) /= 0

!                                                                      *
!***********************************************************************
!                                                                      *
iPrint = iPrintLevel(-1)
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_dScalar('EThr',EThr)
T_Y = min(T_Y,EThr*1.0D-1)
Threshold = min(Threshold,EThr*1.0D-4)
!                                                                      *
!***********************************************************************
!                                                                      *
if ((.not. Reduce_Prt()) .and. (iPrint >= 2)) then
  write(6,*)
  write(6,'(6X,A)') 'Numerical integration parameters'
  write(6,'(6X,A)') '--------------------------------'
  write(6,'(6X,A,21X,A)') 'Radial quadrature type:    ',Quadrature

  if (Quadrature(1:3) == 'LMG') then
    write(6,'(6X,A,E11.4)') 'Radial quadrature accuracy:',Threshold
  else
    write(6,'(6X,A,18X,I5)') 'Size of radial grid:       ',nR
  end if

  if (Check(iOpt_Angular,3)) then
    write(6,'(6X,A,25X,I4)') 'Lebedev angular grid:',L_Quad
  else if (Check(iOpt_Angular,1)) then
    write(6,'(6X,A,I4)') 'Lobatto angular grid, l_max:',L_Quad
  else
    write(6,'(6X,A,I4)') 'Gauss and Gauss-Legendre angular grid, l_max:',L_Quad
  end if

  if (Angular_Prunning == On) then
    write(6,'(6X,A,1X,ES9.2)') 'Angular grid prunned with the crowding factor:',Crowding
    write(6,'(6X,A,1X,ES9.2)') '                            and fading factor:',Fade
  end if
  if (Check(iOpt_Angular,2)) then
    write(6,'(6X,A)') 'The whole atomic grid is scanned for each sub block.'
  end if

  write(6,'(6X,A,2X,ES9.2)') 'Screening threshold for integral computation:',T_Y
  if (Quadrature(1:3) /= 'LMG') then
    write(6,'(6X,A,20X,ES9.2)') 'Radial quadrature accuracy:',Threshold
  end if

  write(6,'(6X,A,17X,I7)') 'Maximum batch size:        ',nGridMax
  if (NQ_Direct == On) then
    write(6,'(6X,A)') 'AO values are recomputed each iteration'
  else
    write(6,'(6X,A)') 'AO values are stored on disk'
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put flag on RUNFILE to indicate that we are doing DFT.

!call Get_iOption(iOpt)
call Get_iScalar('System BitSwitch',iOpt)
iOpt = ior(iOpt,2**6)
!call Put_iOption(iOpt)
call Put_iScalar('System BitSwitch',iOpt)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Funi_Print
