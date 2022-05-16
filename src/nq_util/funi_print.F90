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
use nq_Info, only: Angular_Pruning, Crowding, Fade, iOpt_Angular, L_Quad, NQ_Direct, nR, On, Quadrature, T_Y, Threshold
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iOpt, iPrint
real(kind=wp) :: EThr
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
iPrint = iPrintLevel(-1)
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_dScalar('EThr',EThr)
T_Y = min(T_Y,EThr*0.1_wp)
Threshold = min(Threshold,EThr*1.0e-4_wp)
!                                                                      *
!***********************************************************************
!                                                                      *
if ((.not. Reduce_Prt()) .and. (iPrint >= 2)) then
  write(u6,*)
  write(u6,'(6X,A)') 'Numerical integration parameters'
  write(u6,'(6X,A)') '--------------------------------'
  write(u6,'(6X,A,21X,A)') 'Radial quadrature type:    ',Quadrature

  if (Quadrature(1:3) == 'LMG') then
    write(u6,'(6X,A,E11.4)') 'Radial quadrature accuracy:',Threshold
  else
    write(u6,'(6X,A,18X,I5)') 'Size of radial grid:       ',nR
  end if

  if (btest(iOpt_Angular,2)) then
    write(u6,'(6X,A,25X,I4)') 'Lebedev angular grid:',L_Quad
  else if (btest(iOpt_Angular,0)) then
    write(u6,'(6X,A,I4)') 'Lobatto angular grid, l_max:',L_Quad
  else
    write(u6,'(6X,A,I4)') 'Gauss and Gauss-Legendre angular grid, l_max:',L_Quad
  end if

  if (Angular_Pruning == On) then
    write(u6,'(6X,A,1X,ES9.2)') 'Angular grid prunned with the crowding factor:',Crowding
    write(u6,'(6X,A,1X,ES9.2)') '                            and fading factor:',Fade
  end if
  if (btest(iOpt_Angular,1)) then
    write(u6,'(6X,A)') 'The whole atomic grid is scanned for each sub block.'
  end if

  write(u6,'(6X,A,2X,ES9.2)') 'Screening threshold for integral computation:',T_Y
  if (Quadrature(1:3) /= 'LMG') then
    write(u6,'(6X,A,20X,ES9.2)') 'Radial quadrature accuracy:',Threshold
  end if

  write(u6,'(6X,A,17X,I7)') 'Maximum batch size:        ',nGridMax
  if (NQ_Direct == On) then
    write(u6,'(6X,A)') 'AO values are recomputed each iteration'
  else
    write(u6,'(6X,A)') 'AO values are stored on disk'
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put flag on RUNFILE to indicate that we are doing DFT.

!call Get_iOption(iOpt)
call Get_iScalar('System BitSwitch',iOpt)
iOpt = ibset(iOpt,6)
!call Put_iOption(iOpt)
call Put_iScalar('System BitSwitch',iOpt)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Funi_Print
