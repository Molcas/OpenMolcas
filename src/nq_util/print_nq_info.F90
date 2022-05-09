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

subroutine Print_NQ_Info()

use nq_Info, only: Dens_I, Energy_integrated, Grad_I, nTotGP, Tau_I
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iPL
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPL >= 3) then
  call GAIGOP_SCAL(nTotGP,'+')
  write(u6,*)
  write(u6,'(6X,A,T52,F17.10)') 'Integrated DFT Energy   ',Energy_integrated
  write(u6,'(6X,A,T56,G17.10)') 'Integrated number of electrons',Dens_I
  if (Grad_I /= Zero) write(u6,'(6X,A,T56,G17.10)') 'Integrated |grad|             ',Grad_I
  if (Tau_I /= Zero) write(u6,'(6X,A,T56,G17.10)') 'Integrated tau                ',Tau_I
  write(u6,'(6X,A,T54,I13)') 'Total number of prunned grid points  ',nTotGP
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Add_Info('DFT_Energy',[Energy_integrated],1,6)
call Add_Info('NQ_Density',[Dens_I],1,8)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Print_NQ_Info
