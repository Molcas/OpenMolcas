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
! Copyright (C) Giovanni Li Manni                                      *
!***********************************************************************

subroutine print_MCPDFT(CASDFT_E)
!*****************************************************************
! Purpose:
! This routine is called from RASSCF when performing MC-PDFT jobs,
! for printing out some densities information and some functionals
! as computed in DRVNQ and children routines.
!
! Author:
! G. Li Manni (GLM)
!*****************************************************************

use KSDFT_Info, only: Funcaa, Funcbb, Funccc
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC

implicit none
real*8 CASDFT_E
integer left
character(len=6) Fmt2
character(len=120) Line
integer IPRLEV

IPRLEV = IPRLOC(6)
if (IPRLEV >= USUAL) then

  left = 6
  write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
  Line = ''
  write(Line(left-2:),'(A)') 'MC-PDFT run print out'
  write(6,'(6X,80A)')
  call CollapseOutput(1,Line)
  write(6,Fmt2//'A)') repeat('-',len_trim(Line)-3)
  write(6,'(6X,80A)')

  write(6,'(6X,A25,45X,F10.3)') 'Integrated total density:',Dens_I
  write(6,'(6X,A58,12X,F10.3)') 'Integrated alpha density before functional transformation:',Dens_a1
  write(6,'(6X,A58,12X,F10.3)') 'Integrated  beta density before functional transformation:',Dens_b1
  write(6,'(6X,A58,12X,F10.3)') 'Integrated alpha density  after functional transformation:',Dens_a2
  write(6,'(6X,A58,12X,F10.3)') 'Integrated  beta density  after functional transformation:',Dens_b2
  write(6,'(6X,80A)')
  write(6,'(6X,A32,30X,F18.6)') 'Integrated alpha exchange energy',Funcaa
  write(6,'(6X,A32,30X,F18.6)') 'Integrated beta  exchange energy',Funcbb
  write(6,'(6X,A32,30X,F18.6)') 'Integrated  correlation   energy',Funccc
  write(6,'(6X,80A)')
  write(6,'(6X,A20,42X,F18.8)') 'Total CAS-DFT energy',CASDFT_E

  call CollapseOutput(0,Line)
  write(6,'(6X,80A)')
end if

call Add_Info('dens_tt',[Dens_I],1,6)
call Add_Info('dens_a1',[Dens_a1],1,6)
call Add_Info('dens_b1',[Dens_b1],1,6)
call Add_Info('dens_a2',[Dens_a2],1,6)
call Add_Info('dens_b2',[Dens_b2],1,6)
call Add_Info('excha_a',[Funcaa],1,6)
call Add_Info('excha_b',[Funcbb],1,6)
call Add_Info('corr_e',[Funccc],1,6)
call Add_Info('CASDFTE',[CASDFT_E],1,8)

end subroutine print_MCPDFT
