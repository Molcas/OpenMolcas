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
! History:                                                             *
! 2018 S Dong - added print outs related to scaling                    *
! 2026 Bhavnesh Jangid - Added print outs for LPDFT quadrature results *
!***********************************************************************

subroutine print_pdft_quadrature_info()

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Tau_a1, Tau_a2, Tau_b1, Tau_b2
use Definitions, only: u6

implicit none

write(u6,'(6X,A,45X,F10.3)') 'Integrated total density:',Dens_I
write(u6,'(6X,A,12X,F10.3)') 'Integrated alpha density before functional transformation:',Dens_a1
write(u6,'(6X,A,12X,F10.3)') 'Integrated  beta density before functional transformation:',Dens_b1
write(u6,'(6X,A,12X,F10.3)') 'Integrated alpha density  after functional transformation:',Dens_a2
write(u6,'(6X,A,12X,F10.3)') 'Integrated  beta density  after functional transformation:',Dens_b2
write(u6,'(6X,A,4X,F18.3)') 'Integrated alpha tau     before functional transformation:',Tau_a1
write(u6,'(6X,A,4X,F18.3)') 'Integrated  beta tau     before functional transformation:',Tau_b1
write(u6,'(6X,A,4X,F18.3)') 'Integrated alpha tau      after functional transformation:',Tau_a2
write(u6,'(6X,A,4X,F18.3)') 'Integrated  beta tau      after functional transformation:',Tau_b2
write(u6,'(6X,A)') 'NOTE:'
write(u6,'(6X,2A)') 'Densities after transformation are ','intermediate quantities'
write(u6,'(6X,2A)') 'and should not be interpreted as ','real spin densities'
write(u6,*)
write(u6,'(6X,A,32X,F18.6)') 'Exchange energy scaling factor',CoefX
write(u6,'(6X,A,29X,F18.6)') 'Correlation energy scaling factor',CoefR
write(u6,*)
write(u6,'(6X,A,30X,F18.6)') 'Integrated alpha exchange energy',Funcaa
write(u6,'(6X,A,30X,F18.6)') 'Integrated beta  exchange energy',Funcbb
write(u6,'(6X,A,30X,F18.6)') 'Integrated  correlation   energy',Funccc
write(u6,*)

end subroutine print_pdft_quadrature_info
