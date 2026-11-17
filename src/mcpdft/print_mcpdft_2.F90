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

subroutine print_MCPDFT_2(e_nuc,e_cas,e_ot,state,e_mcscf)

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use mspdft, only: mspdftmethod
use PrintLevel, only: USUAL
use mcpdft_output, only: iPrGlb
use mcpdft_input, only: mcpdft_options
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: e_nuc, e_cas, e_ot, e_mcscf
integer(kind=iwp), intent(in) :: state
real(kind=wp) :: CASDFT_E_1, E_ot_1, e_pdft, e_state, Funcaa1, Funcbb1, Funccc1
character(len=120) :: Line
character(len=6) :: Fmt2

e_state = e_cas+e_ot
e_pdft = e_state
if (mcpdft_options%otfnal%is_hybrid()) e_state = mcpdft_options%otfnal%lambda*e_mcscf+(One-mcpdft_options%otfnal%lambda)*e_state

if (iPrGlb >= USUAL) then
  write(Fmt2,'(A,I3.3,A)') '(6X,'
  Line = ''
  if (mcpdft_options%mspdft) then
    write(Line(4:),'(2A,1X,I4.4)') trim(MSpdftMethod),' INTERMEDIATE STATE',state
  else
    write(Line(4:),'(A,1X,I4.4)') 'MC-PDFT RESULTS, STATE',state
  end if
  write(u6,*)
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') repeat('-',len_trim(Line)-3)
  write(u6,*)

  write(u6,'(6X,A,40X,F18.8)') 'MCSCF reference energy',e_mcscf
  write(u6,*)
  call print_pdft_quadrature_info()

  write(u6,'(6X,A,38X,F18.8)') 'Nuclear Repulsion energy',E_nuc
  write(u6,'(6X,A,36X,F18.8)') 'CASSCF contribution energy',E_cas
  write(u6,'(6X,A,49X,F18.8)') 'On-top energy',E_ot
  write(u6,*)
  if (mcpdft_options%otfnal%is_hybrid()) then
    write(u6,'(6X,A)') 'Information for hybrid PDFT:'
    write(u6,'(6X,A,37X,F6.2)') 'Wave function percentage (Lambda*100)',mcpdft_options%otfnal%lambda*1.0e2_wp
    write(u6,'(6X,A,42X,F18.8)') 'Wave function energy',mcpdft_options%otfnal%lambda*e_mcscf
    write(u6,'(6X,A,51X,F18.8)') 'PDFT energy',(One-mcpdft_options%otfnal%lambda)*e_pdft
    write(u6,*)
  end if

  if ((CoefX*CoefR /= Zero) .and. (CoefX /= One .or. CoefR /= One)) then
    Funcaa1 = Funcaa/CoefX
    Funcbb1 = Funcbb/CoefX
    Funccc1 = Funccc/CoefR
    E_ot_1 = E_ot-Funcaa-Funcbb-Funccc+Funcaa1+Funcbb1+Funccc1
    CASDFT_E_1 = e_state-E_ot+E_ot_1
    write(u6,*)
    write(u6,*)
    write(u6,'(6X,A,19X,F18.6)') 'Integrated alpha exchange energy (unscaled)',Funcaa1
    write(u6,'(6X,A,19X,F18.6)') 'Integrated beta  exchange energy (unscaled)',Funcbb1
    write(u6,'(6X,A,19X,F18.6)') 'Integrated  correlation   energy (unscaled)',Funccc1
    write(u6,'(6X,A,38X,F18.8)') 'On-top energy (unscaled)',E_ot_1
    write(u6,'(6X,A,31X,F18.8)') 'Total MC-PDFT energy (unscaled)',CASDFT_E_1
  end if

  write(u6,'(6X,A,I3,29X,F18.8)') 'Total MCPDFT energy for state ',state,e_state

  call CollapseOutput(0,Line)
  write(u6,*)
end if

call add_pdft_quadrature_info()
call Add_Info('CASDFTE',[e_state],1,8)

end subroutine print_MCPDFT_2
