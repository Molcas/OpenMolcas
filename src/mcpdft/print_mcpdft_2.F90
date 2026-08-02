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
!   2018 S Dong - added print outs related to
!   2026 Bhavnesh Jangid - Added print outs for LPDFT quadrature results*
!***********************************************************************

subroutine print_MCPDFT_2(e_nuc,e_cas,e_ot,state,e_mcscf)

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Tau_a1, Tau_a2, Tau_b1, Tau_b2
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

call Add_Info('dens_tt',[Dens_I],1,6)
call Add_Info('dens_a1',[Dens_a1],1,6)
call Add_Info('dens_b1',[Dens_b1],1,6)
call Add_Info('dens_a2',[Dens_a2],1,6)
call Add_Info('dens_b2',[Dens_b2],1,6)
call Add_Info('exch_f',[CoefX],1,6)
call Add_Info('corr_f',[CoefR],1,6)
call Add_Info('excha_a',[Funcaa],1,6)
call Add_Info('excha_b',[Funcbb],1,6)
call Add_Info('corr_e',[Funccc],1,6)
call Add_Info('CASDFTE',[e_state],1,8)

end subroutine print_MCPDFT_2

subroutine print_lpdft_finalize(e_nuc,e_lpdft_ref,e_ot,e_mcscf,nStates)

use Constants, only: One, Zero
use Definitions, only: wp, iwp, u6
use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use mcpdft_input, only: mcpdft_options
use mcpdft_output, only: iprglb
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Tau_a1, Tau_a2, Tau_b1, Tau_b2
use PrintLevel, only: USUAL

implicit none
integer(kind=iwp), intent(in) :: nStates
real(kind=wp), intent(in) :: e_nuc, e_lpdft_ref, e_ot, e_mcscf(nStates)
real(kind=wp) :: e_lpdft, e_ot_1, Funcaa1, Funcbb1, Funccc1
integer(kind=iwp) :: iState
character(len=120) :: Line
character(len=6) :: Fmt2

e_lpdft = e_lpdft_ref+e_ot

if (iprglb >= USUAL) then
  write(Fmt2,'(A,I3.3,A)') '(6X,'
  Line = ''
  write(Line(4:),'(A)') 'LPDFT QUADRATURE RESULTS'
  write(u6,*)
  call CollapseOutput(1,Line)
  write(u6,Fmt2//'A)') repeat('-',len_trim(Line)-3)
  write(u6,*)

  if (nStates == 1) then
    write(u6,'(6X,A,40X,F18.8)') 'MCSCF reference energy',e_mcscf(1)
  else
    write(u6,'(6X,A)') 'MCSCF reference energies'
    do iState=1,nStates
      write(u6,'(6X,A,I4,T50,F18.8)') 'State',iState,e_mcscf(iState)
    end do
  end if
  write(u6,*)
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

  write(u6,'(6X,A,T70,F18.8)') 'Nuclear Repulsion energy',e_nuc
  write(u6,'(6X,A,T70,F18.8)') 'On-top energy (zeroth-order densities)',e_ot
  write(u6,*)
  if (mcpdft_options%otfnal%is_hybrid()) then
    write(u6,'(6X,A)') 'Information for hybrid PDFT:'
    write(u6,'(6X,A,T70,F6.2)') 'Wave function percentage (Lambda*100)', &
                                mcpdft_options%otfnal%lambda*1.0e2_wp
    do iState=1,nStates
      write(u6,'(6X,A,I4,T70,F18.8)') 'Wave function energy, state',iState, &
                                      mcpdft_options%otfnal%lambda*e_mcscf(iState)
    end do
    write(u6,*)
  end if

  if ((CoefX*CoefR /= Zero) .and. (CoefX /= One .or. CoefR /= One)) then
    Funcaa1 = Funcaa/CoefX
    Funcbb1 = Funcbb/CoefX
    Funccc1 = Funccc/CoefR
    e_ot_1 = e_ot-Funcaa-Funcbb-Funccc+Funcaa1+Funcbb1+Funccc1
    write(u6,*)
    write(u6,*)
    write(u6,'(6X,A,19X,F18.6)') 'Integrated alpha exchange energy (unscaled)',Funcaa1
    write(u6,'(6X,A,19X,F18.6)') 'Integrated beta  exchange energy (unscaled)',Funcbb1
    write(u6,'(6X,A,19X,F18.6)') 'Integrated  correlation   energy (unscaled)',Funccc1
    write(u6,'(6X,A,T70,F18.8)') 'On-top energy (zeroth-order densities, unscaled)',e_ot_1
  end if

  call CollapseOutput(0,Line)
  write(u6,*)
end if

call Add_Info('dens_tt',[Dens_I],1,6)
call Add_Info('dens_a1',[Dens_a1],1,6)
call Add_Info('dens_b1',[Dens_b1],1,6)
call Add_Info('dens_a2',[Dens_a2],1,6)
call Add_Info('dens_b2',[Dens_b2],1,6)
call Add_Info('exch_f',[CoefX],1,6)
call Add_Info('corr_f',[CoefR],1,6)
call Add_Info('excha_a',[Funcaa],1,6)
call Add_Info('excha_b',[Funcbb],1,6)
call Add_Info('corr_e',[Funccc],1,6)
call Add_Info('LPDFTE',[e_lpdft],1,8)

end subroutine print_lpdft_finalize
