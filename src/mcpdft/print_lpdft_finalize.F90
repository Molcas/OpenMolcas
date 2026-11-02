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

subroutine print_lpdft_finalize(e_nuc,e_lpdft_ref,e_ot,e_mcscf,nStates)

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
use mcpdft_input, only: mcpdft_options
use mcpdft_output, only: iprglb
use PrintLevel, only: USUAL
use Constants, only: One, Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nStates
real(kind=wp), intent(in) :: e_nuc, e_lpdft_ref, e_ot, e_mcscf(nStates)
integer(kind=iwp) :: iState
real(kind=wp) :: e_lpdft, e_ot_1, Funcaa1, Funcbb1, Funccc1
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
  call print_pdft_quadrature_info()

  write(u6,'(6X,A,T70,F18.8)') 'Nuclear Repulsion energy',e_nuc
  write(u6,'(6X,A,T70,F18.8)') 'On-top energy (zeroth-order densities)',e_ot
  write(u6,*)
  if (mcpdft_options%otfnal%is_hybrid()) then
    write(u6,'(6X,A)') 'Information for hybrid PDFT:'
    write(u6,'(6X,A,T70,F6.2)') 'Wave function percentage (Lambda*100)',mcpdft_options%otfnal%lambda*1.0e2_wp
    do iState=1,nStates
      write(u6,'(6X,A,I4,T70,F18.8)') 'Wave function energy, state',iState,mcpdft_options%otfnal%lambda*e_mcscf(iState)
    end do
    write(u6,*)
  end if

  if ((CoefX*CoefR /= Zero) .and. ((CoefX /= One) .or. (CoefR /= One))) then
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

call add_pdft_quadrature_info()
call Add_Info('LPDFTREF',[e_lpdft],1,8)

end subroutine print_lpdft_finalize
