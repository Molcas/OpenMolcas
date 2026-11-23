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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine InpPri_m()

use CI_Info, only: CI_St
use Functionals, only: Print_Info
use libxc_parameters, only: Clear_External_Params
use KSDFT_Info, only: CoefR, CoefX
use PrintLevel, only: SILENT, TERSE, USUAL, VERBOSE
use mcpdft_output, only: iPrLoc
use Fock_util_global, only: docholesky
use mcpdft_input, only: mcpdft_options
use rasscf_global, only: header, NAC, NFR, NIN, NROOTS, NSEC
use Orbital_Spaces, only: nash, nbas, ndel, nfro, nish, nssh, nsym
use CI_Input, only: nelec3, nhole1, nrs1, nrs2, nrs3
use rasscf_global, only: nconf
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: i, iPrLev, iSym, left, lPaper, k_rf
character(len=120) :: Line
character(len=8) :: Fmt1, Fmt2
character(len=3) :: lIrrep(8)

IPRLEV = IPRLOC(1)

! This should not be done in this function
call Put_dScalar('DFT exch coeff',CoefX)
call Put_dScalar('DFT corr coeff',CoefR)

if (mcpdft_options%extparam) then
  call CheckFuncParam(mcpdft_options%extparamfile)
else
  call Clear_External_Params()
end if

if (IPRLEV == SILENT) return

! Define the paper width
lPaper = 132

left = (lPaper-len(line))/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

if (iPrLev >= VERBOSE) then
  ! Print the ONEINT file identifier
  write(u6,*)
  write(u6,Fmt1) 'Header of the ONEINT file:'
  write(u6,Fmt1) '--------------------------'
  write(Line,'(36A2)') (Header(i),i=1,36)
  write(u6,Fmt1) trim(adjustl(Line))
  write(Line,'(36A2)') (Header(i),i=37,72)
  write(u6,Fmt1) trim(adjustl(Line))
  write(u6,*)

  ! Print cartesian coordinates of the system
  call PrCoor()

end if

if (iPrLev >= USUAL) then
  ! Print orbital and wavefunction specifications
  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Wave function specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt1) '-----------------------------'
  write(u6,*)
  if (NFR > 0) write(u6,Fmt2//'A,T45,I6)') 'Number of frozen shell electrons',2*NFR
  write(u6,Fmt2//'A,T45,I6)') 'Number of closed shell electrons',2*NIN
  write(u6,Fmt2//'A,T45,I6)') 'Number of electrons in active shells',CI_St(iState)%nActEl
  write(u6,Fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
  write(u6,Fmt2//'A,T45,I6)') 'Max nr of electrons in RAS3 space',NELEC3

  if (NFR > 0) write(u6,Fmt2//'A,T45,I6)') 'Number of frozen orbitals',NFR
  write(u6,Fmt2//'A,T45,I6)') 'Number of inactive orbitals',NIN
  write(u6,Fmt2//'A,T45,I6)') 'Number of active orbitals',NAC
  write(u6,Fmt2//'A,T45,I6)') 'Number of secondary orbitals',NSEC
  write(u6,Fmt2//'A,T45,F6.1)') 'Spin quantum number',Half*real(CI_St(iState)%iSpin-1,kind=wp)
  write(u6,Fmt2//'A,T45,I6)') 'State symmetry',CI_St(iState)%STSym
  write(u6,fmt2//'A,T40,I11)') 'Number of CSFs',NCONF
  write(u6,Fmt2//'A,T45,I6)') 'Number of RASSCF root(s) available',nroots
  call CollapseOutput(0,'Wave function specifications:')

  call Get_cArray('Irreps',lIrrep,24)
  lIrrep(1:nSym) = adjustr(lIrrep(1:nSym))

  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Orbital specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt1) '-----------------------'
  write(u6,*)
  write(u6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',(nRs1(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',(nRs2(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',(nRs3(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  call CollapseOutput(0,'Orbital specifications:')

end if

if (iPrLev >= TERSE) then
  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'MCPDFT specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt1) '----------------------'
  write(u6,*)
  if (DoCholesky) then
    write(u6,Fmt2//'A,T50,A)') 'Cholesky decomposition','On'
  else
    write(u6,Fmt2//'A,T50,A)') 'Cholesky decomposition','Off'
  end if
  if (mcpdft_options%lpdft) then
    write(u6,Fmt2//'A,T50,A)') 'Type of calculation','L-PDFT'
  else if (mcpdft_options%mspdft) then
    write(u6,Fmt2//'A,T50,A)') 'Type of calculation','MS-PDFT'
  else
    write(u6,Fmt2//'A,T50,A)') 'Type of calculation','MC-PDFT'
  end if
  write(u6,Fmt2//'A,T50,A)') 'On-Top Functional',trim(mcpdft_options%otfnal%otxc)
  write(u6,Fmt2//'A,T45,F9.2)') 'Exchange scaling factor',CoefX
  write(u6,Fmt2//'A,T45,F9.2)') 'Correlation scaling factor',CoefR
  write(u6,Fmt2//'A,T45,F9.2)') 'Wave function energy weight',mcpdft_options%otfnal%lambda
  if (mcpdft_options%wjob) write(u6,Fmt1) 'Final energies (and CI vectors) will be written to wave function file'
  if (mcpdft_options%grad) then
    write(u6,Fmt1) 'On-top potentials are computed'
    if (mcpdft_options%nac) then
      write(u6,fmt2//'A,T45,I6,1X,"/",1X,I6)') 'MSPDFT states for NAC',mcpdft_options%nac_states(1),mcpdft_options%nac_states(2)
    else
      write(u6,Fmt2//'A,T45,I6)') 'Root chosen for geometry opt.',mcpdft_options%rlxroot
    end if
  end if
  call CollapseOutput(0,'MCPDFT specifications:')

  ! Print out grid information
  call Funi_Print()

  ! Print reaction field information
  if (mcpdft_options%rfpert .or. mcpdft_options%rf_on) then
    write(u6,*)
    write(u6,fmt1) 'Reaction field specifications:'
    write(u6,fmt1) '------------------------------'
    if (mcpdft_options%rf_scf) then
      if (mcpdft_options%noneq) then
        write(u6,Fmt2//'A)') 'Calculation type: non-equilibrium (slow component from previous equilibrium LPDFT)'
      else
        write(u6,Fmt2//'A)') 'Calculation type: equilibrium'
      end if
      k_rf = mcpdft_options%rfroot
      if (k_rf <= 0) then
        write(u6,Fmt2//'A,T44,A)') 'Reaction field from state:',' State-Averaged'
      else
        write(u6,Fmt2//'A,T44,I2)') 'Reaction field from state:',k_rf
      end if
    else if (mcpdft_options%lpdft) then
      write(u6,'(6X,A)') 'The reaction field of a previous calculation is added to the one-electron part'
      write(u6,'(6X,A)') 'of the L-PDFT Hamiltonian. It is not reevaluated.'
    else
      write(u6,'(6X,A)') 'The reaction field of a previous calculation is used. It is not'
      write(u6,'(6X,A)') 'reevaluated as MC-PDFT has no Hamiltonian.'
    end if
    write(u6,*)
  end if

end if

if (iPrLev >= USUAL) then
  ! Print our DFT functional specifications
  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'DFT functional specifications:'
  call CollapseOutput(1,Line)
  write(u6,Fmt1) '------------------------------'
  call libxc_version()
  call Print_Info()
  call CollapseOutput(0,'DFT functional specifications:')
end if

end subroutine InpPri_m
