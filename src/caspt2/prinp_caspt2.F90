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
! Copyright (C) 1994, Markus P. Fuelscher                              *
!               1994, Per Ake Malmqvist                                *
!***********************************************************************

subroutine prinp_caspt2()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     - echo the input parameters                                      *
!                                                                      *
!     calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P.-AA. Malmqvist                              *
!     University of Lund, Sweden, 1994                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use PrintLevel, only: TERSE, USUAL, VERBOSE
use caspt2_global, only: do_csf, do_grad, do_nac, imag_shift, ipea_shift, iPrGlb, iRoot1, iRoot2, real_shift, sigma_p_epsilon, &
                         sigma_p_exponent
use caspt2_module, only: DWType, Header, HZero, IfDOrtho, IfDW, IfMix, IfMSCoup, IfRMS, IfsadRef, IfXMS, iRlxRoot, iSCF, iSpin, &
                         mState, nActel, nAsh, nAshT, nBas, nConf, nDel, nEle3, nFro, nHole1, nIsh, nIshT, nRoots, nSsh, nSshT, &
                         nState, nSym, Orbin, PT2Method, RFPert, STSym, Zeta
use SC_NEVPT2, only: Do_FIC, Do_SC, SC_prop, SC_thres
#ifdef _DMRG_
use caspt2_global, only: compressMPS
use caspt2_module, only: DMRG
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
use caspt2_module, only: DoCumulant
#endif
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iSym, left, lLine, lPaper
character(len=120) :: Line
character(len=20) :: calctype, FockOpType
character(len=8) :: fmt1, fmt2
character(len=3) :: lIrrep(8)

!----------------------------------------------------------------------*
!     Start and define the paper width,                                *
!     initialize blank and header lines                                *
!----------------------------------------------------------------------*
Line = ' '
lLine = len(Line)
lPaper = 132
left = (lPaper-lLine)/2
write(fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
!     Print the ONEINT file identifier                                 *
!----------------------------------------------------------------------*
if (iprglb >= VERBOSE) then
  write(u6,*)
  write(u6,fmt1) 'Header of the ONEINT file:'
  write(u6,fmt1) '--------------------------'
  write(Line,'(36A2)') (Header(i),i=1,36)
  write(u6,fmt1) trim(adjustl(Line))
  write(Line,'(36A2)') (Header(i),i=37,72)
  write(u6,fmt1) trim(adjustl(Line))
  write(u6,*)
end if
!----------------------------------------------------------------------*
!     Print cartesian coordinates of the system                        *
!----------------------------------------------------------------------*
if (iprglb >= VERBOSE) call prCoor()
!----------------------------------------------------------------------*
!     Print orbital and wavefunction specifications                    *
!----------------------------------------------------------------------*
if (iprglb >= USUAL) then
  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Wave function specifications:'
  call CollapseOutput(1,Line)
  write(u6,fmt1) '-----------------------------'
  write(u6,*)
  write(u6,fmt2//'A,T45,I6)') 'Number of closed shell electrons',2*NISHT
  write(u6,fmt2//'A,T45,I6)') 'Number of electrons in active shells',NACTEL
  write(u6,fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
  write(u6,fmt2//'A,T45,I6)') 'Max number of electrons in RAS3 space',NELE3
  write(u6,fmt2//'A,T45,I6)') 'Number of inactive orbitals',NISHT
  write(u6,fmt2//'A,T45,I6)') 'Number of active orbitals',NASHT
  write(u6,fmt2//'A,T45,I6)') 'Number of secondary orbitals',NSSHT
  write(u6,fmt2//'A,T45,F6.1)') 'Spin quantum number',Half*real(ISPIN-1,kind=wp)
  write(u6,fmt2//'A,T45,I6)') 'State symmetry',STSYM
  write(u6,fmt2//'A,T40,I11)') 'Number of CSFs',NCONF
  write(u6,fmt2//'A,T45,I6)') 'Number of CASSCF root(s) available',NROOTS
  if (do_nac) then
    write(u6,fmt2//'A,T45,I6,1X,"/",1X,I6)') trim(PT2Method)//' states for NAC',iRoot1,iRoot2
  else
    write(u6,fmt2//'A,T45,I6)') trim(PT2Method)//' state passed to geometry opt.',iRlxRoot
  end if
  if (ifmix) write(u6,fmt2//'A,T45,10I3)') 'A file JOBMIX will be created'
  if (nstate > 1) then
    write(u6,fmt1) 'This is a MULTI-STATE CASSCF reference'
    write(u6,fmt2//'A,T45,I6)') 'Number of CI roots used',NSTATE
    write(u6,fmt2//'A,(T47,10I4))') 'These are:',(MSTATE(I),I=1,NSTATE)
    if (ifmscoup) then
      write(u6,fmt1) 'Off-diagonal elements of Heff are computed'
    else
      write(u6,fmt1) 'Heff is assumed diagonal'
    end if
  else
    if (iscf == 0) then
      write(u6,fmt1) 'This is a CASSCF or RASSCF reference function'
#     ifdef _ENABLE_BLOCK_DMRG_
      if (DoCumulant) write(u6,fmt1) 'Using 4-RDM cumulant approximation, activated by 3RDM keyword in RASSCF'
#     elif _ENABLE_CHEMPS2_DMRG_
      if (DoCumulant) write(u6,fmt1) 'This is a DMRG reference with exact 4-RDM, activated by 3RDM keyword in RASSCF'
#     elif _DMRG_
      if (DMRG) then
        write(u6,fmt1) 'This is a DMRG reference wave function, from QCMaquis'
        if (compressMPS > 0) write(u6,fmt1) 'Using compressed MPS for computing 4-RDM.'
      end if
#     endif
    else if (iscf == 1) then
      write(u6,fmt1) 'This is a closed shell RHF reference function'
    else
      write(u6,fmt1) 'This is a high spin open shell RHF reference function'
    end if
  end if
  call CollapseOutput(0,'Wave function specifications:')
end if

call Get_cArray('Irreps',lIrrep,24)
lIrrep(1:nSym) = adjustr(lIrrep(1:nSym))

if (iprglb >= USUAL) then
  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A)') 'Orbital specifications:'
  call CollapseOutput(1,Line)
  write(u6,fmt1) '-----------------------'
  write(u6,*)
  write(u6,fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(u6,fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
  write(u6,fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  call CollapseOutput(0,'Orbital specifications:')
end if
!----------------------------------------------------------------------*
!     Print routing information                                        *
!----------------------------------------------------------------------*
if (iprglb >= TERSE) then
  if (rfpert) then
    write(u6,*)
    write(u6,fmt1) 'Reaction field specifications:'
    write(u6,fmt1) '------------------------------'
    write(u6,*)
    write(u6,'(6X,A)') 'An external reaction field was determined previously and added to the one-electron Hamiltonian'
    write(u6,'(6X,A)') 'It will not be reevaluated even though dynamic correlation may change the density.'
    write(u6,*)
  end if

  write(u6,*)
  Line = ' '
  write(Line(left-2:),'(A,A)') trim(PT2Method),' specifications:'
  call CollapseOutput(1,Line)
  write(u6,fmt1) '----------------------'
  write(u6,*)

  if (IFMSCOUP) then
    if (IFDW) then
      FockOpType = 'dynamically weighted'
      if (IFXMS) then
        calctype = 'XDW-CASPT2'
      else
        calctype = 'DW-CASPT2'
      end if
      if (HZERO == 'DYALL') calctype = 'QD-NEVPT2'
    else if (IFRMS) then
      FockOpType = 'state-specific'
      calctype = 'RMS-CASPT2'
    else if (IFXMS) then
      FockOpType = 'state-average'
      calctype = 'XMS-CASPT2'
    else
      FockOpType = 'state-specific'
      if (IFSADREF) FockOpType = 'state-average'
      if (HZERO /= 'DYALL') then
        calctype = 'MS-CASPT2'
      else
        calctype = 'QD-NEVPT2'
      end if
    end if
  else
    FockOpType = 'state-specific'
    if (IFSADREF) FockOpType = 'state-average'
    if (HZERO /= 'DYALL') then
      calctype = 'SS-CASPT2'
    else
      calctype = 'SS-NEVPT2'
    end if
  end if

  write(u6,fmt2//'A,T50,A)') 'Type of calculation',trim(calctype)

  if (Hzero == 'DYALL') then
    if (do_grad) then
      if (SC_prop) then
        write(u6,fmt2//'A,T50,A)') 'Internal contraction for properties','strongly contracted (SC-NEVPT2)'
      else
        write(u6,fmt2//'A,T50,A)') 'Internal contraction for properties','partially contracted (PC-NEVPT2)'
      end if
    else
      if (Do_FIC .and. Do_SC) then
        write(u6,fmt2//'A,T50,A)') 'Internal contraction for energies','partially and strongly contracted'
        write(u6,fmt2//'A,T50,A)') '','(PC-NEVPT2 and SC-NEVPT2)'
      else if (Do_FIC) then
        write(u6,fmt2//'A,T50,A)') 'Internal contraction for energies','partially contracted'
        write(u6,fmt2//'A,T50,A)') '','(PC-NEVPT2)'
      else if (Do_SC) then
        write(u6,fmt2//'A,T50,A)') 'Internal contraction for energies','partially contracted'
        write(u6,fmt2//'A,T50,A)') '','(PC-NEVPT2)'
      end if
    end if
    if (Do_SC) write(u6,fmt2//'A,T50,ES8.2)') 'Denominator threshold for SC-NEVPT2',SC_thres
  end if

  write(u6,fmt2//'A,T50,A)') 'Fock operator',trim(FockOpType)
  if (IFDW) then
    write(u6,fmt2//'A,T45,I6)') 'DW Type',DWType
    if (zeta >= 0) then
      write(u6,fmt2//'A,T50,ES10.4)') 'DW exponent',zeta
    else
      write(u6,fmt2//'A,T50,A)') 'DW exponent','infinity'
    end if
  end if

  if (Hzero /= 'STANDARD') write(u6,fmt2//'A,T50,A)') '0th-order Hamiltonian',trim(Hzero)

  write(u6,fmt2//'A,T45,F9.2)') 'IPEA shift',ipea_shift
  write(u6,fmt2//'A,T45,F9.2)') 'Real shift',real_shift
  write(u6,fmt2//'A,T45,F9.2)') 'Imaginary shift',imag_shift
  if (sigma_p_epsilon > Zero) then
    if (sigma_p_exponent == 1) then
      write(u6,fmt2//'A,T45,F9.2)') 'Sigma^1 regularizer',sigma_p_epsilon
    else
      write(u6,fmt2//'A,T45,F9.2)') 'Sigma^2 regularizer',sigma_p_epsilon
    end if
  end if

  if (ORBIN == 'TRANSFOR') then
    write(u6,fmt1) 'The input orbitals will be transformed to quasi-canonical'
  else
    write(u6,fmt1) 'The input orbitals will not be transformed to quasi-canonical'
  end if

  if ((IFXMS .or. IFRMS) .and. (HZERO /= 'DYALL')) &
    write(u6,fmt1) 'The input states will be rotated to diagonalize the Fock operator'

  if (IFDORTHO) write(u6,fmt1) 'Canonical orthornormalization will be used for the IC basis'

  if (do_grad) then
    if (do_nac) then
      if (do_csf) then
        write(u6,fmt1) 'Quantities for analytical NAC with CSF term will be calculated'
      else
        write(u6,fmt1) 'Quantities for analytical NAC without CSF term will be calculated'
      end if
    else
      write(u6,fmt1) 'Quantities for analytical gradients will be calculated'
    end if
  end if

  call CollapseOutput(0,Line)
  write(u6,*)
end if

end subroutine prinp_caspt2
