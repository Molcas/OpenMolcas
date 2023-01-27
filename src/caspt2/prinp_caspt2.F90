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
  use definitions, only: iwp, wp
  use caspt2_output, only: iPrGlb, terse, usual, verbose
  use caspt2_global, only: sigma_p_epsilon, sigma_p_exponent, &
                           ipea_shift, imag_shift, real_shift
  use caspt2_gradient, only: do_grad, do_nac, do_csf

  implicit none

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

  integer(kind=iwp)  :: i,iSym,left,lLine,lPaper
  character(len=8)   :: fmt1, fmt2
  character(len=120) :: Line
  character(len=3)   :: lIrrep(8)
  character(len=20)  :: calctype, FockOpType
!----------------------------------------------------------------------*
!     Start and define the paper width,                                *
!     initialize blank and header lines                                *
!----------------------------------------------------------------------*
  Line = ' '
  lLine = Len(Line)
  lPaper = 132
  left = (lPaper-lLine)/2
  write(fmt1,'(A,I3.3,A)') '(',left,'X,A)'
  write(fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
!     Print the ONEINT file identifier                                 *
!----------------------------------------------------------------------*
  if (iprglb >= verbose) then
    write(6,*)
    write(6,fmt1) 'Header of the ONEINT file:'
    write(6,fmt1) '--------------------------'
    write(Line,'(36A2)')(Header(i),i=1,36)
    write(6,fmt1) trim(adjustl(Line))
    write(Line,'(36A2)')(Header(i),i=37,72)
    write(6,fmt1) trim(adjustl(Line))
    write(6,*)
  end if
!----------------------------------------------------------------------*
!     Print cartesian coordinates of the system                        *
!----------------------------------------------------------------------*
  if (iprglb >= verbose) then
    call prCoor()
  end if
!----------------------------------------------------------------------*
!     Print orbital and wavefunction specifications                    *
!----------------------------------------------------------------------*
  if (iprglb >= usual) then
    write(6,*)
    Line = ' '
    write(Line(left-2:),'(A)') 'Wave function specifications:'
    call CollapseOutput(1,Line)
    write(6,fmt1) '-----------------------------'
    write(6,*)
    write(6,fmt2//'A,T45,I6)') 'Number of closed shell electrons', 2*NISHT
    write(6,fmt2//'A,T45,I6)') 'Number of electrons in active shells', NACTEL
    write(6,fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space', NHOLE1
    write(6,fmt2//'A,T45,I6)') 'Max number of electrons in RAS3 space', NELE3
    write(6,fmt2//'A,T45,I6)') 'Number of inactive orbitals', NISHT
    write(6,fmt2//'A,T45,I6)') 'Number of active orbitals', NASHT
    write(6,fmt2//'A,T45,I6)') 'Number of secondary orbitals', NSSHT
    write(6,fmt2//'A,T45,F6.1)') 'Spin quantum number', 0.5_wp * real(ISPIN-1, kind=wp)
    write(6,fmt2//'A,T45,I6)') 'State symmetry', STSYM
    write(6,fmt2//'A,T40,I11)') 'Number of CSFs', NCONF
    write(6,fmt2//'A,T45,I6)') 'Number of CASSCF root(s) available', NROOTS
    write(6,fmt2//'A,T45,I6)') 'CASPT2 state passed to geometry opt.', iRlxRoot
    if (ifmix) then
      write(6,fmt2//'A,T45,10I3)') 'A file JOBMIX will be created'
    end if
    if (nstate > 1) then
      write(6,fmt1) 'This is a MULTI-STATE CASSCF reference'
      write(6,fmt2//'A,T45,I6)') 'Number of CI roots used', NSTATE
      write(6,fmt2//'A,(T47,10I4))') 'These are:', (MSTATE(I),I=1,NSTATE)
      if (ifmscoup) then
        write(6,fmt1) 'Off-diagonal elements of Heff are computed'
      else
        write(6,fmt1) 'Heff is assumed diagonal'
      end if
    else
      if (iscf == 0) then
        write(6,fmt1) 'This is a CASSCF or RASSCF reference function'
#ifdef _ENABLE_BLOCK_DMRG_
        if (DoCumulant) then
          write(6,fmt1) 'Using 4-RDM cumulant approximation,'// &
                        ' activated by 3RDM keyword in RASSCF'
        end if
#elif _ENABLE_CHEMPS2_DMRG_
        if (DoCumulant) then
          write(6,fmt1) 'This is a DMRG reference with exact 4-RDM,'// &
                        ' activated by 3RDM keyword in RASSCF'
        end if
#endif
      else if (iscf == 1) then
        write(6,fmt1) 'This is a closed shell RHF reference function'
      else
        write(6,fmt1) 'This is a high spin open shell RHF reference function'
      end if
    end if
    call CollapseOutput(0,'Wave function specifications:')
  end if

  call Get_cArray('Irreps',lIrrep,24)
  do iSym = 1,nSym
    lIrrep(iSym) = adjustr(lIrrep(iSym))
  end do

  if (iprglb >= usual) then
    write(6,*)
    Line = ' '
    write(Line(left-2:),'(A)') 'Orbital specifications:'
    call CollapseOutput(1,Line)
    write(6,fmt1) '-----------------------'
    write(6,*)
    write(6,fmt2//'A,T47,8I4)') 'Symmetry species', (iSym,iSym=1,nSym)
    write(6,fmt2//'A,T47,8(1X,A))') '                ', (lIrrep(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Frozen orbitals', (nFro(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Inactive orbitals', (nIsh(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Active orbitals', (nAsh(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Secondary orbitals', (nSsh(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Deleted orbitals', (nDel(iSym),iSym=1,nSym)
    write(6,fmt2//'A,T47,8I4)') 'Number of basis functions', (nBas(iSym),iSym=1,nSym)
    call CollapseOutput(0,'Orbital specifications:')
  end if
!----------------------------------------------------------------------*
!     Print routing information                                        *
!----------------------------------------------------------------------*
  if (iprglb >= terse) then
    if (rfpert) then
      write(6,*)
      write(6,fmt1) 'Reaction field specifications:'
      write(6,fmt1) '------------------------------'
      write(6,*)
      write(6,'(6X,A)') 'An external reaction field was determined'// &
            ' previously and added to the one-electron Hamiltonian'
      write(6,'(6X,A)') 'It will not be reevaluated even though'//    &
                  ' dynamic correlation may change the density.'
      write(6,*)
    end if

    write(6,*)
    Line = ' '
    write(Line(left-2:),'(A)') 'CASPT2 specifications:'
    call CollapseOutput(1,Line)
    write(6,fmt1) '----------------------'
    write(6,*)

    if (IFMSCOUP) then
      if (IFDW) then
        FockOpType = 'dynamically weighted'
        if (IFXMS) then
          calctype = 'XDW-CASPT2'
        else
          calctype = 'DW-CASPT2'
        end if
      else if (IFRMS) then
        FockOpType = 'state-specific'
        calctype = 'RMS-CASPT2'
      else
        if (IFXMS) then
          FockOpType = 'state-average'
          calctype = 'XMS-CASPT2'
        else
          FockOpType = 'state-specific'
          if (IFSADREF) FockOpType = 'state-average'
          calctype = 'MS-CASPT2'
        end if
      end if
    else
      FockOpType = 'state-specific'
      if (IFSADREF) FockOpType = 'state-average'
      calctype = 'SS-CASPT2'
    end if

    write(6,fmt2//'A,T50,A)') 'Type of calculation', trim(calctype)

    write(6,fmt2//'A,T50,A)') 'Fock operator', trim(FockOpType)
    if (IFDW) then
      write(6,fmt2//'A,T45,I6)') 'DW Type', DWType
      if (zeta >= 0) then
        write(6,fmt2//'A,T50,E10.4)') 'DW exponent', zeta
      else
        write(6,fmt2//'A,T50,A)') 'DW exponent','infinity'
      end if
    end if

    if (Hzero /= 'STANDARD') then
      write(6,fmt2//'A,T50,A)') '0th-order Hamiltonian', trim(Hzero)
    end if

    write(6,fmt2//'A,T45,F9.2)') 'IPEA shift', ipea_shift
    write(6,fmt2//'A,T45,F9.2)') 'Real shift', real_shift
    write(6,fmt2//'A,T45,F9.2)') 'Imaginary shift', imag_shift
    if (sigma_p_epsilon > 0.0_wp) then
      if (sigma_p_exponent == 1) then
        write(6,fmt2//'A,T45,F9.2)') 'Sigma^1 regularizer', sigma_p_epsilon
      else
        write(6,fmt2//'A,T45,F9.2)') 'Sigma^2 regularizer', sigma_p_epsilon
      end if
    end if

    if (ORBIN == 'TRANSFOR') then
      write(6,fmt1) 'The input orbitals will be transformed to quasi-canonical'
    else
      write(6,fmt1) 'The input orbitals will not be transformed to quasi-canonical'
    end if

    if (IFXMS .or. IFRMS) then
      write(6,fmt1) 'The input states will be rotated to diagonalize the Fock operator'
    end if

    if (IFDORTHO) then
      write(6,fmt1) 'Unscaled orthornormalization will be used for the IC basis'
    end if

    if (do_grad .and. (.not. do_nac)) then
      write(6,fmt1) 'Quantities for analytical gradients will be calculated'
    end if

    if (do_nac) then
      if (do_csf) then
        write(6,fmt1) 'Quantities for analytical NAC with CSF term will be calculated'
      else
        write(6,fmt1) 'Quantities for analytical NAC without CSF term will be calculated'
      end if
    end if

    call CollapseOutput(0,'CASPT2 specifications:')
    write(6,*)
  end if

end subroutine prinp_caspt2
