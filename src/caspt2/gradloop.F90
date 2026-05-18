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
! Copyright (C) 2023, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine GradLoop(Heff,Ueff,H0,U0,H0Sav,nState)
! Gradient loop for MS-CASPT2 variants
! Usually, we do not solve the CASPT2 equation again; the excitation
! amplitude etc. have been stored on disk (should be avoided, though)
! in the first loop (for energy) and are restored in the second loop
! below (by calling SavGradParams)

use caspt2_global, only: iPrGlb
use caspt2_global, only: do_grad, IDSAVGRD, iStpGrd
use PrintLevel, only: USUAL, VERBOSE
use EQSOLV, only: iRHS, iVecC, iVecC2, iVecR, iVecW, iVecX
use caspt2_module, only: CPUEIG, CPUFG3, CPUFMB, CPUGIN, CPUGRD, CPUINT, CPULCS, CPUNAD, CPUOVL, CPUPCG, CPUSCA, CPUSER, CPUSGM, &
                         CPUSIN, CPUVEC, CPUPRP, CPUPT2, CPURHS, CPUSBM, TIOEIG, TIOFG3, TIOFMB, TIOGIN, TIOGRD, TIOINT, TIOLCS, &
                         TIONAD, TIOOVL, TIOPCG, TIOSCA, TIOSER, TIOSGM, TIOSIN, TIOVEC, TIOPRP, TIOPT2, TIORHS, TIOSBM, Energy, &
                         IfChol, IfDens, IfDW, IfMSCoup, IfProp, IfRMS, IfXMS, iRlxRoot, jState, nGroup, nGroupState, mState
use constants, only: Zero
use definitions, only: iwp, wp, u6

implicit none
#include "warnings.h"
integer(kind=iwp), intent(in) :: nState
! Effective Hamiltonian
!real(kind=wp), allocatable :: Heff(:,:), Ueff(:,:)
real(kind=wp), intent(inout) :: Heff(nState,nState), Ueff(nState,nState)
! Zeroth-order Hamiltonian
!real(kind=wp), allocatable :: H0(:,:), U0(:,:)
real(kind=wp), intent(inout) :: H0(nState,nState), U0(nState,nState), H0Sav(nState,nState)
character(len=60) :: STLNE2
! Timers
real(kind=wp) :: CPTF0, CPTF11, CPTF12, CPTF13, CPTF14, TIOTF0, TIOTF11, TIOTF12, TIOTF13, TIOTF14, CPE, CPUTOT, TIOE, TIOTOT
! Indices
integer(kind=iwp) :: ISTATE, IGROUP, JSTATE_OFF
! Convergence check
integer(kind=iwp) :: ICONV
! For verification only
integer(kind=iwp) LAXITY
integer(kind=iwp), external :: Cho_X_GetTol

IDSAVGRD = 0

if (iStpGrd == 0) then
  ! just for verification
  LAXITY = 8
  if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
  call Add_Info('E_MSPT2',ENERGY,nState,LAXITY)
end if

iStpGrd = 2

! For (X)Multi-State, a long loop over root states.
! The states are ordered by group, with each group containing a number
! of group states for which GRPINI is called.
JSTATE_OFF = 0
stateloop2: do IGROUP=1,NGROUP

  if (IPRGLB >= USUAL) then
    write(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
    call CollapseOutput(1,trim(STLNE2))
    write(u6,*)
  end if

  call GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0,nState)

  if (do_grad) call CNSTFIFAFIMO(1)

  do ISTATE=1,NGROUPSTATE(IGROUP)
    JSTATE = JSTATE_OFF+ISTATE

    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    !call STINI(JSTATE)
    call RHS_INIT() !! somehow
    call SavGradParams(2,IDSAVGRD)
    call SavGradParams2(2,UEFF,U0,H0,nState)
    if ((IFXMS .and. IFDW) .or. IFRMS) H0(:,:) = H0Sav(:,:)
    call TIMING(CPTF11,CPE,TIOTF11,TIOE)
    CPUSIN = CPTF11-CPTF0
    TIOSIN = TIOTF11-TIOTF0

    ! Solve CASPT2 equation system and compute corr energies.
    if (IPRGLB >= USUAL) then
      write(u6,'(A)') repeat('*',80)
      write(u6,*) ' CASPT2 EQUATION SOLUTION (SECOND RUN)'
      write(u6,'(A)') repeat('-',80)
    end if

    write(STLNE2,'(A,I0)') 'Solve CASPT2 eqs. for state ',MSTATE(JSTATE)
    call StatusLine('CASPT2: ',STLNE2)
    call EQCTL2(ICONV)

    ! Save the final caspt2 energy in the global array ENERGY():
    !ENERGY(JSTATE) = E2TOT

    call TIMING(CPTF12,CPE,TIOTF12,TIOE)
    CPUPT2 = CPTF12-CPTF11
    TIOPT2 = TIOTF12-TIOTF11

    ! Orbitals, properties:
    ! if the dens keyword is used, need accurate density and
    ! for that the serial LUSOLV file is needed, in that case copy
    ! the distributed LURHS() to LUSOLV here.
    if (IFDENS .or. (do_grad .and. ((iRlxRoot == MSTATE(JSTATE)) .or. IFMSCOUP))) then
      call PCOLLVEC(IRHS,0)
      call PCOLLVEC(IVECX,0)
      call PCOLLVEC(IVECR,0)
      call PCOLLVEC(IVECC,1)
      call PCOLLVEC(IVECC2,1)
      call PCOLLVEC(IVECW,1)
    end if

    if (IFPROP .or. (do_grad .and. ((IRLXroot == MSTATE(JSTATE)) .or. IFMSCOUP))) then

      call PRPCTL(0,UEFF,U0,nState)

    else
      if (IPRGLB >= USUAL) then
        write(u6,*)
        write(u6,*) '  (Skipping property calculation,'
        write(u6,*) '   use PROP keyword to activate)'
      end if
    end if

    call TIMING(CPTF13,CPE,TIOTF13,TIOE)
    CPUPRP = CPTF13-CPTF12
    TIOPRP = TIOTF13-TIOTF12

    call TIMING(CPTF14,CPE,TIOTF14,TIOE)
    CPUGRD = CPTF14-CPTF13
    TIOGRD = TIOTF14-TIOTF13

    CPUTOT = CPTF14-CPTF0
    TIOTOT = TIOTF14-TIOTF0

    call Iter_Timing()

    ! End of long loop over states in the group
  end do

  if (IPRGLB >= USUAL) then
    call CollapseOutput(0,'CASPT2 computation for group ')
    write(u6,*)
  end if
  ! End of long loop over groups
  JSTATE_OFF = JSTATE_OFF+NGROUPSTATE(IGROUP)
end do stateloop2

contains

subroutine Iter_Timing()
  if (ISTATE == 1) then
    CPUTOT = CPUTOT+CPUGIN
    TIOTOT = TIOTOT+TIOGIN
  else
    CPUGIN = Zero
    TIOGIN = Zero
    CPUFMB = Zero
    TIOFMB = Zero
    CPUINT = Zero
    TIOINT = Zero
  end if

  if (IPRGLB >= VERBOSE) then
    write(u6,*)
    write(u6,'(A,I6)') '  CASPT2 TIMING INFO FOR STATE ',MSTATE(JSTATE)
    write(u6,*)
    write(u6,'(A)') '                         cpu time  (s)  wall time (s) '
    write(u6,'(A)') '                         -------------  ------------- '
    write(u6,*)
    write(u6,'(A,2F14.2)') '  Group initialization  ',CPUGIN,TIOGIN
    write(u6,'(A,2F14.2)') '  - Fock matrix build   ',CPUFMB,TIOFMB
    write(u6,'(A,2F14.2)') '  - integral transforms ',CPUINT,TIOINT
    write(u6,'(A,2F14.2)') '  State initialization  ',CPUSIN,TIOSIN
    write(u6,'(A,2F14.2)') '  - density matrices    ',CPUFG3,TIOFG3
    write(u6,'(A,2F14.2)') '  CASPT2 equations      ',CPUPT2,TIOPT2
    write(u6,'(A,2F14.2)') '  - H0 S/B matrices     ',CPUSBM,TIOSBM
    write(u6,'(A,2F14.2)') '  - H0 S/B diag         ',CPUEIG,TIOEIG
    write(u6,'(A,2F14.2)') '  - H0 NA diag          ',CPUNAD,TIONAD
    write(u6,'(A,2F14.2)') '  - RHS construction    ',CPURHS,TIORHS
    write(u6,'(A,2F14.2)') '  - PCG solver          ',CPUPCG,TIOPCG
    write(u6,'(A,2F14.2)') '    - scaling           ',CPUSCA,TIOSCA
    write(u6,'(A,2F14.2)') '    - lin. comb.        ',CPULCS,TIOLCS
    write(u6,'(A,2F14.2)') '    - inner products    ',CPUOVL,TIOOVL
    write(u6,'(A,2F14.2)') '    - basis transforms  ',CPUVEC,TIOVEC
    write(u6,'(A,2F14.2)') '    - sigma routines    ',CPUSGM,TIOSGM
    write(u6,'(A,2F14.2)') '  - array collection    ',CPUSER,TIOSER
    write(u6,'(A,2F14.2)') '  Properties            ',CPUPRP,TIOPRP
    write(u6,'(A,2F14.2)') '  MS coupling           ',CPUGRD,TIOGRD
    write(u6,'(A,2F14.2)') ' Total time             ',CPUTOT,TIOTOT
    write(u6,*)
  end if
end subroutine Iter_Timing

end subroutine GradLoop
