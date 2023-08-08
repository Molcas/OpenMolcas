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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

Subroutine GradLoop(IRETURN,Heff,Ueff,H0,U0,H0Sav)

  USE SUPERINDEX
  USE PT2WFN
  use caspt2_output, only: iPrGlb, terse, usual, verbose
  use caspt2_gradient, only: do_grad, IDSAVGRD, iStpGrd
  use definitions, only: iwp,wp

  Implicit None

#include "rasdim.fh"
#include "warnings.h"
#include "constants.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "stdalloc.fh"
#include "caspt2_grad.fh"

  Integer(kind=iwp), intent(out) :: IRETURN
  character(len=60) :: STLNE2

! Timers
  Real(kind=wp) :: CPTF0, CPTF10, CPTF11, CPTF12, CPTF13, CPTF14, &
     &            TIOTF0,TIOTF10,TIOTF11,TIOTF12,TIOTF13,TIOTF14, &
     &               CPE,CPUTOT,TIOE,TIOTOT
! Indices
  Integer(kind=iwp) :: I,ISTATE,IGROUP,JSTATE_OFF
! Convergence check
  Integer(kind=iwp) :: ICONV
! Effective Hamiltonian
! Real(kind=wp), Allocatable :: Heff(:,:), Ueff(:,:)
  Real(kind=wp) :: Heff(*), Ueff(*)
! Zeroth-order Hamiltonian
! Real(kind=wp), Allocatable :: H0(:,:), U0(:,:)
  Real(kind=wp) :: H0(*), U0(*), H0Sav(*)

  IDSAVGRD = 0
  iStpGrd = 2_iwp

! For (X)Multi-State, a long loop over root states.
! The states are ordered by group, with each group containing a number
! of group states for which GRPINI is called.
    JSTATE_OFF=0
    STATELOOP2: DO IGROUP=1,NGROUP

    IF (IPRGLB.GE.USUAL) THEN
      WRITE(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
      CALL CollapseOutput(1,TRIM(STLNE2))
      WRITE(6,*)
    END IF

    CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
    CALL GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0)
    If ((IFXMS.and.IFDW).OR.IFRMS) Call DCopy_(nState*nState,H0Sav,1,H0,1)
    CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
    CPUGIN=CPTF10-CPTF0
    TIOGIN=TIOTF10-TIOTF0

    If (do_grad) CALL CNSTFIFAFIMO(1)

    DO ISTATE=1,NGROUPSTATE(IGROUP)
      JSTATE = JSTATE_OFF + ISTATE

      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      !CALL STINI
      CALL RHS_INIT !! somehow
      Call SavGradParams(2,IDSAVGRD)
      Call SavGradParams2(2,UEFF,U0,H0)
      CALL TIMING(CPTF11,CPE,TIOTF11,TIOE)
      CPUSIN=CPTF11-CPTF0
      TIOSIN=TIOTF11-TIOTF0

! Solve CASPT2 equation system and compute corr energies.
      IF (IPRGLB.GE.USUAL) THEN
         WRITE(6,'(20A4)')('****',I=1,20)
         WRITE(6,*)' CASPT2 EQUATION SOLUTION (SECOND RUN)'
         WRITE(6,'(20A4)')('----',I=1,20)
      END IF

      Write(STLNE2,'(A27,I3)')'Solve CASPT2 eqs for state ', MSTATE(JSTATE)
      Call StatusLine('CASPT2:',TRIM(STLNE2))
      CALL EQCTL2(ICONV)

! Save the final caspt2 energy in the global array ENERGY():
!     ENERGY(JSTATE)=E2TOT

      CALL TIMING(CPTF12,CPE,TIOTF12,TIOE)
      CPUPT2=CPTF12-CPTF11
      TIOPT2=TIOTF12-TIOTF11

      IF (ICONV .NE. 0) THEN
! No convergence. Skip the rest of the calculation.
        IRETURN = _RC_NOT_CONVERGED_
        EXIT STATELOOP2
      END IF

! Orbitals, properties:
      ! if the dens keyword is used, need accurate density and
      ! for that the serial LUSOLV file is needed, in that case copy
      ! the distributed LURHS() to LUSOLV here.
      IF (IFDENS.OR.(do_grad.and.(iRlxRoot.eq.MSTATE(JSTATE).or.IFMSCOUP))) THEN
        CALL PCOLLVEC(IRHS,0)
        CALL PCOLLVEC(IVECX,0)
        CALL PCOLLVEC(IVECR,0)
        CALL PCOLLVEC(IVECC,1)
        CALL PCOLLVEC(IVECC2,1)
        CALL PCOLLVEC(IVECW,1)
      END IF

      IF (IFPROP.OR.(do_grad.and.(IRLXroot.eq.MSTATE(JSTATE).or.IFMSCOUP))) THEN
        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,'(20A4)')('****',I=1,20)
          WRITE(6,*)' CASPT2 PROPERTY SECTION'
        END IF
        CALL PRPCTL(UEFF,U0)
      ELSE
        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,*)'  (Skipping property calculation,'
          WRITE(6,*)'   use PROP keyword to activate)'
        END IF
      END IF

      CALL TIMING(CPTF13,CPE,TIOTF13,TIOE)
      CPUPRP=CPTF13-CPTF12
      TIOPRP=TIOTF13-TIOTF12

      CALL TIMING(CPTF14,CPE,TIOTF14,TIOE)
      CPUGRD=CPTF14-CPTF13
      TIOGRD=TIOTF14-TIOTF13

      CPUTOT=CPTF14-CPTF0
      TIOTOT=TIOTF14-TIOTF0

      IF (ISTATE.EQ.1) THEN
        CPUTOT=CPUTOT+CPUGIN
        TIOTOT=TIOTOT+TIOGIN
      ELSE
        CPUGIN=0.0D0
        TIOGIN=0.0D0
        CPUFMB=0.0D0
        TIOFMB=0.0D0
        CPUINT=0.0D0
        TIOINT=0.0D0
      END IF

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(A,I6)')    '  CASPT2 TIMING INFO FOR STATE ',MSTATE(JSTATE)
        WRITE(6,*)
        WRITE(6,'(A)')       '                        '// &
     &                       ' cpu time  (s) '// &
     &                       ' wall time (s) '
        WRITE(6,'(A)')       '                        '// &
     &                       ' ------------- '// &
     &                       ' ------------- '
        WRITE(6,*)
        WRITE(6,'(A,2F14.2)')'  Group initialization  ',CPUGIN,TIOGIN
        WRITE(6,'(A,2F14.2)')'  - Fock matrix build   ',CPUFMB,TIOFMB
        WRITE(6,'(A,2F14.2)')'  - integral transforms ',CPUINT,TIOINT
        WRITE(6,'(A,2F14.2)')'  State initialization  ',CPUSIN,TIOSIN
        WRITE(6,'(A,2F14.2)')'  - density matrices    ',CPUFG3,TIOFG3
        WRITE(6,'(A,2F14.2)')'  CASPT2 equations      ',CPUPT2,TIOPT2
        WRITE(6,'(A,2F14.2)')'  - H0 S/B matrices     ',CPUSBM,TIOSBM
        WRITE(6,'(A,2F14.2)')'  - H0 S/B diag         ',CPUEIG,TIOEIG
        WRITE(6,'(A,2F14.2)')'  - H0 NA diag          ',CPUNAD,TIONAD
        WRITE(6,'(A,2F14.2)')'  - RHS construction    ',CPURHS,TIORHS
        WRITE(6,'(A,2F14.2)')'  - PCG solver          ',CPUPCG,TIOPCG
        WRITE(6,'(A,2F14.2)')'    - scaling           ',CPUSCA,TIOSCA
        WRITE(6,'(A,2F14.2)')'    - lin. comb.        ',CPULCS,TIOLCS
        WRITE(6,'(A,2F14.2)')'    - inner products    ',CPUOVL,TIOOVL
        WRITE(6,'(A,2F14.2)')'    - basis transforms  ',CPUVEC,TIOVEC
        WRITE(6,'(A,2F14.2)')'    - sigma routines    ',CPUSGM,TIOSGM
        WRITE(6,'(A,2F14.2)')'  - array collection    ',CPUSER,TIOSER
        WRITE(6,'(A,2F14.2)')'  Properties            ',CPUPRP,TIOPRP
        WRITE(6,'(A,2F14.2)')'  MS coupling           ',CPUGRD,TIOGRD
        WRITE(6,'(A,2F14.2)')' Total time             ',CPUTOT,TIOTOT
        WRITE(6,*)
      END IF

! End of long loop over states in the group
    END DO

    IF (IPRGLB.GE.USUAL) THEN
      CALL CollapseOutput(0,'CASPT2 computation for group ')
      WRITE(6,*)
    END IF
! End of long loop over groups
    JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
  END DO STATELOOP2

End Subroutine GradLoop
