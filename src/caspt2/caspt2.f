************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE CASPT2(IRETURN)
      USE SUPERINDEX
      USE INPUTDATA
      USE PT2WFN
      IMPLICIT NONE
      INTEGER IRETURN
*----------------------------------------------------------------------*
*     1998  PER-AAKE MALMQUIST                                         *
*     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
*     UNIVERSITY OF LUND, SWEDEN                                       *
*----------------------------------------------------------------------*

C     SECOND ORDER PERTURBATION CALCULATIONS WITH A CASSCF
C     REFERENCE FUNCTION.

C     ORIGINAL CASPT2 PROGRAM WRITTEN 890526 BY:
C     KERSTIN ANDERSSON (ALL CASPT2 CODE WITH FOLLOWING EXCEPTIONS)
C     PER-AKE MALMQVIST (THE GUGA PART)
C     BJORN O ROOS      (THE INTEGRAL TRANSFORMATION)

C     MODIFIED 1991-02-23 BY PER-AAKE MALMQUIST, FOR USE WITH THE
C     MOLCAS RASSCF PROGRAM.

C     ALL PROGRAM REWRITTEN (MALMQVIST 1993) FOR MOLCAS VERSION 3,
C     EXCEPT TRACTL AND TRA2 KEPT BUT SLIGHTLY MODIFIED, AND ALSO:
C     INTERFACING WITH MOLCAS-3 BY MARCUS FUELSCHER

C     REWRITTEN FOR 1) EXACT ACTIVE DENSITY MATRIX
C     2) QUANTITIES NEEDED FOR GRADIENTS
C     3) MULTI-STATE CASPT2
C     BY MALMQVIST 1998.
C
C     SINCE THEN, THE FOLLOWING MODIFICATIONS HAVE BEEN MADE:
C     IPEA HAMILTONIAN BY G. GHIGO & P. MALMQVIST
C       Chem. Phys. Lett. 396, pp 142 (2004)
C       http://dx.doi.org/10.1016/j.cplett.2004.08.032
C     LOVCASPT2/FNOCASPT2/AFREEZE BY B. ROOS & F. AQUILANTE
C       J. Chem. Phys. 131, 034113 (2009)
C       http://dx.doi.org/10.1063/1.3157463
C     CHOLESKY SUPPORT BY P. MALMQVIST AND F. AQUILANTE
C       J. Chem. Theory Comput. 4 (5), pp 694-702 (2008)
C       http://dx/doi.org/10.1021/ct700263h
C     RASPT2 BY P. MALMQVIST
C       J. Chem. Phys. 128, 204109 (2008)
C       http://dx.doi.org/10.1063/1.2920188
C     PARALLELIZATION BY S. VANCOILLIE
C       J. Comp. Chem. 34, pp 1937-1948 (2013)
C       http://dx.doi.org/10.1002/jcc.23342
C     MAJOR REFACTORING OF INITIALIZATION PHASE AND ADDITION OF
C     HDF5 SUPPORT BY S. VANCOILLIE (2016)
C
************************************************************************
#include "rasdim.fh"
#include "warnings.fh"
#include "constants.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "stdalloc.fh"
      CHARACTER(60) STLNE2
#ifdef _MOLCAS_MPP_
      LOGICAL KING, Is_Real_Par
#endif
* Timers
      REAL*8  CPTF0, CPTF10, CPTF11, CPTF12, CPTF13, CPTF14,
     &       TIOTF0,TIOTF10,TIOTF11,TIOTF12,TIOTF13,TIOTF14,
     &          CPE,CPUTOT,TIOE,TIOTOT
* Indices
      INTEGER I
      INTEGER ISTATE
      INTEGER IGROUP,JSTATE_OFF
* Convergence check
      INTEGER ICONV
* Relative energies
      REAL*8  RELAU,RELEV,RELCM,RELKJ

* Effective Hamiltonian
      REAL*8, ALLOCATABLE :: Heff(:,:), Ueff(:,:)

* Zeroth-order Hamiltonian
      REAL*8, ALLOCATABLE :: H0(:,:), U0(:,:)

      Call StatusLine('CASPT2:','Just starting')

      IRETURN = 0
      CALL QENTER('CASPT2')

      CALL SETTIM
      ! CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)

* Probe the environment to globally set the IPRGLB value
      Call Set_Print_Level

*======================================================================*
*
      Call StatusLine('CASPT2:','Initializing')
      CALL PT2INI
* Initialize effective Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(Heff,Nstate,Nstate)
      CALL MMA_ALLOCATE(Ueff,Nstate,Nstate)
      Heff=0.0D0
      Ueff=0.0D0
* Initialize zeroth-order Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(H0,Nstate,Nstate)
      CALL MMA_ALLOCATE(U0,Nstate,Nstate)
      H0=0.0D0
* U0 is initialized as the identity matrix, in the case of a
* standard MS-CASPT2 calculation it will not be touched anymore
      U0=0.0D0
      call dcopy_(Nstate,[1.0d0],0,U0,Nstate+1)
*
*======================================================================*
* If the EFFE keyword has been used, we already have the multi state
* coupling Hamiltonian effective matrix, just copy the energies and
* proceed to the MS coupling section.
* Otherwise, put the CASSCF energies on the diagonal, i.e. form the
* first-order corrected effective Hamiltonian:
*     Heff[1] = PHP
* and later on we will add the second-order correction
* Heff(2) = PH \Omega_1 P to Heff[1]
      IF (INPUT%JMS) THEN
        DO I=1,NSTATE
          ENERGY(I)=INPUT%HEFF(I,I)
        END DO
        HEFF(:,:)=INPUT%HEFF(:,:)
        GOTO 1000
      ELSE
        DO I=1,NSTATE
          HEFF(I,I) = REFENE(I)
        END DO
        IF (IPRGLB.GE.DEBUG) THEN
          write(6,*)' Heff[1] in the original model space basis:'
          call prettyprint(Heff,Nstate,Nstate)
        END IF
      END IF

* In case of a XDW-CASPT2 calculation we first rotate the CASSCF
* states according to the XMS prescription in xdwinit
      if (IFXMS.and.IFDW) then
        call xdwinit(Heff,H0,U0)
        if (IFEFOCK) then
          call wgtinit(H0)
        else
          call wgtinit(Heff)
        end if
      else
        call wgtinit(Heff)
      end if

* Before entering the long loop over groups and states, precompute
* the 1-RDMs for all states and mix them according to the type of
* calculation: MS, XMS, DW or XDW.
      call rdminit

* For (X)Multi-State, a long loop over root states.
* The states are ordered by group, with each group containing a number
* of group states for which GRPINI is called.
      JSTATE_OFF=0
      STATELOOP: DO IGROUP=1,NGROUP

       IF ((NLYGROUP.NE.0).AND.(IGROUP.NE.NLYGROUP)) THEN
         JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
         CYCLE
       END IF

       IF (IPRGLB.GE.USUAL) THEN
        If(.not.IFNOPT2) Then
         WRITE(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
         CALL CollapseOutput(1,TRIM(STLNE2))
         WRITE(6,*)
        End If
       END IF

       CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
       CALL GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0)
       CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
       CPUGIN=CPTF10-CPTF0
       TIOGIN=TIOTF10-TIOTF0
CProducing XMS Rotated States
       If(IFNOPT2) then
        GOTO 9999
       END IF

       DO ISTATE=1,NGROUPSTATE(IGROUP)
         JSTATE = JSTATE_OFF + ISTATE


* Skip this state if we only need 1 state and it isn't this "one".
         IF ((NLYROOT.NE.0).AND.(JSTATE.NE.NLYROOT)) CYCLE

         CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
         CALL STINI
         CALL TIMING(CPTF11,CPE,TIOTF11,TIOE)
         CPUSIN=CPTF11-CPTF0
         TIOSIN=TIOTF11-TIOTF0

* Solve CASPT2 equation system and compute corr energies.
         IF (IPRGLB.GE.USUAL) THEN
            WRITE(6,'(20A4)')('****',I=1,20)
            WRITE(6,*)' CASPT2 EQUATION SOLUTION'
            WRITE(6,'(20A4)')('----',I=1,20)
         END IF

         Write(STLNE2,'(A27,I3)')'Solve CASPT2 eqs for state ',
     &                               MSTATE(JSTATE)
         Call StatusLine('CASPT2:',TRIM(STLNE2))
         CALL EQCTL2(ICONV)

* Save the final caspt2 energy in the global array ENERGY():
         ENERGY(JSTATE)=E2TOT

         CALL TIMING(CPTF12,CPE,TIOTF12,TIOE)
         CPUPT2=CPTF12-CPTF11
         TIOPT2=TIOTF12-TIOTF11

         IF (ICONV .NE. 0) THEN
C     No convergence. Skip the rest of the calculation.
            IRETURN = _RC_NOT_CONVERGED_
            EXIT STATELOOP
         END IF

C     Orbitals, properties:

         ! if the dens keyword is used, need accurate density and
         ! for that the serial LUSOLV file is needed, in that case copy
         ! the distributed LURHS() to LUSOLV here.
         IF(IFDENS) THEN
           CALL PCOLLVEC(IRHS,0)
           CALL PCOLLVEC(IVECX,0)
           CALL PCOLLVEC(IVECR,0)
           CALL PCOLLVEC(IVECC,1)
           CALL PCOLLVEC(IVECC2,1)
           CALL PCOLLVEC(IVECW,1)
         END IF

         IF (IFPROP) THEN
           IF (IPRGLB.GE.USUAL) THEN
             WRITE(6,*)
             WRITE(6,'(20A4)')('****',I=1,20)
             WRITE(6,*)' CASPT2 PROPERTY SECTION'
           END IF
           CALL PRPCTL
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

C     Gradients.
C     Note: Quantities computed in gradients section can also
C     be used efficiently for computing Multi-State HEFF.
C     NOTE: atm the MS-CASPT2 couplings computed here are wrong!
         IF(IFDENS) THEN
           IF (IPRGLB.GE.VERBOSE) THEN
              WRITE(6,*)
              WRITE(6,'(20A4)')('****',I=1,20)
              IF(NSTATE.GT.1) THEN
              WRITE(6,*)' CASPT2 GRADIENT/MULTI-STATE COUPLINGS SECTION'
              ELSE
                 WRITE(6,*)' CASPT2 GRADIENT SECTION'
              END IF
           END IF
           Call StatusLine('CASPT2:','Multi-State couplings')
C-SVC: for now, this part is only performed on the master node
#ifdef _MOLCAS_MPP_
           IF (Is_Real_Par()) THEN
             Call Set_Do_Parallel(.False.)
             IF (KING()) CALL GRDCTL(HEFF)
             Call Set_Do_Parallel(.True.)
             CALL GASync
           ELSE
             CALL GRDCTL(HEFF)
           END IF
#else
           CALL GRDCTL(HEFF)
#endif
         END IF

         IF((.NOT.IFDENS) .AND. IFMSCOUP) THEN
C     If this was NOT a gradient, calculation, then the multi-state
C     couplings are more efficiently computed via three-body
C     transition density matrices.
           IF (IPRGLB.GE.VERBOSE) THEN
              WRITE(6,*)
              WRITE(6,'(20A4)')('****',I=1,20)
              WRITE(6,*)' CASPT2 MULTI-STATE COUPLINGS SECTION'
           END IF
           Call StatusLine('CASPT2:','Multi-State couplings')
           CALL MCCTL(HEFF)
         END IF


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
          WRITE(6,'(A,I6)')    '  CASPT2 TIMING INFO FOR STATE ',JSTATE
          WRITE(6,*)
          WRITE(6,'(A)')       '                        '//
     &                         ' cpu time  (s) '//
     &                         ' wall time (s) '
          WRITE(6,'(A)')       '                        '//
     &                         ' ------------- '//
     &                         ' ------------- '
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
          WRITE(6,'(A,2F14.2)')'  Gradient/MS coupling  ',CPUGRD,TIOGRD
          WRITE(6,'(A,2F14.2)')' Total time             ',CPUTOT,TIOTOT
          WRITE(6,*)
        END IF

C End of long loop over states in the group
       END DO
       IF (IPRGLB.GE.USUAL) THEN
        If(.not.IFNOPT2) Then
         CALL CollapseOutput(0,'CASPT2 computation for group ')
         WRITE(6,*)
        End If
       END IF
C End of long loop over groups
        JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
9999    write (6,*)
      END DO STATELOOP

1000  CONTINUE

      IF (IRETURN.NE.0) GOTO 9000
       If(IFNOPT2) then  !XMS Skip multistate calculation.
        write(6,*)'PT2 calculation skipped with XROH keyword'
        write(6,*)
        CALL MMA_DEALLOCATE(UEFF)
        CALL MMA_DEALLOCATE(U0)
       Else
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,*)' Total CASPT2 energies:'
       DO I=1,NSTATE
        IF ((NLYROOT.NE.0).AND.(I.NE.NLYROOT)) CYCLE
        CALL PrintResult(6,'(6x,A,I3,5X,A,F16.8)',
     &    'CASPT2 Root',I,'Total energy:',ENERGY(I),1)
       END DO
       WRITE(6,*)
       IF (IFXMS) THEN
        WRITE(6,*)' Note that these CASPT2 energies are obtained using'
        WRITE(6,*)' the XMS Fock operator and thus do not correspond'
        WRITE(6,*)' to the true single-state CASPT2 ones.'
       END IF
       WRITE(6,*)
      END IF
      IF(IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
       WRITE(6,*)' Relative CASPT2 energies:'
       WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &   'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
       DO I=1,NSTATE
        RELAU = ENERGY(I)-ENERGY(1)
        RELEV = RELAU * CONV_AU_TO_EV_
        RELCM = RELAU * CONV_AU_TO_CM1_
        RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
        WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &   I, RELAU, RELEV, RELCM, RELKJ
       END DO
       WRITE(6,*)
      END IF

      IF(NLYROOT.NE.0) IFMSCOUP=.FALSE.
      IF(IFMSCOUP) THEN
        Call StatusLine('CASPT2:','Effective Hamiltonian')
        CALL MLTCTL(HEFF,UEFF,U0)
      END IF

      IF(IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
       WRITE(6,*)' Relative (X)MS-CASPT2 energies:'
       WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &   'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
       DO I=1,NSTATE
        RELAU = ENERGY(I)-ENERGY(1)
        RELEV = RELAU * CONV_AU_TO_EV_
        RELCM = RELAU * CONV_AU_TO_CM1_
        RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
        WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &   I, RELAU, RELEV, RELCM, RELKJ
       END DO
       WRITE(6,*)
      END IF

* create a JobMix file
* (note that when using HDF5 for the PT2 wavefunction, IFMIX is false)
      CALL CREIPH_CASPT2(Heff,Ueff,U0)

* Store the PT2 energy and effective hamiltonian on the wavefunction file
      CALL PT2WFN_ESTORE(HEFF)

* store information on runfile for geometry optimizations
      Call Put_iScalar('NumGradRoot',iRlxRoot)
      Call Store_Energies(NSTATE,ENERGY,iRlxRoot)

      CALL MMA_DEALLOCATE(UEFF)
      CALL MMA_DEALLOCATE(U0)
      End If  !Skipping MultiState calculation when IFNOPT2=true
9000  CONTINUE

C Free resources, close files
      CALL PT2CLS

      CALL MMA_DEALLOCATE(HEFF)
      CALL MMA_DEALLOCATE(H0)

C     PRINT I/O AND SUBROUTINE CALL STATISTICS
      IF ( IPRGLB.GE.USUAL ) THEN
        CALL FASTIO('STATUS')
        CALL QSTAT(' ')
      END IF

      Call StatusLine('CASPT2:','Finished.')
      CALL QEXIT('CASPT2')
      RETURN
      END
