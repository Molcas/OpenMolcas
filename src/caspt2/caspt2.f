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
      USE INPUTDATA, ONLY: INPUT
      USE PT2WFN
      use fciqmc_interface, only: DoFCIQMC
      use caspt2_output, only: iPrGlb, terse, usual, verbose
      use caspt2_gradient, only: do_grad, nStpGrd, iStpGrd, IDSAVGRD
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King, Set_Do_Parallel
#endif
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
#include "warnings.h"
#include "constants.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "stdalloc.fh"
#include "caspt2_grad.fh"
      CHARACTER(len=60) STLNE2
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

* Gradient stuff
      REAL*8, ALLOCATABLE :: UeffSav(:,:),U0Sav(:,:),H0Sav(:,:),ESav(:)
      LOGICAL :: IFGRDT0 = .False.

      Call StatusLine('CASPT2:','Just starting')

      IRETURN = 0

      CALL SETTIM

* Probe the environment to globally set the IPRGLB value
      Call Set_Print_Level

*======================================================================*
*
      Call StatusLine('CASPT2:','Initializing')
      CALL PT2INI
* Initialize effective Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(Heff,Nstate,Nstate,Label='Heff')
      CALL MMA_ALLOCATE(Ueff,Nstate,Nstate,Label='Ueff')
      Heff=0.0D0
      Ueff=0.0D0
* Initialize zeroth-order Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(H0,Nstate,Nstate,Label='H0')
      CALL MMA_ALLOCATE(U0,Nstate,Nstate,Label='U0')
      H0(:,:)=0.0D0
* U0 is initialized as the identity matrix, in the case of a
* standard MS-CASPT2 calculation it will not be touched anymore
      U0(:,:)=0.0D0
      call dcopy_(Nstate,[1.0d0],0,U0,Nstate+1)
*
* Some preparations for gradient calculation
      IF (do_grad) Then
        CALL MMA_ALLOCATE(UeffSav,Nstate,Nstate)
        CALL MMA_ALLOCATE(U0Sav,Nstate,Nstate)
        IDSAVGRD = 0
      End If
*======================================================================*
* Put the CASSCF energies on the diagonal of Heff, i.e. form the
* first-order corrected effective Hamiltonian:
*     Heff[1] = PHP
* and later on we will add the second-order correction
* Heff(2) = PH \Omega_1 P to Heff[1]
      DO I=1,NSTATE
        HEFF(I,I) = REFENE(I)
      END DO
      IF (IPRGLB.GE.VERBOSE) THEN
        write(6,*)' Heff[1] in the original model space basis:'
        call prettyprint(Heff,Nstate,Nstate)
      END IF
* If the EFFE keyword has been used, we already have the multi state
* coupling Hamiltonian effective matrix, just copy the energies and
* proceed to the MS coupling section.
      IF (INPUT%JMS) THEN
        ! in case of XMS, XDW, RMS, we need to rotate the states
        if (IFXMS .or. IFRMS) then
          call xdwinit(Heff,H0,U0)
        end if
        DO I=1,NSTATE
          ENERGY(I)=INPUT%HEFF(I,I)
        END DO
        HEFF(:,:)=INPUT%HEFF(:,:)
        GOTO 1000
      END IF

* In case of a XDW-CASPT2 calculation we first rotate the CASSCF
* states according to the XMS prescription in xdwinit
      if ((IFXMS .and. IFDW) .or. (IFRMS)) then
        call xdwinit(Heff,H0,U0)
        call wgtinit(Heff)
      else
        call wgtinit(Heff)
      end if

* Before entering the long loop over groups and states, precompute
* the 1-RDMs for all states and mix them according to the type of
* calculation: MS, XMS, DW or XDW.
      call rdminit

      !! loop for multistate CASPT2 gradient
      !! In the first step, the effective Hamiltonian and the
      !! rotation vector are computed. In the second step,
      !! quantities needed for gradient are computed
      nStpGrd = 1
      IFGRDT0 = do_grad
      ! Why do we need to do this here?
      If (do_grad.AND.IfChol) CALL CNSTFIFAFIMO(0)
      IF (do_grad.AND.IFMSCOUP) Then
        nStpGrd = 2
        !! avoid doing a lot of calculations in dens.f in the first loop
        do_grad = .false.
        CALL MMA_ALLOCATE(ESav,Nstate)
        CALL MMA_ALLOCATE(H0Sav,Nstate,Nstate)
        If ((IFXMS.and.IFDW) .or. IFRMS) Then
          ! at this stage, H0 is just a zero matrix
          Call DCopy_(Nstate*Nstate,H0,1,H0Sav,1)
        End IF
      End If

      !! iStpGrd = 0 means PT2GRD is found and program requests gradient
      !! for the second root or NAC vectors. Skip the energy calculation
      !! and directly goes to the gradient part
      If (iStpGrd == 0) Then
        Call SavGradParams2(2,UEFF,U0,H0)
        Call DCopy_(Nstate,ENERGY,1,Esav,1)
        Call DCopy_(Nstate*Nstate,UEFF,1,UEFFSav,1)
        Call DCopy_(Nstate*Nstate,U0,1,U0Sav,1)
        If ((IFXMS .and. IFDW) .or. IFRMS)
     *    Call DCopy_(Nstate*Nstate,H0,1,H0Sav,1)
        iStpGrd = 2
C       Call EQCTL2(ICONV)
        Go To 8999
      End If

* FIRST GRAD LOOP ITER

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
         WRITE(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
         CALL CollapseOutput(1,TRIM(STLNE2))
         WRITE(6,*)
       END IF

       CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
       CALL GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0)
       CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
       CPUGIN=CPTF10-CPTF0
       TIOGIN=TIOTF10-TIOTF0

       If (do_grad) CALL CNSTFIFAFIMO(1)

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
* No convergence. Skip the rest of the calculation.
           IRETURN = _RC_NOT_CONVERGED_
           EXIT STATELOOP
         END IF

* Orbitals, properties:

         ! if the dens keyword is used, need accurate density and
         ! for that the serial LUSOLV file is needed, in that case copy
         ! the distributed LURHS() to LUSOLV here.
         IF (IFDENS.OR.(do_grad.and.(iRlxRoot.eq.MSTATE(JSTATE)))) THEN
           CALL PCOLLVEC(IRHS,0)
           CALL PCOLLVEC(IVECX,0)
           CALL PCOLLVEC(IVECR,0)
           CALL PCOLLVEC(IVECC,1)
           CALL PCOLLVEC(IVECC2,1)
           CALL PCOLLVEC(IVECW,1)
         END IF

         IF (IFPROP.OR.(do_grad.and.(iRlxRoot.eq.MSTATE(JSTATE)))) THEN
           IF (IPRGLB.GE.USUAL) THEN
             WRITE(6,*)
             WRITE(6,'(20A4)')('****',I=1,20)
             WRITE(6,*)' CASPT2 PROPERTY SECTION'
           END IF
           CALL PRPCTL(0,UEFF,U0)
         ELSE
           IF (IPRGLB.GE.USUAL) THEN
             WRITE(6,*)
             WRITE(6,*)'  (Skipping property calculation,'
             WRITE(6,*)'   use PROP keyword to activate)'
           END IF
         END IF
         !! Save many quantities needed for gradient calculations
         If (nStpGrd == 2) Call SavGradParams(1,IDSAVGRD)

         CALL TIMING(CPTF13,CPE,TIOTF13,TIOE)
         CPUPRP=CPTF13-CPTF12
         TIOPRP=TIOTF13-TIOTF12

        if (.not. DoFCIQMC) then
* Gradients.
* Note: Quantities computed in gradients section can also
* be used efficiently for computing Multi-State HEFF.
* NOTE: atm the MS-CASPT2 couplings computed here are wrong!
         IF (IFDENS.AND..NOT.IFGRDT0) THEN
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
* SVC: for now, this part is only performed on the master node
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

         IF ((.NOT.IFDENS.OR.IFGRDT0).AND.IFMSCOUP) THEN
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
        end if

        IF (IPRGLB.GE.VERBOSE) THEN
          WRITE(6,*)
          WRITE(6,'(A,I6)')    '  CASPT2 TIMING INFO FOR STATE ',
     &                         MSTATE(JSTATE)
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
          if (.not. DoFCIQMC) then ! MS-CASPT2 currently not possible
           WRITE(6,'(A,2F14.2)')'  MS coupling           ',CPUGRD,TIOGRD
          end if
          WRITE(6,'(A,2F14.2)')' Total time             ',CPUTOT,TIOTOT
          WRITE(6,*)
        END IF

* End of long loop over states in the group
       END DO
       IF (IPRGLB.GE.USUAL) THEN
         CALL CollapseOutput(0,'CASPT2 computation for group ')
         WRITE(6,*)
       END IF
* End of long loop over groups
        JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
      END DO STATELOOP

1000  CONTINUE

      IF (IRETURN.NE.0) GOTO 9000
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,*)' Total CASPT2 energies:'
       IF (IFXMS.or.IFRMS) THEN
        WRITE(6,*)
        WRITE(6,*)' State-specific CASPT2 energies obtained using'
        WRITE(6,*)' rotated states do not have a well-defined physical'
        WRITE(6,*)' meaning, however they can be extracted from the'
        WRITE(6,*)' diagonal of the effective Hamiltonian.'
       ELSE
        DO I=1,NSTATE
         IF ((NLYROOT.NE.0).AND.(I.NE.NLYROOT)) CYCLE
         CALL PrintResult(6,'(6x,A,I3,5X,A,F16.8)',
     &    'CASPT2 Root',MSTATE(I),'Total energy:',ENERGY(I),1)
        END DO
        IF(IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
         WRITE(6,*)
         WRITE(6,*)' Relative CASPT2 energies:'
         WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &     'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
         ISTATE=1
         DO I=2,NSTATE
           IF (ENERGY(I).LT.ENERGY(ISTATE)) ISTATE=I
         END DO
         DO I=1,NSTATE
          RELAU = ENERGY(I)-ENERGY(ISTATE)
          RELEV = RELAU * CONV_AU_TO_EV_
          RELCM = RELAU * CONV_AU_TO_CM1_
          RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
          WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &     MSTATE(I), RELAU, RELEV, RELCM, RELKJ
         END DO
        END IF
       END IF
       WRITE(6,*)
      END IF

 8999 Continue

      if (.not. doFCIQMC) then
        if (iStpGrd.ne.2) then
          IF (NLYROOT.NE.0) IFMSCOUP=.FALSE.
          IF (IFMSCOUP) THEN
            Call StatusLine('CASPT2:','Effective Hamiltonian')
            CALL MLTCTL(HEFF,UEFF,U0)

            IF (IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
             WRITE(6,*)' Relative (X)MS-CASPT2 energies:'
             WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &         'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
             DO I=1,NSTATE
              RELAU = ENERGY(I)-ENERGY(1)
              RELEV = RELAU * CONV_AU_TO_EV_
              RELCM = RELAU * CONV_AU_TO_CM1_
              RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
              WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &         I, RELAU, RELEV, RELCM, RELKJ
             END DO
             WRITE(6,*)
            END IF
          END IF
        end if

! Beginning of second step, in case gradient of (X)MS

        If (nStpGrd == 2) Then
          IF (IPRGLB.GE.TERSE) THEN
            WRITE(6,'(20A4)')('****',I=1,20)
            WRITE(6,'(A)')
     &      ' SECOND RUN to compute analytical gradients/NAC quantities'
            WRITE(6,'(20A4)')('----',I=1,20)
            WRITE(6,*)
            CALL XFlush(6)
          END IF
          do_grad = .true.
          Call DCopy_(nState,ENERGY,1,Esav,1)
          Call DCopy_(nState**2,Ueff,1,UeffSav,1)
          If (IFXMS .or. IFRMS) Call DCopy_(nState**2,U0,1,U0Sav,1)
!
          !!Somehow H0 is wrong for XDW-CASPT2
          !!Maybe, H0(1,1) is computed with rotated basis with
          !!DW-density, while the true value is computed with SA-density
          If (do_grad .AND. IFMSCOUP) Then
            If ((IFXMS .and. IFDW) .or. IFRMS) Then
              Call DCopy_(nState*nState,H0Sav,1,H0,1)
            End If
          End If
          Call SavGradParams2(1,UEFF,U0,H0)
!
          Call GradLoop(Heff,Ueff,H0,U0,H0Sav)
        End If

* Back-transform the effective Hamiltonian and the transformation matrix
* to the basis of original CASSCF states
        If (nStpGrd == 2) Then
          CALL Backtransform(Heff,UeffSav,U0sav)
          Call DCopy_(nState*nState,UeffSav,1,Ueff,1)
          Call DCopy_(nState,ESav,1,ENERGY,1)
        Else
          CALL Backtransform(Heff,Ueff,U0)
        End If

* create a JobMix file
* (note that when using HDF5 for the PT2 wavefunction, IFMIX is false)
        CALL CREIPH_CASPT2(Heff,Ueff,U0)

* Store the PT2 energy and effective Hamiltonian on the wavefunction file
        CALL PT2WFN_ESTORE(HEFF)

* Store rotated states if XMUL + NOMUL
        IF ((IFXMS .or. IFRMS) .AND. (.NOT.IFMSCOUP)) CALL PT2WFN_DATA

* store information on runfile for geometry optimizations
        Call Put_iScalar('NumGradRoot',iRlxRoot)
        Call Store_Energies(NSTATE,ENERGY,iRlxRoot)
      end if


9000  CONTINUE

      !! Finishing for gradient calculation
      IF (IFGRDT0) Then
        Call GrdCls(IRETURN,UEFFSav,U0Sav,H0)
        CALL MMA_DEALLOCATE(UeffSav)
        CALL MMA_DEALLOCATE(U0Sav)
        IF (IFMSCOUP) Then
          CALL MMA_DEALLOCATE(ESav)
          CALL MMA_DEALLOCATE(H0Sav)
        End If
      End If
C Free resources, close files
      CALL PT2CLS()

      CALL MMA_DEALLOCATE(UEFF)
      CALL MMA_DEALLOCATE(U0)
      CALL MMA_DEALLOCATE(HEFF)
      CALL MMA_DEALLOCATE(H0)

C     PRINT I/O AND SUBROUTINE CALL STATISTICS
      IF ( IPRGLB.GE.USUAL ) THEN
        CALL FASTIO('STATUS')
      END IF

      Call StatusLine('CASPT2:','Finished.')
      RETURN
      END
