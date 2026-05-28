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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!               2019, Stefano Battaglia                                *
!***********************************************************************

subroutine CASPT2(IRETURN)
!----------------------------------------------------------------------*
!     1998  PER-AAKE MALMQUIST                                         *
!     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
!     UNIVERSITY OF LUND, SWEDEN                                       *
!----------------------------------------------------------------------*

! SECOND ORDER PERTURBATION CALCULATIONS WITH A CASSCF
! REFERENCE FUNCTION.

! ORIGINAL CASPT2 PROGRAM WRITTEN 890526 BY:
! KERSTIN ANDERSSON (ALL CASPT2 CODE WITH FOLLOWING EXCEPTIONS)
! PER-AKE MALMQVIST (THE GUGA PART)
! BJORN O ROOS      (THE INTEGRAL TRANSFORMATION)

! MODIFIED 1991-02-23 BY PER-AAKE MALMQUIST, FOR USE WITH THE
! MOLCAS RASSCF PROGRAM.

! ALL PROGRAM REWRITTEN (MALMQVIST 1993) FOR MOLCAS VERSION 3,
! EXCEPT TRACTL AND TRA2 KEPT BUT SLIGHTLY MODIFIED, AND ALSO:
! INTERFACING WITH MOLCAS-3 BY MARCUS FUELSCHER

! REWRITTEN FOR 1) EXACT ACTIVE DENSITY MATRIX
! 2) QUANTITIES NEEDED FOR GRADIENTS
! 3) MULTI-STATE CASPT2
! BY MALMQVIST 1998.

! SINCE THEN, THE FOLLOWING MODIFICATIONS HAVE BEEN MADE:
! IPEA HAMILTONIAN BY G. GHIGO & P. MALMQVIST
!   Chem. Phys. Lett. 396, pp 142 (2004)
!   http://dx.doi.org/10.1016/j.cplett.2004.08.032
! LOVCASPT2/FNOCASPT2/AFREEZE BY B. ROOS & F. AQUILANTE
!   J. Chem. Phys. 131, 034113 (2009)
!   http://dx.doi.org/10.1063/1.3157463
! CHOLESKY SUPPORT BY P. MALMQVIST AND F. AQUILANTE
!   J. Chem. Theory Comput. 4 (5), pp 694-702 (2008)
!   http://dx/doi.org/10.1021/ct700263h
! RASPT2 BY P. MALMQVIST
!   J. Chem. Phys. 128, 204109 (2008)
!   http://dx.doi.org/10.1063/1.2920188
! PARALLELIZATION BY S. VANCOILLIE
!   J. Comp. Chem. 34, pp 1937-1948 (2013)
!   http://dx.doi.org/10.1002/jcc.23342
! MAJOR REFACTORING OF INITIALIZATION PHASE AND ADDITION OF
! HDF5 SUPPORT BY S. VANCOILLIE (2016)

!***********************************************************************

use INPUTDATA, only: INPUT
use PT2WFN, only: PT2WFN_DATA, PT2WFN_ESTORE
use fciqmc_interface, only: DoFCIQMC
use caspt2_global, only: do_grad, IDSAVGRD, iPrGlb, iStpGrd, nStpGrd
use caspt2_module, only: CPUEIG, CPUFG3, CPUFMB, CPUGIN, CPUGRD, CPUINT, CPULCS, CPUNAD, CPUOVL, CPUPCG, CPUPRP, CPUPT2, CPURHS, &
                         CPUSBM, CPUSCA, CPUSER, CPUSGM, CPUSIN, CPUVEC, E2ToT, Energy, IfChol, IfDens, IfDW, IfMSCoup, IfProp, &
                         IfRMS, IfXMS, iRlxRoot, jState, mState, nGroup, nGroupState, nLyGroup, nLyRoot, nState, RefEne, TIOEIG, &
                         TIOFG3, TIOFMB, TIOGIN, TIOGRD, TIOINT, TIOLCS, TIONAD, TIOOVL, TIOPCG, TIOPRP, TIOPT2, TIORHS, TIOSBM, &
                         TIOSCA, TIOSER, TIOSGM, TIOSIN, TIOVEC
use PrintLevel, only: TERSE, USUAL, VERBOSE
use EQSOLV, only: iRHS, iVecC, iVecC2, iVecR, iVecW, iVecX
#ifdef _DMRG_
use, intrinsic :: iso_c_binding, only: c_bool, c_int
use qcmaquis_interface, only: qcmaquis_interface_compute_and_store_123rdm_full, &
                              qcmaquis_interface_compute_and_store_trans_123rdm_full
use caspt2_module, only: DMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auTocm, auToeV, auTokJmol
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: IRETURN
integer(kind=iwp) :: ICONV, IGROUP, ISTATE, JSTATE_OFF
#ifdef _DMRG_
integer(kind=iwp) :: I, J
#endif
real(kind=wp) :: CPE, CPTF0, CPTF11, CPTF12, CPTF13, CPTF14, CPUTOT, RELAU, RELCM, RELEV, RELKJ, TIOE, TIOTF0, TIOTF11, TIOTF12, &
                 TIOTF13, TIOTF14, TIOTOT
logical(kind=iwp) :: IFGRDT0
character(len=60) :: STLNE2
real(kind=wp), allocatable :: ESav(:), H0(:,:), H0Sav(:,:), Heff(:,:), U0(:,:), U0Sav(:,:), Ueff(:,:), UeffSav(:,:)
#include "warnings.h"

call StatusLine('CASPT2: ','Just starting')

IRETURN = 0
IFGRDT0 = .false.

!=======================================================================
!
! Miscellaneous setup, reading of the input and start of writing
! input parameters in the log file.
call PT2INI()

! Set up effective Hessian for multi-state calculations. Note that
! single state calculations are treated as an one-dimensional case.
call HEFF_INI()

! If the EFFE keyword has been used proceed to the MS coupling section.
if (INPUT%JMS) then
  call Print_Truff()
  call Post_Process()
  call CASPT2_TERM()
  return
end if

! Before entering the long loop over groups and states, precompute
! the 1-RDMs for all states and mix them according to the type of
! calculation: MS, XMS, DW or XDW. This is done in the natural
! orbital basis.

call rdminit()

!! loop for multistate CASPT2 gradient
!! In the first step, the effective Hamiltonian and the
!! rotation vector are computed. In the second step,
!! quantities needed for gradient are computed
nStpGrd = 1
IFGRDT0 = do_grad
! Why do we need to do this here?
if (do_grad .and. IfChol) call CNSTFIFAFIMO(0)
if (do_grad .and. IFMSCOUP) then
  nStpGrd = 2
  !! avoid doing a lot of calculations in dens in the first loop
  do_grad = .false.
  call MMA_ALLOCATE(ESav,Nstate)
  call MMA_ALLOCATE(H0Sav,Nstate,Nstate)
  ! at this stage, H0 is just a zero matrix
  if ((IFXMS .and. IFDW) .or. IFRMS) H0Sav(:,:) = H0(:,:)
end if

!! iStpGrd = 0 means PT2GRD is found and program requests gradient
!! for the second root or NAC vectors. Skip the energy calculation
!! and directly goes to the gradient part
if (iStpGrd == 0) then
  call SavGradParams2(2,UEFF,U0,H0,Nstate)
  Esav(:) = ENERGY(1:Nstate)
  UEFFSav(:,:) = UEFF(:,:)
  U0Sav(:,:) = U0(:,:)
  if ((IFXMS .and. IFDW) .or. IFRMS) H0Sav(:,:) = H0(:,:)
  iStpGrd = 2
  call Post_Process()
  call CASPT2_TERM()
  return
end if

! FIRST GRAD LOOP ITER
#ifdef _DMRG_
if (DMRG) then
  do I=1,NSTATE
    call qcmaquis_interface_compute_and_store_123rdm_full(int(I-1,c_int),logical(.true.,c_bool))
    do J=1,NSTATE
      if (I /= J) call qcmaquis_interface_compute_and_store_trans_123rdm_full(int(I-1,c_int),int(J-1,c_int),logical(.true.,c_bool))
    end do
  end do
end if
#endif

!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
! For (X)Multi-State, a long loop over root states.
! The states are ordered by group, with each group containing a number
! of group states for which GRPINI is called.
JSTATE_OFF = 0
STATELOOP: do IGROUP=1,NGROUP

  call TIMING(CPTF0,CPE,TIOTF0,TIOE)

  if ((NLYGROUP /= 0) .and. (IGROUP /= NLYGROUP)) then
    JSTATE_OFF = JSTATE_OFF+NGROUPSTATE(IGROUP)
    cycle
  end if

  if (IPRGLB >= USUAL) then
    write(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
    call CollapseOutput(1,trim(STLNE2))
    write(u6,*)
  end if

  ! GRPINI (group init?) does a number of things as follows (note this
  ! list is not complete).

  ! 1. Reads the natural orbitals (NO) MOs from LUONEM. Stored in
  !    CMO.
  ! 2. transforms HONE from AO to NO basis (traone).
  ! 3. transforms the two-electron integrals to the NO basis
  ! 4. Form the 1-particle density matrix (DREF) from DMIX as
  !    generated RDMINIT above. In the NO basis. This is needed
  !    to form the Fock matrix in the NO basis.
  ! 5. Forms the Fock matrix in NO basis by a call to MkOrb.
  !    For conventional integrals for it directly in the NO basis,
  !    while for the Cholesky decomposition approach the Fock matrix
  !    is first formed in the AO basis and then transformed to NO
  !    basis. This routine also reads the CI coefficients in the
  !    NO basis from the file LUCIEX. Transforms the to PCO basis and
  !    writes them back to LUCIEX on another place.
  ! 6. Forms the pseudo canonical orbitals (PCO), in OrbCtl.
  !    6a. For the PCO basis and the transformation matrix between
  !        the NO and the PCO basis, the TOrb matrix.
  !    6b. Transform the FIMO matrix to the PCO basis
  !    6c. Transform the FIFA matrix to the PCO basis
  ! 7. Transforms the Cholesky vectors or the two-electron integrals
  !    to the PCO basis.
  ! 8. Drops the NO CMO orbitals as stored in CMO.

  call GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0,nState)

  if (do_grad) call CNSTFIFAFIMO(1)

  do ISTATE=1,NGROUPSTATE(IGROUP)
    JSTATE = JSTATE_OFF+ISTATE

    ! Skip this state if we only need 1 state and it isn't this "one".
    if ((NLYROOT /= 0) .and. (JSTATE /= NLYROOT)) cycle

    ! STINI (state init?) does the following

    ! 1. It calls POLY3 to form the 1-, 2-, and 3-particle density
    !    matrix in the PCO basis. This is based on the CI vector
    !    expressed in the PCO basis. Using pt2_put the result is
    !    stored on LUDMAT.
    ! 2. Using GETDPREF it pulls the one- and two-particle density
    !    matrices (gamma 1 and gamma 2), using pt2_get, of the
    !    LUDMAT file and sticks them in to DREF and PREF.
    ! 3. Sets the EREF value
    ! 4. Recomputes EASUM.

    call STINI(JSTATE)

    ! Solve CASPT2 equation system and compute corr energies.
    if (IPRGLB >= USUAL) then
      write(u6,'(A)') repeat('*',80)
      write(u6,*) ' CASPT2 EQUATION SOLUTION'
      write(u6,'(A)') repeat('-',80)
    end if

    write(STLNE2,'(A,I0)') 'Solve CASPT2 eqs. for state ',MSTATE(JSTATE)

    call TIMING(CPTF11,CPE,TIOTF11,TIOE)

    call StatusLine('CASPT2: ',STLNE2)

    call EQCTL2(ICONV)

    ! Save the final caspt2 energy in the global array ENERGY():
    ENERGY(JSTATE) = E2TOT

    call TIMING(CPTF12,CPE,TIOTF12,TIOE)
    CPUPT2 = CPTF12-CPTF11
    TIOPT2 = TIOTF12-TIOTF11

    if (ICONV /= 0) then
      ! No convergence. Skip the rest of the calculation.
      IRETURN = _RC_NOT_CONVERGED_
      exit STATELOOP
    end if

    ! Orbitals, properties:

    ! if the dens keyword is used, need accurate density and
    ! for that the serial LUSOLV file is needed, in that case copy
    ! the distributed LURHS() to LUSOLV here.
    if (IFDENS .or. (do_grad .and. (iRlxRoot == MSTATE(JSTATE)))) then
      call PCOLLVEC(IRHS,0)
      call PCOLLVEC(IVECX,0)
      call PCOLLVEC(IVECR,0)
      call PCOLLVEC(IVECC,1)
      call PCOLLVEC(IVECC2,1)
      call PCOLLVEC(IVECW,1)
    end if

    if (IFPROP .or. (do_grad .and. (iRlxRoot == MSTATE(JSTATE)))) then

      call PRPCTL(0,UEFF,U0,nState)

    else
      if (IPRGLB >= USUAL) then
        write(u6,*)
        write(u6,*) '  (Skipping property calculation,'
        write(u6,*) '   use PROP keyword to activate)'
      end if
    end if
    !! Save many quantities needed for gradient calculations
    if (nStpGrd == 2) call SavGradParams(1,IDSAVGRD)

    call TIMING(CPTF13,CPE,TIOTF13,TIOE)
    CPUPRP = CPTF13-CPTF12
    TIOPRP = TIOTF13-TIOTF12

    if (.not. DoFCIQMC) then

      if (IFMSCOUP) then
        ! If this was NOT a gradient, calculation, then the multi-state
        ! couplings are more efficiently computed via three-body
        ! transition density matrices.
        if (IPRGLB >= VERBOSE) then
          write(u6,*)
          write(u6,'(A)') repeat('*',80)
          write(u6,*) ' CASPT2 MULTI-STATE COUPLINGS SECTION'
        end if
        call StatusLine('CASPT2: ','Multi-State couplings')
        call MCCTL(HEFF,NSTATE,JSTATE)
      end if

    end if

    call Iter_Timing()

    ! End of long loop over states in the group
  end do

  if (IPRGLB >= USUAL) then
    call CollapseOutput(0,'CASPT2 computation for group ')
    write(u6,*)
  end if

  ! End of long loop over groups
  JSTATE_OFF = JSTATE_OFF+NGROUPSTATE(IGROUP)
end do STATELOOP

!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************

call Print_Truff()
call Post_Process()
call CASPT2_TERM()

contains

subroutine Print_Truff()

  integer(kind=iwp) :: I

  if (IRETURN /= 0) then
    call CASPT2_TERM()
    return
  end if

  if (IPRGLB >= TERSE) then
    write(u6,*) ' Total CASPT2 energies:'
    if (IFXMS .or. IFRMS) then
      write(u6,*)
      write(u6,*) ' State-specific CASPT2 energies obtained using'
      write(u6,*) ' rotated states do not have a well-defined physical'
      write(u6,*) ' meaning, however they can be extracted from the'
      write(u6,*) ' diagonal of the effective Hamiltonian.'
    else
      do I=1,NSTATE
        if ((NLYROOT /= 0) .and. (I /= NLYROOT)) cycle
        call PrintResult(u6,'(6x,A,I3,5X,A,F16.8)','CASPT2 Root',MSTATE(I),'Total energy:',ENERGY(I),1)
      end do
      if ((IPRGLB >= VERBOSE) .and. (NLYROOT == 0)) then
        write(u6,*)
        write(u6,*) ' Relative CASPT2 energies:'
        write(u6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)') 'Root','(a.u.)','(eV)','(cm^-1)','(kJ/mol)'
        ISTATE = minloc(ENERGY(1:NSTATE),dim=1)
        do I=1,NSTATE
          RELAU = ENERGY(I)-ENERGY(ISTATE)
          RELEV = RELAU*auToeV
          RELCM = RELAU*auTocm
          RELKJ = RELAU*auTokJmol
          write(u6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)') MSTATE(I),RELAU,RELEV,RELCM,RELKJ
        end do
      end if
    end if
    write(6,*)
  end if

end subroutine Print_Truff

subroutine Post_Process()

  integer(kind=iwp) :: I

  if (.not. doFCIQMC) then
    if (iStpGrd /= 2) then
      if (NLYROOT /= 0) IFMSCOUP = .false.
      if (IFMSCOUP) then
        call StatusLine('CASPT2: ','Effective Hamiltonian')
        call MLTCTL(HEFF,UEFF,U0)

        if ((IPRGLB >= VERBOSE) .and. (NLYROOT == 0)) then
          write(u6,*) ' Relative (X)MS-CASPT2 energies:'
          write(u6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)') 'Root','(a.u.)','(eV)','(cm^-1)','(kJ/mol)'
          do I=1,NSTATE
            RELAU = ENERGY(I)-ENERGY(1)
            RELEV = RELAU*auToeV
            RELCM = RELAU*auTocm
            RELKJ = RELAU*auTokJmol
            write(u6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)') I,RELAU,RELEV,RELCM,RELKJ
          end do
          write(u6,*)
        end if
      end if
    end if

    ! Beginning of second step, in case gradient of (X)MS

    if (nStpGrd == 2) then
      if (IPRGLB >= TERSE) then
        write(u6,'(A)') repeat('*',80)
        write(u6,'(A)') ' SECOND RUN to compute analytical gradients/NAC quantities'
        write(u6,'(A)') repeat('-',80)
        write(u6,*)
        call XFlush(u6)
      end if
      do_grad = .true.
      Esav(:) = ENERGY(1:nState)
      UeffSav(:,:) = Ueff(:,:)
      if (IFXMS .or. IFRMS) U0Sav(:,:) = U0(:,:)

      !! Somehow H0 is wrong for XDW-CASPT2
      !! Maybe, H0(1,1) is computed with rotated basis with
      !! DW-density, while the true value is computed with SA-density
      if (do_grad .and. IFMSCOUP) then
        if ((IFXMS .and. IFDW) .or. IFRMS) H0(:,:) = H0Sav(:,:)
      end if
      call SavGradParams2(1,UEFF,U0,H0,nState)

      call GradLoop(Heff,Ueff,H0,U0,H0Sav,nState)
    end if

    ! Back-transform the effective Hamiltonian and the transformation matrix
    ! to the basis of original CASSCF states
    if (nStpGrd == 2) then
      call Backtransform(Heff,UeffSav,U0sav,nState)
      Ueff(:,:) = UeffSav(:,:)
      ENERGY(1:nState) = ESav(:)
    else
      call Backtransform(Heff,Ueff,U0,nState)
    end if

    ! create a JobMix file
    ! (note that when using HDF5 for the PT2 wavefunction, IFMIX is false)
    call CREIPH_CASPT2(Heff,Ueff,U0,nState)

    ! Store the PT2 energy and effective Hamiltonian on the wavefunction file
    call PT2WFN_ESTORE(HEFF,nState)

    ! Store rotated states if XMUL + NOMUL
    if ((IFXMS .or. IFRMS) .and. (.not. IFMSCOUP)) call PT2WFN_DATA()

    ! store information on runfile for geometry optimizations
    call Put_iScalar('NumGradRoot',iRlxRoot)
    call Store_Energies(NSTATE,ENERGY,iRlxRoot)
  end if

end subroutine Post_Process

subroutine CASPT2_TERM()

  !! Finishing for gradient calculation
  if (IFGRDT0) then
    call GrdCls(IRETURN,nState,UEFFSav,U0Sav,H0)
    call MMA_DEALLOCATE(UeffSav)
    call MMA_DEALLOCATE(U0Sav)
    if (IFMSCOUP) then
      call MMA_DEALLOCATE(ESav)
      call MMA_DEALLOCATE(H0Sav)
    end if
  end if

  ! Free resources, close files
  call PT2CLS()

  call MMA_DEALLOCATE(UEFF)
  call MMA_DEALLOCATE(U0)
  call MMA_DEALLOCATE(HEFF)
  call MMA_DEALLOCATE(H0)

  ! PRINT I/O AND SUBROUTINE CALL STATISTICS
  if (IPRGLB >= USUAL) call FASTIO('STATUS')

  call StatusLine('CASPT2: ','Finished.')

end subroutine CASPT2_TERM

subroutine HEFF_INI()

  integer(kind=iwp) :: I

  ! Initialize effective Hamiltonian and eigenvectors
  call MMA_ALLOCATE(Heff,Nstate,Nstate,Label='Heff')
  call MMA_ALLOCATE(Ueff,Nstate,Nstate,Label='Ueff')
  Heff(:,:) = Zero
  Ueff(:,:) = Zero
  ! Initialize zeroth-order Hamiltonian and eigenvectors
  call MMA_ALLOCATE(H0,Nstate,Nstate,Label='H0')
  call MMA_ALLOCATE(U0,Nstate,Nstate,Label='U0')
  H0(:,:) = Zero
  ! U0 is initialized as the identity matrix, in the case of a
  ! standard MS-CASPT2 calculation it will not be touched anymore
  call unitmat(U0,Nstate)

  ! Some preparations for gradient calculation
  if (do_grad) then
    call MMA_ALLOCATE(UeffSav,Nstate,Nstate)
    call MMA_ALLOCATE(U0Sav,Nstate,Nstate)
    IDSAVGRD = 0
  end if
  !=====================================================================
  ! Put the CASSCF energies on the diagonal of Heff, i.e. form the
  ! first-order corrected effective Hamiltonian:
  !     Heff[1] = PHP
  ! and later on we will add the second-order correction
  ! Heff(2) = PH \Omega_1 P to Heff[1]
  do I=1,NSTATE
    HEFF(I,I) = REFENE(I)
  end do
  if (IPRGLB >= VERBOSE) then
    write(u6,*) ' Heff[1] in the original model space basis:'
    call prettyprint(Heff,Nstate,Nstate)
  end if
  ! If the EFFE keyword has been used, we already have the multi state
  ! coupling Hamiltonian effective matrix, just copy the energies.
  if (INPUT%JMS) then
    ! in case of XMS, XDW, RMS, we need to rotate the states
    if (IFXMS .or. IFRMS) call xdwinit(Heff,H0,U0,nState)
    do I=1,NSTATE
      ENERGY(I) = INPUT%HEFF(I,I)
    end do
    HEFF(:,:) = INPUT%HEFF(:,:)
  else

    ! In case of a XDW-CASPT2 calculation we first rotate the CASSCF
    ! states according to the XMS prescription in xdwinit
    if ((IFXMS .and. IFDW) .or. IFRMS) call xdwinit(Heff,H0,U0,nState)
    call wgtinit(Heff,nState)
  end if

end subroutine HEFF_INI

subroutine Iter_Timing()

  if (.not. DoFCIQMC) then
    call TIMING(CPTF14,CPE,TIOTF14,TIOE)
    CPUGRD = CPTF14-CPTF13
    TIOGRD = TIOTF14-TIOTF13
    CPUTOT = CPTF14-CPTF0
    TIOTOT = TIOTF14-TIOTF0

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
  end if

  if (IPRGLB >= VERBOSE) then
    write(u6,*)
    write(u6,'(A,I6)') '  CASPT2 TIMING INFO FOR STATE ',MSTATE(JSTATE)
    write(u6,*)
    write(u6,'(A)') '                         cpu time  (s)  wall time (s)'
    write(u6,'(A)') '                         -------------  -------------'
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
    if (.not. DoFCIQMC) then ! MS-CASPT2 currently not possible
      write(u6,'(A,2F14.2)') '  MS coupling           ',CPUGRD,TIOGRD
    end if
    write(u6,'(A,2F14.2)') ' Total time             ',CPUTOT,TIOTOT
    write(u6,*)
  end if

end subroutine Iter_Timing

end subroutine CASPT2
