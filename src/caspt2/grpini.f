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
* Copyright (C) 2012, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE GRPINI(IGROUP,NGRP,JSTATE_OFF,HEFF,H0,U0,nState)
      use caspt2_global, only:iPrGlb
      use caspt2_global, only: CMO, CMO_Internal, FIFA, DREF, DMIX,
     &                       CMOPT2, NCMO, Weight, TORB, FIMO
      use caspt2_global, only: LUONEM
      use fciqmc_interface, only: DoFCIQMC
#ifdef _DMRG_
      use qcmaquis_interface
#endif
      use PrintLevel, only: DEBUG, USUAL, VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: CPUFMB, CPUINT, DMRG, DoCumulant,
     &                         IEOF1M, IfDW, IfsadRef, IfXMS, jState,
     &                         nConf, STSym, TIOFMB, TIOINT, mState,
     &                         iAd1m, IfChol, CPUGIN, TIOGIN, NoTri
      use pt2_guga, only: CIThr
      use Constants, only: Zero, One
      use definitions, only: iwp, wp, u6
      IMPLICIT None
      Integer(kind=iwp), intent(in):: IGROUP,NGRP,JSTATE_OFF,nState
* 2012  PER-AKE MALMQVIST
* Multi-State and XMS initialization phase
* Purpose: For a selected set IGROUP, create a set of CMO coefficients
* and orbital energies which are to be used in common for the
* calculations belonging to this set. Also, change the CI arrays in this
* group such that they diagonalize the H0 matrix.
* The states in the group can be obtained from the ordered MSTATE array,
* for which a group offset JSTATE_OFF is passed in.
#include "warnings.h"
      real(kind=wp), intent(inout):: Heff(Nstate,Nstate)
      real(kind=wp), intent(inout):: H0(Nstate,Nstate)
      real(kind=wp), intent(inout):: U0(Nstate,Nstate)

      CHARACTER(LEN=27)  STLNE2
      real(kind=wp), allocatable:: CIRef(:,:), CIXMS(:), HONE(:)
      Integer(kind=iwp) I,J,iDisk,K,iState
      Real(kind=wp) CPU1,CPU0,TIO1,TIO0,CPU,TIO
      Real(kind=wp) CPE,TIOE,CPTF0,TIOTF0,CPTF10,TIOTF10
      logical(kind=iwp) Initiate

* ---------------------------------------------------------------------
      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
* ---------------------------------------------------------------------
* Number of states in this group.
      IF (IPRGLB.EQ.DEBUG) THEN
        write(u6,*)' Entered GRPINI.'
        write(u6,*)' NSTATE=',NSTATE
        write(u6,*)' The MSTATE array:'
        write(u6,'(1x,20I4)')(MSTATE(J),J=1,NSTATE)
        write(u6,*)' IGROUP,NGRP=',IGROUP,NGRP
      END IF

      IF (NGRP.EQ.0) THEN
        WRITE(u6,*) ' Number of states in the (X)MS group is 0!'
        WRITE(u6,*) ' This should never happen, aborting...'
        CALL ABEND()
      END IF

      Write(STLNE2,'(A,I0)')'Initial phase for group ',IGROUP
      Call StatusLine('CASPT2: ',STLNE2)
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(u6,'(20A4)')('****',I=1,20)
        WRITE(u6,'(A,I3)')
     &  ' Multi-State initialization phase begins for group ',IGROUP
        WRITE(u6,'(20A4)')('----',I=1,20)
        CALL XFlush(6)
      END IF
*
* ---------------------------------------------------------------------
*
* Load the CASSCF naural orbital MO coefficients
      call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
      CMO=>CMO_Internal

      IDISK=IAD1M(1)
      call ddafile(LUONEM,2,CMO,NCMO,IDISK)

!     IEOF1M is the next free disk address on LUOneM
      IEOF1M=IDISK

!     Allocate memory for HONE, used exclusivly in MkFock.
!     Since MkFock might be called several times with the same
!     CMOs but different density matrices we controll the computation
!     of the one- and two-electron integrals from the outside. This to
!     remove redundant work.

      Call mma_allocate(HONE,NoTri,Label='HONE')
      Initiate=.TRUE.

* Loop over states, selecting those belonging to this group.
* For each such state, compute the Fock matrix in original MO basis,
* and then the zeroth-order Hamiltonian elements between states.

* Timer for Fock matrix build
      call timing(CPU0,CPU,TIO0,TIO)
* Loop over all states in group
      do J=1,Ngrp
        Jstate=J+JSTATE_OFF

* Copy the 1-RDM of Jstate from DMIX into DREF
        ! this might be obsolete if we remove sadref
        IF (IFSADREF) Then
          !! This DREF is used only for constructing the Fock in H0.
          !! DREF used in other places will be constructed in elsewhere
          !! (STINI).
          DREF(:)=Zero
          Do K = 1, Nstate
            DREF(:)=DREF(:)+Weight(K)*DMIX(:,K)
          End Do
        Else
         DREF(:)=DMIX(:,Jstate)
        End If

* Compute the Fock matrix in MO basis for state Jstate
        Call MkFock(CMO,nCMO,FIMO,SIZE(FIMO),
     &              FIFA,SIZE(FIFA),DREF,SIZE(DREF),
     &              HONE,SIZE(HONE),Initiate)

* NN.15, TODO:
* the following transformation are skipped in DMRG-CASPT2 run
* for the time, this will be fixed later to implement DMRG-MS-CASPT2
        IF (DoCumulant .or. DoFCIQMC .or. DMRG) THEN
           Call  GPRINI_FINISH()
           Return
        End If

* Loop over bra functions
        do I=1,Ngrp
          Istate=I+JSTATE_OFF
* Compute matrix element and put it into H0
          call FOPAB(FIFA,SIZE(FIFA),Istate,Jstate,H0(Istate,Jstate))
        end do

        if (IPRGLB.ge.VERBOSE) then
* In case of MS- and XDW-CASPT2 calculations, compute off-diagonal
* elements of the Fock matrix as a sanity check of the diagonal
* approximation within the generalized Bloch equation.
          if (IFDW.or.(.not.IFXMS)) then
            write(6,*)
            write(6,*) 'Fock matrix couplings'
            write(6,*) '---------------------'
            write(6,*)
            write(6,'(10X,6X,A3,I4,A3)') ' | ', MSTATE(Jstate), ' > '
            do Istate=1,Nstate
              if (Istate.ne.Jstate) then
* Compute matrix element and print it out
                call FOPAB(FIFA,SIZE(FIFA),Istate,Jstate,
     &                     H0(Istate,Jstate))
                write(6,'(A3,I4,A3,F16.8)')
     &                  ' < ',MSTATE(Istate),' | ', H0(Istate,Jstate)
* Then set it to zero because we are within the diagonal approximation
                H0(Istate,Jstate) = 0.0d0
              else
* Just print out the already computed diagonal element
                write(6,'(A3,I4,A3,F16.8)')
     &                  ' < ',MSTATE(Istate),' | ', H0(Istate,Jstate)
              end if
            end do
            write(6,*)
          end if
        end if

* End of long loop over Jstate
      end do
* End timer Fock matrix build
      call timing(CPU1,CPU,TIO1,TIO)
      CPUFMB=CPU1-CPU0
      TIOFMB=TIO1-TIO0

* In case of a XMS calculation, i.e. Ngrp > 1 and not DW, transform
* the CI arrays of this group of states to make the Fock matrix
* diagonal in the model space.
* Note that this is only done for XMS here. For XDW and RMS it is
* done in the xdwinit subroutine. This code duplication is silly,
* but we will get rid of it once we drop the groups of states
      if (Ngrp.gt.1.and.IFXMS.and.(.not.IFDW)) then

* In case of XMS-CASPT2, printout H0 in original basis
        if (IPRGLB.ge.USUAL) then
          write(u6,*)
          write(u6,*)' H0 in the original model space basis:'
          call prettyprint(H0,Ngrp,Ngrp)
        end if
* Diagonalize H0 and save eigenvectors in U0
        call eigen(H0,U0,Ngrp)

* Transform the Fock matrix in the new basis
        call transmat(H0,U0,Ngrp)
        if (IPRGLB.ge.USUAL) then
          write(u6,*)' H0 eigenvectors:'
          call prettyprint(U0,Ngrp,Ngrp)
        end if
        if (IPRGLB.ge.DEBUG) then
          write(u6,*)' H0 in the rotated model space basis:'
          call prettyprint(H0,Ngrp,Ngrp)
        end if

* As well as Heff
        call transmat(Heff,U0,Ngrp)
        if (IPRGLB.ge.VERBOSE) then
          write(u6,*)' Heff[1] in the rotated model space basis:'
          call prettyprint(Heff,Ngrp,Ngrp)
        end if

* Mix the CI arrays according to the H0 eigenvectors. Assume we can
* put all the original ones in memory, but put the resulting vectors
* one by one in a buffer.
        if (IPRGLB.ge.VERBOSE) then
          write(u6,'(A)')' The CASSCF states are now rotated'//
     &                  ' according to the H0 eigenvectors'
          write(u6,*)
        end if

        call mma_allocate(CIref,Nconf,Ngrp,Label='CIRef')
* Load the CI arrays into memory
        do I=1,Ngrp
          call loadCI(CIref(:,I),I)
        end do

        call mma_allocate(CIXMS,Nconf,Label='CIXMS')
        do J=1,Ngrp
* Transform the states
          call dgemm_('N','N',Nconf,1,Ngrp,
     &               One,CIREF,Nconf,U0(:,J),Ngrp,
     &               Zero,CIXMS,Nconf)

* Write the rotated CI coefficients back into LUCIEX and REPLACE the
* original unrotated CASSCF states. Note that the original states
* are still available in the JobIph file
          call writeCI(CIXMS,J)

          if (IPRGLB.ge.VERBOSE) then
            write(u6,'(1x,a,i3)')
     &      ' The CI coefficients of rotated model state nr. ',MSTATE(J)
            call PRWF_CP2(STSYM,NCONF,CIXMS,CITHR)
          end if
        end do

        call mma_deallocate(CIREF)
        call mma_deallocate(CIXMS)

      end if

      Call  GPRINI_FINISH()

      Contains

      Subroutine  GPRINI_FINISH()

      Call mma_deallocate(HONE)
* We now know FIFA as expressed in initial RAS (natural) orbitals.
* Transform it to a new basis in which the non-diagonal couplings
* between subspaces (inactive, ras1, etc) are zero. As a by-product,
* the CI arrays will be transformed so they still represent the
* model functions, but using the new orbitals.
* Note that the matrices FIFA, FIMO, etc are transformed as well

      CALL ORBCTL(CMO,NCMO,TORB,SIZE(TORB),FIFA,SIZE(FIFA),
     &            FIMO,SIZE(FIMO))

* In subroutine stini, the individual RHS, etc, arrays will be computed
* for the states. If this is a true XMS calculation (Ngrp > 1) then
* there is one data set that is in common for these calculations,
* namely the transformed MO integrals (if conventional), or the
* transformed Cholesky vectors (if IfChol), so these are computed here

      CALL TIMING(CPU0,CPU,TIO0,TIO)

* TRACHO3 computes MO-transformed Cholesky vectors without computing
* Fock matrices
* TRACTL(nCMO,CMO,0) computes transformed 2-body MO integrals
      if (IfChol) then
          call TRACHO3(CMO,NCMO)
      else
          if (.not. DoFCIQMC) call TRACTL(nCMO,CMO,0)
      end if

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUINT=CPU1-CPU0
      TIOINT=TIO1-TIO0
      CMOPT2(:)=CMO(:)

      call mma_deallocate(CMO_Internal)
      nullify(CMO)

#ifdef _DMRG_
      if (DMRG) then
************************************************************************
* load back two-electron integrals (pu|vx)
************************************************************************
        ! Load the integrals in memory
        call read_integrals()

        ! set to compute 2-, 3- and 4-rdm
        call qcmaquis_interface_set_param('MEASURE[1rdm]','1')
        call qcmaquis_interface_set_param('MEASURE[2rdm]','1')
        call qcmaquis_interface_set_param('MEASURE[3rdm]','1')
        ! call qcmaquis_interface_set_param('MEASURE[4rdm]','1')
      end if
#endif
* ---------------------------------------------------------------------
      CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
      CPUGIN=CPTF10-CPTF0
      TIOGIN=TIOTF10-TIOTF0
* ---------------------------------------------------------------------
      End Subroutine  GPRINI_FINISH
      end SUBROUTINE GRPINI
