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
      SUBROUTINE GRPINI(IGROUP,NGRP,JSTATE_OFF,HEFF,H0,U0)
      IMPLICIT REAL*8 (A-H,O-Z)
* 2012  PER-AKE MALMQVIST
* Multi-State and XMS initialization phase
* Purpose: For a selected set IGROUP, create a set of CMO coefficients
* and orbital energies which are to be used in common for the
* calculations belonging to this set. Also, change the CI arrays in this
* group such that they diagonalize the H0 matrix.
* The states in the group can be obtained from the ordered MSTATE array,
* for which a group offset JSTATE_OFF is passed in.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "warnings.fh"
#include "stdalloc.fh"
      LOGICAL IF_TRNSF
      CHARACTER(27)  STLNE2
      real(8) Heff(Nstate,Nstate)
      real(8) H0(Nstate,Nstate)
      real(8) U0(Nstate,Nstate)

      CALL QENTER('GRPINI')
* ---------------------------------------------------------------------
* Number of states in this group.
      IF (IPRGLB.EQ.DEBUG) THEN
        write(6,*)' Entered GRPINI.'
        write(6,*)' NSTATE=',NSTATE
        write(6,*)' The MSTATE array:'
        write(6,'(1x,20I4)')(MSTATE(J),J=1,NSTATE)
        write(6,*)' IGROUP,NGRP=',IGROUP,NGRP
      END IF

      IF (NGRP.EQ.0) THEN
        WRITE(6,*) ' Number of states in the (X)MS group is 0!'
        WRITE(6,*) ' This should never happen, aborting...'
        CALL ABEND
      END IF

      Write(STLNE2,'(A,I3)')'Initial phase for group ',IGROUP
      Call StatusLine('CASPT2:',STLNE2)
      IF(IPRGLB.GE.USUAL) THEN
       If(.not.IFNOPT2) Then
        WRITE(6,'(20A4)')('****',I=1,20)
        WRITE(6,'(A,I3)')
     &  ' Multi-State initialization phase begins for group ',IGROUP
        WRITE(6,'(20A4)')('----',I=1,20)
        CALL XFlush(6)
       End If
      END IF
* ---------------------------------------------------------------------

* Load CASSCF MO coefficients
      call getmem('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      call ddafile(LUONEM,2,WORK(LCMO),NCMO,IDISK)
      IAD1M(2)=IDISK
      call ddafile(LUONEM,1,WORK(LCMO),NCMO,IDISK)
      IEOF1M=IDISK

* Loop over states, selecting those belonging to this group.
* For each such state, compute the Fock matrix in original MO basis,
* and then the zeroth-order Hamiltonian elements between states.

* Timer for Fock matrix build
      call timing(CPU0,CPU,TIO0,TIO)
* Loop over all states in group
      do J=1,Ngrp
        Jstate=J+JSTATE_OFF

* Copy the 1-RDM of Jstate from LDMIX into LDREF
        CALL DCOPY_(NDREF,WORK(LDMIX+(Jstate-1)*NDREF),1,WORK(LDREF),1)

* Compute the Fock matrix in MO basis for state Jstate
* INTCTL1/INTCTL2 call TRACTL(0) and other routines to compute the
* Fock matrix in MO basis: FIMO, FAMO, FIFA and orbital energies
        if (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis
          IF_TRNSF=.FALSE.
          call INTCTL2(IF_TRNSF)
        else
* INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis
          call INTCTL1(WORK(LCMO))
          call dcopy_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
        end If

c Modify the Fock matrix if needed
        IF (FOCKTYPE.NE.'STANDARD') THEN
           CALL NEWFOCK(WORK(LFIFA))
        END IF

* NN.15, TODO:
* MKFOP and following transformation are skipped in DMRG-CASPT2 run
* for the time, this will be fixed later to implement DMRG-MS-CASPT2
        IF (DoCumulant) GoTo 100

* Loop over bra functions
        do I=1,Ngrp
          Istate=I+JSTATE_OFF
* Compute matrix element and put it into H0
          call FOPAB(WORK(LFIFA),Istate,Jstate,H0(Istate,Jstate))
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
            write(6,'(10X,6X,A3,I4,A3)') ' | ', Jstate, ' > '
            do Istate=1,Nstate
              if (Istate.ne.Jstate) then
* Compute matrix element and print it out
                call FOPAB(WORK(LFIFA),Istate,Jstate,H0(Istate,Jstate))
                write(6,'(A3,I4,A3,F16.8)')
     &                  ' < ',Istate,' | ', H0(Istate,Jstate)
* Then set it to zero becuase we are within the diagonal approximation
                H0(Istate,Jstate) = 0.0d0
              else
* Just print out the already computed diagonal element
                write(6,'(A3,I4,A3,F16.8)')
     &                  ' < ',Istate,' | ', H0(Istate,Jstate)
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
* diagonal in the model space
      if (Ngrp.gt.1.and.IFXMS.and.(.not.IFDW)) then

* In case of XMS-CASPT2, printout H0 in original basis
        if (IPRGLB.ge.VERBOSE) then
          write(6,*)
          write(6,*)' H0 in the original model space basis:'
          call prettyprint(H0,Ngrp,Ngrp)
        end if
* Diagonalize H0 and save eigenvectors in U0
        call eigen(H0,U0,Ngrp)

* Transform the Fock matrix in the new basis
        call transmat(H0,U0,Ngrp)
        if (IPRGLB.ge.VERBOSE) then
          write(6,*)' H0 eigenvectors:'
          call prettyprint(U0,Ngrp,Ngrp)
        end if
        if (IPRGLB.ge.DEBUG) then
          write(6,*)' H0 in the rotated model space basis:'
          call prettyprint(H0,Ngrp,Ngrp)
        end if

* As well as Heff
        call transmat(Heff,U0,Ngrp)
        if (IPRGLB.ge.DEBUG) then
          write(6,*)' Heff[1] in the rotated model space basis:'
          call prettyprint(Heff,Ngrp,Ngrp)
        end if

       if(IFXMS) then
        call prrotmat(NGRP,U0,HEFF,NSTATE,IFSILPrRot)
       end if

* Mix the CI arrays according to the H0 eigenvectors. Assume we can
* put all the original ones in memory, but put the resulting vectors
* one by one in a buffer.
        if (IPRGLB.ge.VERBOSE) then
          write(6,'(A)')' The CASSCF states are now rotated'//
     &                  ' according to the H0 eigenvectors'
          write(6,*)
        end if

        call getmem('CIREF','ALLO','REAL',LCIref,Ngrp*Nconf)
* Load the CI arrays into memory
        do I=1,Ngrp
          call loadCI(WORK(LCIref+Nconf*(I-1)),I)
        end do

        call getmem('CIXMS','ALLO','REAL',LCIXMS,Nconf)
        do J=1,Ngrp
* Transform the states
          call dgemm_('N','N',Nconf,1,Ngrp,
     &               1.0D0,WORK(LCIREF),Nconf,U0(:,J),Ngrp,
     &               0.0D0,WORK(LCIXMS),Nconf)

* Write the rotated CI coefficients back into LUCIEX and REPLACE the
* original unrotated CASSCF states. Note that the original states
* are still available in the JobIph file
          call writeCI(WORK(LCIXMS),J)

          if (IPRGLB.ge.VERBOSE) then
            write(6,'(1x,a,i3)')
     &      ' The CI coefficients of rotated model state nr. ',MSTATE(J)
            call PRWF_CP2(LSYM,NCONF,WORK(LCIXMS),CITHR)
          end if
        end do

        call getmem('CIREF','FREE','REAL',LCIREF,Ngrp*Nconf)
        call getmem('CIXMS','FREE','REAL',LCIXMS,Nconf)

      end if

 100  continue

* We now know FIFA as expressed in initial RAS (natural) orbitals.
* Transform it to a new basis in which the non-diagonal couplings
* between subspaces (inactive, ras1, etc) are zero. As a by-product,
* the CI arrays will be transformed so they still represent the
* model functions, but using the new orbitals.
* Note that the matrices FIFA, FIMO, etc are transformed as well

      CALL ORBCTL(WORK(LCMO))

* In subroutine stini, the individual RHS, etc, arrays will be computed
* for the states. If this is a true XMS calculation (Ngrp > 1) then
* there is one data set that is in common for these calculations,
* namely the transformed MO integrals (if conventional), or the
* transformed Cholesky vectors (if IfChol), so these are computed here

      CALL TIMING(CPU0,CPU,TIO0,TIO)
      if (IfChol) then
* TRACHO3 computes MO-transformed Cholesky vectors without computing
* Fock matrices
        call TRACHO3(WORK(LCMO))
      else
* TRACTL(0) computes transformed 2-body MO integrals
        call TRACTL(0)
      end if
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUINT=CPU1-CPU0
      TIOINT=TIO1-TIO0
      call dcopy_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)

      call getmem('LCMO','FREE','REAL',LCMO,NCMO)

      CALL QEXIT('GRPINI')
      return
      end
