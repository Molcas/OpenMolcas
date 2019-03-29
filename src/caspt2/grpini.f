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
      SUBROUTINE GRPINI(IGROUP,NGRP,JSTATE_OFF,HEFF,H0)
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
      LOGICAL IF_TRNSF
      CHARACTER(27)  STLNE2
      REAL*8 HEFF(NSTATE,NSTATE)
      REAL*8 H0(NSTATE,NSTATE)

      CALL QENTER('GRPINI')
* ---------------------------------------------------------------------
* Number of states in this group.
      IF (IPRGLB.EQ.DEBUG) THEN
        write(6,*)' Entered GRPINI.'
        write(6,*)' NSTATE=',NSTATE
        write(6,*)' The MSTATE array:'
        write(6,'(1x,20I4)')( MSTATE(J),J=1,NSTATE)
        write(6,*)' IGROUP,NGRP=',IGROUP,NGRP
      END IF

      IF (NGRP.EQ.0) THEN
        WRITE(6,*) ' Number of states in the (X)MS group is 0!'
        WRITE(6,*) ' This should never happen, aborting...'
        CALL ABEND
      END IF

* ---------------------------------------------------------------------
      Write(STLNE2,'(A,I3)')'Initial phase for group ',IGROUP
      Call StatusLine('CASPT2:',STLNE2)
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,'(20A4)')('****',I=1,20)
        WRITE(6,'(A,I3)')
     &  ' Multi-State initialization phase begins for group ',IGROUP
        WRITE(6,'(20A4)')('----',I=1,20)
        CALL XFlush(6)
      END IF

* ---------------------------------------------------------------------
* GET ORIGINAL CASSCF CMO COEFFICIENTS.
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WORK(LCMO),NCMO,IDISK)
* Also (for temporary back-compatibility with older code) save as
*  'current' CMO data on LUONEM:
      IAD1M(2)=IDISK
      CALL DDAFILE(LUONEM,1,WORK(LCMO),NCMO,IDISK)
      IEOF1M=IDISK

* ---------------------------------------------------------------------
* Loop over states, selecting those belonging to this group.
* For each such state, compute the Fock matrix in original MO basis,
* and then the zeroth-order Hamiltonian elements between states.

      CALL GETMEM('FOPXMS','ALLO','REAL',LFOPXMS,NGRP**2)
      CALL DCOPY_(NGRP**2,0.0D0,0,WORK(LFOPXMS),1)

* Loop over all states in group
      DO J2=1,NGRP
        JSTATE=J2+JSTATE_OFF

* Copy the 1-RDM of JSTATE in LDREF
        CALL DCOPY_(NDREF,WORK(LDMIX+(JSTATE-1)*NDREF),1,WORK(LDREF),1)

* INTCTL1/INTCTL2 call TRACTL(0), and other routines, for FIMO, FAMO,
* FIFA and orbital energies.
        If (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
          IF_TRNSF=.FALSE.
          CALL INTCTL2(IF_TRNSF)
        Else
* INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis.
          CALL INTCTL1(WORK(LCMO))
          CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
        End If

c Modify the Fock matrix, if needed:
        IF(FOCKTYPE.NE.'STANDARD') THEN
           CALL NEWFOCK(WORK(LFIFA))
        END IF

* NN.15
* TODO : MKFOP and following transformation are skipped in DMRG-CASPT2 run
*        for the time, this will be fixed later to implement DMRG-MS-CASPT2.
        IF (DoCumulant) GoTo 100

* Compute elements of zeroth-order Hamiltonian matrix obtained as
* <BRA|FOP|KET> where FOP is the (average) Fock operator (FIFA)
        IKET=JSTATE
* Loop over bra functions:
        DO J1=1,NGRP
          IBRA=J1+JSTATE_OFF
* Compute matrix element and put it into FOPXMS:
          FOPEL = 0.0D0
          CALL FOPAB(WORK(LFIFA),IBRA,IKET,FOPEL)
          WORK(LFOPXMS + (J1-1) + (NGRP*(J2-1))) = FOPEL
        END DO

* Save f_pq(D) into shared memory LFIFA_ALL, required for
* DW-XMS later on in the loop over states to obtain quasi-canonical
* orbitals for each state independently
        WRITE(6,'(2x,A,I1,A)')'Saving F(D_',JSTATE,
     &                      ') into shared memory LFIFA_ALL...'
        CALL DCOPY_(NFIFA,WORK(LFIFA),1,WORK(LFIFA_ALL+
     &                                ((J2-1)*NFIFA)),1)

      END DO
* End of loop over states

      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' FOPXMS (Asymmetric):'
        DO I=1,NGRP
          WRITE(6,'(1x,5F16.8)')(WORK(LFOPXMS-1+I+NGRP*(J-1)),J=1,NGRP)
        END DO
      END IF

* Symmetrize FOPXMS (it really does something only for DW-XMS)
      DO I=1,NGRP
        DO J=I,NGRP
          FIJ = WORK(LFOPXMS+I-1+NGRP*(J-1))
          FJI = WORK(LFOPXMS+J-1+NGRP*(I-1))
          WORK(LFOPXMS+I-1+NGRP*(J-1)) = 0.5D0*(FIJ+FJI)
          WORK(LFOPXMS+J-1+NGRP*(I-1)) = 0.5D0*(FIJ+FJI)
        END DO
      END DO

      IF(IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' FOPXMS (Symmetric):'
        DO I=1,NGRP
          WRITE(6,'(1x,5F16.8)')(WORK(LFOPXMS+I-1+NGRP*(J-1)),J=1,NGRP)
        END DO
      END IF

* Store zeroth order energies
      DO I=1,NGRP
        DO J=1,NGRP
          H0(I+JSTATE_OFF,J+JSTATE_OFF) = WORK(LFOPXMS+I-1+NGRP*(J-1))
        END DO
      END DO

* Transform the CI arrays of this group of states, to make the FOP matrix diagonal.
* Note that the Fock matrix, etc are still assumed to be valid -- this seems
* illogical, but is the way XMS is defined -- else we would need to repeat the
* whole thing iteratively.
      IF(NGRP.GT.1.AND.IFXMS) THEN

        CALL GETMEM('EVEC','ALLO','REAL',LEVEC,NGRP**2)
        NSCR=(NGRP*(NGRP+1))/2
        CALL GETMEM('SCR','ALLO','REAL',LSCR,NSCR)
        CALL DIAFOP(NGRP,WORK(LFOPXMS),WORK(LSCR),WORK(LEVEC))
        CALL GETMEM('SCR','FREE','REAL',LSCR,NSCR)

* Transform H0 (= FOPXMS) and HEFF (which at this stage corresponds to the 1st-order
* corrected effective Hamiltonian) in the new basis that diagonalize FOPXMS
        CALL TRANSMAT(H0,WORK(LEVEC),NGRP)
        CALL TRANSMAT(HEFF,WORK(LEVEC),NGRP)

        IF(IPRGLB.GE.DEBUG) THEN
          WRITE(6,*) ' HEFF IN THE NEW "XMS" BASIS:'
          DO J1=1,NSTATE
            WRITE(6,'(1x,5F16.8)')(HEFF(J1,J2),J2=1,NSTATE)
          END DO
        END IF

        IF(IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,'(6X,A)')' Eigenvectors:'
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            DO J1=0,NSTATE-1
              WRITE(6,'(6x,5F16.8)')(WORK(LEVEC+J1+NGRP*(I-1)),
     &                               I=ISTA,IEND)
            END DO
            WRITE(6,*)
          END DO
        END IF

* And then, transform the CI arrays. Assume we can put all the
* original ones in memory, but put the resulting vectors one by
* one in a buffer.
        WRITE(6,*)
        WRITE(6,*)' Mixing the CASSCF states according to FOPXMS...'
        WRITE(6,*)
        CALL GETMEM('CIREF','ALLO','REAL',LCIREF,NGRP*NCONF)
        DO J=1,NGRP
          IK=JSTATE_OFF+J
          ID=IDCIEX
          DO I=1,IK-1
            CALL DDAFILE(LUCIEX,0,WORK(LCIREF),NCONF,ID)
          END DO
          CALL DDAFILE(LUCIEX,2,WORK(LCIREF+NCONF*(J-1)),NCONF,ID)
        END DO
        CALL GETMEM('CIXMS','ALLO','REAL',LCIXMS,NCONF)
        DO J=1,NGRP
          CALL DGEMM_('N','N',NCONF,1,NGRP,
     &                1.0D0,WORK(LCIREF),NCONF,
     &                WORK(LEVEC+NGRP*(J-1)),NGRP,
     &                0.0D0,WORK(LCIXMS),NCONF)
          IK=JSTATE_OFF+J
          ID=IDCIEX
          DO I=1,IK-1
            CALL DDAFILE(LUCIEX,0,WORK(LCIXMS),NCONF,ID)
          END DO
          CALL DDAFILE(LUCIEX,1,WORK(LCIXMS),NCONF,ID)

          IF(IPRGLB.GE.VERBOSE) THEN
            WRITE(6,'(1x,a,i3)')
     &      ' The CI coefficients for the MIXED state nr. ',MSTATE(IK)
            CALL PRWF_CP2(LSYM,NCONF,WORK(LCIXMS),CITHR)
          END IF

        END DO

        CALL GETMEM('CIREF','FREE','REAL',LCIREF,NGRP*NCONF)
        CALL GETMEM('CIXMS','FREE','REAL',LCIXMS,NCONF)
        CALL GETMEM('EVEC','FREE','REAL',LEVEC,NGRP**2)

      END IF

      CALL GETMEM('FOPXMS','FREE','REAL',LFOPXMS,NGRP**2)

 100  CONTINUE
* We now know FIFA, as expressed in initial RAS orbitals. Transform to use new
* orbitals, in which non-diagonal couplings within subspaces (inactive, ras1, etc)
* are zero. As a by-product, the CI arrays will be transformed so they still
* represent the XMS root functions, using the new orbitals.
* Also, the matrices FIFA, etc, are themselves transformed:

* This transformations are done here for all cases but DW-XMS.
* In that case, we need to obtain quasi-canonical orbitals for
* state independently, and we do that in the main caspt2 subroutine
* when looping over the states in the group.
      IF (IFDW.AND.IFXMS) THEN
        WRITE(6,*)
        WRITE(6,'(2x,A,A)')'Transformation to quasi-canonical MOs',
     &                    ' will be done in CASPT2 loop'
        WRITE(6,*)
        CALL DCOPY_(NFIMO,WORK(LFIMO),1,WORK(LFIMO_CMO),1)
        CALL DCOPY_(NHONE,WORK(LHONE),1,WORK(LHONE_CMO),1)
      ELSE
        WRITE(6,*)
        WRITE(6,'(2x,A,A)')'Transformation to quasi-canonical MOs',
     &                    ' will be done in GRPINI'
        WRITE(6,*)

        CALL ORBCTL(WORK(LCMO))

* In subroutine stini, the individual RHS, etc, arrays will be computed for
* the states. If this is a true XMS calculation (NGRP.gt.1) then there is one
* data set that is in common for these calculations, namely the transformed
* MO integrals (if conventional), or the transformed Cholesky vectors (if
* IfChol), so these are computed here:
        IF (IfChol) then
* TRACHO3 computes MO-transformed Cholesky vectors without computing Fock matrices.
          CALL TRACHO3(WORK(LCMO))
        ELSE
* TRACTL(0) computes transformed 2-body MO integrals.
          CALL TRACTL(0)
        END IF
        CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)

      END IF

      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)

      CALL QEXIT('GRPINI')
      RETURN
      END
