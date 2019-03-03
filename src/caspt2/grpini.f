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
************************************************************************
      SUBROUTINE GRPINI(IGROUP,NGRP,JSTATE_OFF,HEFF)
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
      REAL*8 NORMFAC

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
* For each such state, compute the one-electron Hamilonian to be used
* in the CASPT2 H0, in original MO basis, and finally replace it with
* the average over the group.
* Note that, in principle, also FAMO and DREF should be averaged over
* the states, but since we never use them during the XMS initialization
* we don't compute them.

      NFIFA_AVE=NFIFA
      CALL GETMEM('FIFA_AVE','ALLO','REAL',LFIFA_AVE,NFIFA_AVE)
      CALL DCOPY_(NFIFA_AVE,0.0D0,0,WORK(LFIFA_AVE),1)
      SCL=1.0D0/DBLE(NGRP)

      CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)

      IF (IFDW) THEN
        NGRP = NSTATE
        WRITE(6,*)'    ZETA = ',NZETA
      END IF

      DO ISTATE=1,NGRP
        JSTATE=JSTATE_OFF+ISTATE
        IF (IFDW) THEN
          JSTATE=ISTATE
        END IF

* If it is a dw-caspt2 calculation, compute the relative weight
        IF (IFDW) THEN
          WRITE(6,*)
          WRITE(6,*)'    ALPHA = ',IGROUP
          EALPHA = REFENE(IGROUP)
          WRITE(6,*)'   EALPHA = ',EALPHA
          WRITE(6,*)'  ----------------------------'
          WRITE(6,*)'     BETA = ',JSTATE
          EBETA = REFENE(JSTATE)
          WRITE(6,*)'    EBETA = ',EBETA
          WRITE(6,*)'  ----------------------------'

* compute normalization factor
          NORMFAC = 0.0D0
          DO I=1,NSTATE
            EGAMMA  = REFENE(I)
            NORMFAC = NORMFAC + EXP(-NZETA*(EALPHA - EGAMMA)**2)
          END DO
          WRITE(6,*)'  NORMFAC = ',NORMFAC
          WRITE(6,*)'  ----------------------------'

* compute the weight Wab
          SCL = EXP(-NZETA*(EALPHA - EBETA)**2)/NORMFAC
          WRITE(6,*)'      SCL = ',SCL
        END IF


* Accumulate the average active density matrix over this group.
        IF(ISCF.NE.0) THEN
* Then we still need the "CI array": It is used in subroutine calls
         WORK(LCI)=1.0D0
        ELSE IF(DoCumulant) THEN
*          write(6,*) 'Cumulant approximated 4RDM'
         WORK(LCI)=0.0D0
        ELSE
* Get the CI array:
         ID=IDCIEX
* This loop is just to move ID to the right place in the file
* so we can read the CI coeffs for the correct state.
* Basically, we want to read from JSTATE
         DO I=1,JSTATE-1
           CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,ID)
         END DO
         CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,ID)
        END IF
* We may want to write out the CI array.
        IF(IPRGLB.GE.VERBOSE .AND. ORBIN.EQ.'NO TRANS') THEN
          WRITE(6,*)
          WRITE(6,*)' CI array of CASSCF state nr ',MSTATE(JSTATE)
          CALL PRWF_CP2(LSYM,NCONF,WORK(LCI),CITHR)
        END IF

* POLY2: Computing 1- and 2-particle active density matrices GAMMA1 and GAMMA2
        CALL POLY2(WORK(LCI))
* GETDPREF: Restructure GAMMA1 and GAMMA2, as DREF and PREF arrays.
        CALL GETDPREF(WORK(LDREF),WORK(LPREF))

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
* Compute average Fock matrix:
        CALL DAXPY_(NFIFA_AVE,SCL,WORK(LFIFA),1,WORK(LFIFA_AVE),1)

      END DO
      CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)

* Replace FIFA with average Fock matrix:
      CALL DCOPY_(NFIFA,WORK(LFIFA_AVE),1,WORK(LFIFA),1)
      CALL GETMEM('FIFA_AVE','FREE','REAL',LFIFA_AVE,NFIFA_AVE)

* NN.15
* TODO : MKFOP and following transformation are skipped in DMRG-CASPT2 run
*        for the time, this will be fixed later to implement DMRG-MS-CASPT2.
      IF(DoCumulant) GoTo 100

* Resetting NGRP to 1
      IF (IFDW) THEN
        NGRP = 1
      END IF

* Compute elements of Hamiltonian matrix obtained as
* <BRA|FOP|KET> where FOP is the average Fock operator (FIFA)

      CALL GETMEM('FOPXMS','ALLO','REAL',LFOPXMS,NGRP**2)
      CALL DCOPY_(NGRP**2,0.0D0,0,WORK(LFOPXMS),1)

      CALL MKFOP(WORK(LFIFA),NGRP,JSTATE_OFF,WORK(LFOPXMS))

      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(6,*)' GRPINI computed FOPXMS:'
       DO I=1,NGRP
        WRITE(6,'(1x,5F16.8)')(WORK(LFOPXMS-1+I+NGRP*(J-1)),J=1,NGRP)
       END DO
      END IF

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
* Also change HEFF accordingly. Originally, it was diagonal:
       CALL GETMEM('HTMP1','ALLO','REAL',LHTMP1,NGRP**2)
       CALL GETMEM('HTMP2','ALLO','REAL',LHTMP2,NGRP**2)
       DO J1=1,NGRP
        IK1=JSTATE_OFF+J1
        DO J2=1,NGRP
         IK2=JSTATE_OFF+J2
         WORK(LHTMP1-1+J1+NGRP*(J2-1))=HEFF(IK1,IK2)
        END DO
       END DO
       CALL DGEMM_('T','N',NGRP,NGRP,NGRP,
     &              1.0d0,WORK(LEVEC),NGRP,WORK(LHTMP1),NGRP,
     &              0.0d0,WORK(LHTMP2),NGRP)
       CALL DGEMM_('N','N',NGRP,NGRP,NGRP,
     &              1.0d0,WORK(LHTMP2),NGRP,WORK(LEVEC),NGRP,
     &              0.0d0,WORK(LHTMP1),NGRP)
       DO J1=1,NGRP
        IK1=JSTATE_OFF+J1
        DO J2=1,NGRP
         IK2=JSTATE_OFF+J2
         HEFF(IK1,IK2)=WORK(LHTMP1-1+J1+NGRP*(J2-1))
        END DO
       END DO
       CALL GETMEM('HTMP1','FREE','REAL',LHTMP1,NGRP**2)
       CALL GETMEM('HTMP2','FREE','REAL',LHTMP2,NGRP**2)

       IF(IPRGLB.GE.DEBUG) THEN
        WRITE(6,*) 'HEFF AFTER TRANSFORMATION IN THE NEW "XMS" BASIS:'
        DO J1=1,NSTATE
         WRITE(6,'(5F16.8)')(HEFF(J1,J2),J2=1,NSTATE)
        END DO
       END IF
* and then, transform the CI arrays. Assume we can put all the
* original ones in memory, but put the resulting vectors one by
* one in a buffer.
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
     &             1.0D0,WORK(LCIREF),NCONF,
     &             WORK(LEVEC+NGRP*(J-1)),NGRP,
     &             0.0D0,WORK(LCIXMS),NCONF)
        IK=JSTATE_OFF+J
        ID=IDCIEX
        DO I=1,IK-1
         CALL DDAFILE(LUCIEX,0,WORK(LCIXMS),NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,1,WORK(LCIXMS),NCONF,ID)
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
      CALL ORBCTL(WORK(LCMO))

* In subroutine stini, the individual RHS, etc, arrays will be computed for
* the states. If this is a true XMS calculation (NGRP.gt.1) then there is one
* data set that is in common for these calculations, namely the transformed
* MO integrals (if conventional), or the transformed Cholesky vectors (if
* IfChol), so these are computed here:
      If (IfChol) then
* TRACHO3 computes MO-transformed Cholesky vectors without computing Fock matrices.
        CALL TRACHO3(WORK(LCMO))
      Else
* TRACTL(0) computes transformed 2-body MO integrals.
        Call TRACTL(0)
      End If
      CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)

      CALL QEXIT('GRPINI')
      RETURN
      END
