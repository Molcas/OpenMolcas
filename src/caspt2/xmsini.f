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
      SUBROUTINE XMSINI(HEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "warnings.fh"
      REAL*8 HEFF(NSTATE,NSTATE)
      LOGICAL IF_TRNSF

      CALL QENTER('XMSINI')

      NGRP = NSTATE

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

      CALL GETMEM('FOPXMS','ALLO','REAL',LFOPXMS,NSTATE**2)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LFOPXMS),1)

* Loop over all states in group
      DO J2=1,NSTATE
        JSTATE=J2

* Copy the 1-RDM of JSTATE in LDREF
        CALL DCOPY_(NDREF,WORK(LDMIX+(JSTATE-1)*NDREF),1,WORK(LDREF),1)

* Compute sa-FIFA
        If (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
          IF_TRNSF=.FALSE.
          CALL INTCTL2(IF_TRNSF)
        Else
* INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis.
          CALL INTCTL1(WORK(LCMO))
          CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
        End If

* Compute elements of zeroth-order Hamiltonian matrix obtained as
* <BRA|FOP|KET> where FOP is the (average) Fock operator (FIFA)
        IKET=JSTATE
* Loop over bra functions:
        DO J1=1,NSTATE
          IBRA=J1
* Compute matrix element and put it into FOPXMS:
          FOPEL = 0.0D0
          CALL FOPAB(WORK(LFIFA),IBRA,IKET,FOPEL)
          WORK(LFOPXMS + (J1-1) + (NSTATE*(J2-1))) = FOPEL
        END DO

      END DO
* End of loop over states

      IF (.TRUE.) THEN
        NGRP = NSTATE
        WRITE(6,*)
        WRITE(6,*)' H0 (Asymmetric in XMSini):'
        DO I=1,NGRP
          WRITE(6,'(1x,5F16.8)')(WORK(LFOPXMS-1+I+NGRP*(J-1)),J=1,NGRP)
        END DO
      END IF

      IF (IPRGLB.GE.USUAL.AND.(NGRP.GT.1)) THEN
        WRITE(6,*)
        WRITE(6,*)' Zeroth-order Hamiltonian matrix (H0) in XMSini:'
        DO ISTA=1,NGRP,5
          IEND=MIN(ISTA+4,NGRP)
          WRITE(6,*)
          WRITE(6,'(1x,5I16)')(I,I=ISTA,IEND)
          DO I=ISTA,NGRP
            II0=(I*(I-1))/2
            WRITE(6,'(1x,I3,2X,5F16.8)')
     &            I,(WORK(LFOPXMS+I-1+NGRP*(J-1)),J=ISTA,IEND)
          END DO
        END DO
      END IF

* Transform the CI arrays of this group of states, to make the FOP matrix diagonal.
* Note that the Fock matrix, etc are still assumed to be valid -- this seems
* illogical, but is the way XMS is defined -- else we would need to repeat the
* whole thing iteratively.
      IF (.TRUE.) THEN

        CALL GETMEM('EVEC','ALLO','REAL',LEVEC,NGRP**2)
        CALL DIAFOP(NGRP,WORK(LFOPXMS),WORK(LEVEC))

* Transform H0 (= FOPXMS) and HEFF (which at this stage corresponds to the 1st-order
* corrected effective Hamiltonian) in the new basis that diagonalize FOPXMS
        ! CALL TRANSMAT(H0,WORK(LEVEC),NGRP)
        CALL TRANSMAT(HEFF,WORK(LEVEC),NGRP)

        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,'(6X,A)')' H0 eigenvectors:'
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            DO J1=0,NSTATE-1
              WRITE(6,'(6x,5F16.8)')(WORK(LEVEC+J1+NGRP*(I-1)),
     &                               I=ISTA,IEND)
            END DO
            WRITE(6,*)
          END DO
        END IF

        IF (IPRGLB.GE.VERBOSE) THEN
          WRITE(6,*)
          WRITE(6,*)' Heff[1] in H0 basis:'
          DO J1=1,NSTATE
            WRITE(6,'(1x,5F16.8)')(HEFF(J1,J2),J2=1,NSTATE)
          END DO
        END IF

* And then, transform the CI arrays. Assume we can put all the
* original ones in memory, but put the resulting vectors one by
* one in a buffer.
        WRITE(6,*)
        WRITE(6,*)' Mixing the CASSCF states according to XMS...'
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

      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)

      CALL QEXIT('XMSINI')
      RETURN
      END

