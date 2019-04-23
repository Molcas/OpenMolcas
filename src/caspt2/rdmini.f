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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE RDMINI

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

      CALL QENTER('RDMINI')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)' Entered RDMINI.'
      END IF

* Get CASSCF MO coefficients
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WORK(LCMO),NCMO,IDISK)

      CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)

* Initialize array of weights with all zeros
      CALL GETMEM('LWGT','ALLO','REAL',LWGT,NSTATE**2)
      CALL DCOPY_(NSTATE**2,0.0D0,0,WORK(LWGT),1)

* Initialize array of 1-RDMs with zeros
      CALL DCOPY_(NSTATE*NDREF,0.0D0,0,WORK(LDMIX),1)
* Initialize array for temporary storage of DMIX
      ! CALL GETMEM('LDTMP','ALLO','REAL',LDTMP,NDREF*NSTATE)

      DO ISTATE=1,NSTATE

        IF(ISCF.NE.0) THEN
* Then we still need the "CI array": It is used in subroutine calls
          WORK(LCI)=1.0D0
        ELSE IF(DoCumulant) THEN
          WORK(LCI)=0.0D0
        ELSE
* Get the CI array
          ID=IDCIEX
          DO I=1,ISTATE-1
            CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,ID)
          END DO
          CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,ID)
        END IF

        IF(IPRGLB.GE.VERBOSE) THEN
          WRITE(6,*)
          WRITE(6,*)' CI array of CASSCF state nr. ',MSTATE(ISTATE)
          CALL PRWF_CP2(LSYM,NCONF,WORK(LCI),CITHR)
        END IF

* Compute the density matrices weights
        IF (IFDW) THEN
          EBETA = REFENE(ISTATE)
* Compute normalization factor
          DO JSTATE=1,NSTATE
            EALPHA = REFENE(JSTATE)
            FAC = 0.0D0
            DO KSTATE=1,NSTATE
              EGAMMA = REFENE(KSTATE)
              FAC = FAC + EXP(-NZETA*(EALPHA - EGAMMA)**2)
            END DO
            IJ = NSTATE*(ISTATE-1) + (JSTATE-1)
            WORK(LWGT+IJ) = EXP(-NZETA*(EALPHA - EBETA)**2)/FAC
          END DO
* If it is an XMS-CASPT2 calculation, all the weights are 1/NSTATE
        ELSE IF ( IFXMS .AND. (.NOT.IFDW) ) THEN
          II = NSTATE*(ISTATE-1)
          CALL DCOPY_(NSTATE,1.0D0/NSTATE,0,WORK(LWGT+II),1)
        ELSE
* It is a normal MS-CASPT2 and the weight vectors are the standard e_1, e_2, ...
          WORK(LWGT + (NSTATE*(ISTATE-1)) + (ISTATE-1)) = 1.0D0
        END IF

* Compute 1-particle active density matrix GAMMA1
        CALL POLY1(WORK(LCI))
* Restructure GAMMA1 as DREF array, but keep it in DMIX
        CALL GETDREF(WORK(LDREF))
        ! CALL GETDREF(WORK(LDMIX+NDREF*(ISTATE-1)))

        IFTEST = 1
        IF ( IFTEST.NE.0 ) THEN
          WRITE(6,*)' DREF for state nr. ',MSTATE(ISTATE)
          DO I=1,NASHT
            WRITE(6,'(1x,14f10.6)')(WORK(LDREF+(I*(I-1))/2+J-1),J=1,I)
          END DO
          WRITE(6,*)
        END IF

* Averaging the density
        DO JSTATE=1,NSTATE
          SCL = WORK(LWGT + (NSTATE*(ISTATE-1)) + (JSTATE-1))
          JOFF = NDREF*(JSTATE-1)
          CALL DAXPY_(NDREF,SCL,WORK(LDREF),1,WORK(LDMIX+JOFF),1)

          IFTEST = 1
          IF ( IFTEST.NE.0 ) THEN
            JJOFF = LDMIX+JOFF
            WRITE(6,*)' DMIX for state nr. ',MSTATE(JSTATE)
            DO I=1,NASHT
              WRITE(6,'(1x,14f10.6)')(WORK(JJOFF+(I*(I-1))/2+J-1),J=1,I)
            END DO
            WRITE(6,*)
          END IF
        END DO

      END DO

      IFTEST = 1
      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)
        DO I=1,NSTATE
          WRITE(6,'(1x,10f8.4)')(WORK(LWGT + (NSTATE*(I-1)) + (J-1)),
     &    J=1,NSTATE)
        END DO
      END IF

      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
      CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)
      CALL GETMEM('LWGT','FREE','REAL',LWGT,NSTATE**2)

      CALL QEXIT('RDMINI')
      RETURN
      END