************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SIGMA(SGM,AREF,CI,INTSY,INDX,BMN,IBMN,BIAC2,
     &              BICA2,BFIN3,ISAB,AC1,AC2,BFIN4,ABIJ,
     &              AIBJ,AJBI,ASCR1,BSCR1,FSCR1,FSEC,FOCK,
*PAM04     &              BFIN5,ASCR2,BSCR2,FSCR2,DBK,CSPCK)
     &              BFIN5,ASCR2,BSCR2,FSCR2,DBK,ICSPCK)
      IMPLICIT REAL*8 (A-H,O-Z)
*PAM04      Integer INTSY(*),INDX(*),IBMN(*),ISAB(*)
      Integer INTSY(*),INDX(*),IBMN(*),ISAB(*),ICSPCK(*)
      Real*8 SGM(*),AREF(*),CI(*),BMN(*),BIAC2(*),BICA2(*)
      Real*8 BFIN3(*),AC1(*),AC2(*),BFIN4(*),ABIJ(*)
      Real*8 AIBJ(*),AJBI(*),ASCR1(*),BSCR1(*),FSCR1(*)
      Real*8 FSEC(*),FOCK(*),BFIN5(*),ASCR2(*),BSCR2(*),FSCR2(*)
*PAM04      Real*8 DBK(*),CSPCK(*)
      Real*8 DBK(*)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "mrci.fh"
      CALL DCOPY_(NCONF,0.0D0,0,SGM,1)

      CALL CSFTRA(' CSF',CI,AREF)
      SQGP=1.0D00
      SQG =1.0D00
      IF(ICPF.EQ.1) THEN
        SQGP=SQRT(GFAC)
        SQG=1.0D00/SQGP
        DO IREF=1,NREF
          ICSF=IREFX(IREF)
           CI(ICSF)=SQGP*CI(ICSF)
        END DO
      END IF

      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      CALL DIAGC(INTSY,CI,SGM)
      IF(IFIRST.EQ.0 .AND. ((IREST.EQ.1).OR.(ITER.GT.1)) ) THEN
        CALL ABCI_MRCI(INTSY,INDX,CI,SGM,BMN,IBMN,BIAC2,BICA2,BFIN3)
        CALL ABCD_MRCI(INTSY,INDX,ISAB,CI,SGM,AC1,AC2,BFIN4)
      END IF
*PAM04      CALL IJKL(INTSY,INDX,CI,SGM,FIJKL)
      CALL IJKL(INTSY,INDX,CI,SGM,WORK(LFIJKL))
c
      CALL GETMEM('BUF','ALLO','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','ALLO','INTE',LIBUF,NBITM3+2)
*PAM04      CALL FAIBJ(INTSY,INDX,CI,SGM,ABIJ,AIBJ,AJBI,BFIN1,BFIN1,
*PAM04     &           ASCR1,BSCR1,FSCR1,FSEC)
      CALL FAIBJ(INTSY,INDX,CI,SGM,ABIJ,AIBJ,AJBI,WORK(LBUF),
     &            IWORK(LIBUF),ASCR1,BSCR1,FSCR1,FSEC)
      CALL GETMEM('BUF','FREE','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','FREE','INTE',LIBUF,NBITM3+2)

      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IF(ITER.GT.0) THEN
        KTYP=1
* Switch KTYP=1 means AI is actually handling AIJK integrals.
      CALL GETMEM('BUF','ALLO','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','ALLO','INTE',LIBUF,NBITM3+2)
*PAM04        CALL AI(INTSY,INDX,CI,SGM,FOCK,BFIN5,BFIN5,ASCR2,BSCR2,
*PAM04     &          FSCR2,DBK,KTYP)
        CALL AI_MRCI(INTSY,INDX,CI,SGM,FOCK,WORK(LBUF),IWORK(LIBUF),
     &          ASCR2,BSCR2,FSCR2,DBK,KTYP)
      CALL GETMEM('BUF','FREE','REAL',LBUF,NBITM3)
      CALL GETMEM('IBUF','FREE','INTE',LIBUF,NBITM3+2)
      END IF
*PAM04      CALL FIJ(CSPCK,INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK)
      CALL FIJ_MRCI(ICSPCK,INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK)

      CALL DAXPY_(NCONF,POTNUC-ESHIFT,CI,1,SGM,1)
      IF(ICPF.EQ.1) THEN
        GINV=1.0D00/GFAC
        CALL DSCAL_(NCONF,GINV,SGM,1)
        DO IREF=1,NREF
          ICSF=IREFX(IREF)
           CI(ICSF)=SQG*CI(ICSF)
           SGM(ICSF)=SQGP*SGM(ICSF)
        END DO
      END IF

      CALL CSFTRA('MCSF',SGM,AREF)

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real_array(BFIN5)
      END
