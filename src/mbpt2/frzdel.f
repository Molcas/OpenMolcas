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
      SUBROUTINE FrzDel(NREMO,IREMO,EOCC,E1,
     *                  NREME,IREME,EEXT,E2,
     *                  CMO,CMO1,ISEQ)
C
C     MOLPT2(MOLCAS) SUBROUTINE
C     CODED AJS, MAR. 15, 1990
C
C     THIS SUBROUTINE IS USED TO MOVE THE ADDITIONALLY FROZEN AND
C     ADDITIONALLY DELETED ORBITALS TO THE BOTTOM OR TO THE TOP
C     RESPECTIVELY, OF THE EIGENVECTOR LIST. THE ORBITAL ENERGIES='PRIN'
C     ARE REARRANGED ACCORDINGLY
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "files_mbpt2.fh"
C
      DIMENSION NREMO(*),IREMO(8,*),EOCC(*),E1(*)
      DIMENSION NREME(*),IREME(8,*),EEXT(*),E2(*)
      DIMENSION CMO(*),CMO1(*),ISEQ(*)
C
      Call qEnter('FrzDel')
C
      IAD0=1
      IADR=1
      IEO=0
      IENO=0
      IENE=0
      IEE=0
      DO 100 ISYM=1,NSYM
        IADF=IADR
        IADO=IADF+NBAS(ISYM)*(NFRO(ISYM)+NREMO(ISYM))
        IADE=IADF+NBAS(ISYM)*(NFRO(ISYM)+NOCC(ISYM))
        IADD=IADF+NBAS(ISYM)*(NBAS(ISYM)-NDEL(ISYM)-NREME(ISYM))
        DO 1 I=1,NBAS(ISYM)
          ISEQ(I)=I
1       CONTINUE
        DO 2 I=1,NFRO(ISYM)
          ISEQ(I)=0
2       CONTINUE
        DO 3 I=NBAS(ISYM),NBAS(ISYM)-NDEL(ISYM)+1,-1
          ISEQ(I)=0
3       CONTINUE
        DO 4 I=1,NREMO(ISYM)
           J=IREMO(ISYM,I)
           ISEQ(J)=0
4       CONTINUE
        DO 5 I=1,NREME(ISYM)
           J=IREME(ISYM,I)+NFRO(ISYM)+NOCC(ISYM)
           ISEQ(J)=0
5       CONTINUE
C
        DO 10 I=1,NFRO(ISYM)+NOCC(ISYM)
          IAD=IAD0+(I-1)*NBAS(ISYM)
          IF (ISEQ(I).EQ.0) THEN
*           frozen, move to appropriate place in symmetry block
*           observe the sequence: frozen/occupied/external/deleted
            CALL DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADF),1)
            IADF=IADF+NBAS(ISYM)
          ELSE
*           occupied...
            CALL DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADO),1)
            IADO=IADO+NBAS(ISYM)
            IENO=IENO+1
            EOCC(IENO)=E1(IEO+I-NFRO(ISYM))
          ENDIF
10      CONTINUE
        DO 20 I=NFRO(ISYM)+NOCC(ISYM)+1,NBAS(ISYM)
          IAD=IAD0+(I-1)*NBAS(ISYM)
          IF (ISEQ(I).EQ.0) THEN
*           delete, move to appropriate place in symmetry block
*           observe the sequence: frozen/occupied/external/deleted
            CALL DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADD),1)
            IADD=IADD+NBAS(ISYM)
          ELSE
*           external...
            CALL DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADE),1)
            IADE=IADE+NBAS(ISYM)
            IENE=IENE+1
            EEXT(IENE)=E2(IEE+I-NFRO(ISYM)-NOCC(ISYM))
          ENDIF
20      CONTINUE
        IADR=IADR+NBAS(ISYM)**2
        IAD0=IADR
        IEO=IEO+NOCC(ISYM)
        IEE=IEE+NEXT(ISYM)
C
C       UPDATE NFRO, NOCC, NEXT, NDEL, AND NORB
C
        NFRO(ISYM)=NFRO(ISYM)+NREMO(ISYM)
        NOCC(ISYM)=NOCC(ISYM)-NREMO(ISYM)
        NDEL(ISYM)=NDEL(ISYM)+NREME(ISYM)
        NEXT(ISYM)=NEXT(ISYM)-NREME(ISYM)
        NORB(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
100   CONTINUE
C
      Call qExit('FrzDel')
      RETURN
      END
