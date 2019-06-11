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
      SUBROUTINE MKTDAO(CMO,TDMO,TDAO,SCR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO(NCMO),TDMO(NBAST,NBAST),TDAO(NBAST,NBAST)
      DIMENSION SCR(NBMAX,NBMAX)

#include "SysDef.fh"

#include "mrci.fh"
C REORDER TDMO (USE TDAO AS TEMPORARY STORAGE):
      CALL FZERO(TDAO,NBAST**2)
      DO 10 I=1,NORBT
        II=ICH(I)
        IF(II.LE.0) GOTO 10
        DO 5 J=1,NORBT
          JJ=ICH(J)
          IF(JJ.LE.0) GOTO 5
          TDAO(I,J)=TDMO(II,JJ)
5       CONTINUE
10    CONTINUE
      CALL DCOPY_(NBAST**2,TDAO,1,TDMO,1)
      CALL FZERO(TDAO,NBAST**2)
      IECMO=0
      IEO=0
      IEB=0
      DO 100 ISYM=1,NSYM
        ISO=IEO+1
        ISB=IEB+1
        NO=NORB(ISYM)
        NB=NBAS(ISYM)
        IEO=IEO+NO
        IEB=IEB+NB
C ORBITALS PRE-FROZEN IN MOTRA, OR FROZEN IN MRCI:
        NF=NFMO(ISYM)+NFRO(ISYM)
        NBF=NB*NF
        ISCMO=IECMO+1
        IECMO=IECMO+NBF
C ORBITALS EXPLICITLY USED IN CI:
        NCO=NO-NFRO(ISYM)-NDEL(ISYM)
        ISCO=ISO+NFRO(ISYM)
        NBCO=NB*NCO
        ISCMO=IECMO+1
        IECMO=IECMO+NBCO
        IF(NCO.GT.0) THEN
          CALL DGEMM_('N','N',
     &                NB,NCO,NCO,
     &                1.0d0,CMO(ISCMO),NB,
     &                TDMO(ISCO,ISCO),NBAST,
     &                0.0d0,SCR,NB)
          CALL DGEMM_('N','T',NB,NB,NCO,1.0D00,SCR,NB,
     *               CMO(ISCMO),NB,1.0D00,TDAO(ISB,ISB),NBAST)
        END IF
C ORBITALS PRE-DELETED IN MOTRA OR DELETED IN MRCI:
        ND=NDMO(ISYM)+NDEL(ISYM)
        NBD=NB*ND
        ISCMO=IECMO+1
        IECMO=IECMO+NBD
100   CONTINUE
      RETURN
      END
