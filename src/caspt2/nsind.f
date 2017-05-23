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
      SUBROUTINE NSIND(INS,ISYM,ICASE,IP,IQ,IR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"

      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13) ICASE

   1  CONTINUE
      IIABS=INS+NIES(ISYM)
      IP=IINAIS(IIABS)
      IQ=0
      IR=0
      RETURN
   2  CONTINUE
  12  CONTINUE
        IIJABS=INS+NIGEJES(ISYM)
        IIABS=MIGEJ(1,IIJABS)
        IJABS=MIGEJ(2,IIJABS)
      GOTO 1213
   3  CONTINUE
  13  CONTINUE
        IIJABS=INS+NIGTJES(ISYM)
        IIABS=MIGTJ(1,IIJABS)
        IJABS=MIGTJ(2,IIJABS)
 1213 CONTINUE
      IP=IINAIS(IIABS)
      IQ=IINAIS(IJABS)
      IR=0
      RETURN
   4  CONTINUE
      IAABS=INS+NSES(ISYM)
CFUE  IP=IINAIS(IAABS)
      IP=IEXTIS(IAABS)
      IQ=0
      IR=0
      RETURN
   5  CONTINUE
      IIA=INS
      DO ISYMA=1,NSYM
        ISYMI=MUL(ISYMA,ISYM)
        NIA=NISH(ISYMI)*NSSH(ISYMA)
        IF(IIA.LE.NIA) THEN
          IA=1+(IIA-1)/NISH(ISYMI)
          II=IIA-NISH(ISYMI)*(IA-1)
          IAABS=IA+NSES(ISYMA)
          IIABS=II+NIES(ISYMI)
          IP=IINAIS(IIABS)
          IQ=IEXTIS(IAABS)
          IR=0
          RETURN
        END IF
        IIA=IIA-NIA
      END DO
      WRITE(6,*)'NSIND AIVX: Impossible situation.'
      CALL ABEND()
   6  CONTINUE
   7  CONTINUE
      IAIJ=INS
      NIJ = 0 ! dummy initialize
      DO ISYMA=1,NSYM
        ISYMIJ=MUL(ISYMA,ISYM)
        IF(ICASE.EQ.6) NIJ=NIGEJ(ISYMIJ)
        IF(ICASE.EQ.7) NIJ=NIGTJ(ISYMIJ)
        NA=NSSH(ISYMA)
        NAIJ=NA*NIJ
        IF(IAIJ.LE.NAIJ) THEN
          IIJ=1+(IAIJ-1)/NA
          IA=IAIJ-NA*(IIJ-1)
          IAABS=IA+NSES(ISYMA)
          IF(ICASE.EQ.6) THEN
            IIJABS=IIJ+NIGEJES(ISYMIJ)
            IIABS=MIGEJ(1,IIJABS)
            IJABS=MIGEJ(2,IIJABS)
          ELSE
            IIJABS=IIJ+NIGTJES(ISYMIJ)
            IIABS=MIGTJ(1,IIJABS)
            IJABS=MIGTJ(2,IIJABS)
          END IF
          IP=IEXTIS(IAABS)
          IQ=IINAIS(IIABS)
          IR=IINAIS(IJABS)
          RETURN
        END IF
        IAIJ=IAIJ-NAIJ
      END DO
      WRITE(6,*)'NSIND VJAI: Impossible situation.'
      CALL ABEND()
   8  CONTINUE
        IABABS=INS+NAGEBES(ISYM)
        IAABS=MAGEB(1,IABABS)
        IBABS=MAGEB(2,IABABS)
      GOTO 89
   9  CONTINUE
        IABABS=INS+NAGTBES(ISYM)
        IAABS=MAGTB(1,IABABS)
        IBABS=MAGTB(2,IABABS)
  89  CONTINUE
      IP=IEXTIS(IAABS)
      IQ=IEXTIS(IBABS)
      IR=0
      RETURN
  10  CONTINUE
  11  CONTINUE
      IIAB=INS
      NAB = 0 ! dummy initialize
      DO ISYMI=1,NSYM
        ISYMAB=MUL(ISYMI,ISYM)
        IF(ICASE.EQ.10) NAB=NAGEB(ISYMAB)
        IF(ICASE.EQ.11) NAB=NAGTB(ISYMAB)
        NI=NISH(ISYMI)
        NIAB=NI*NAB
        IF(IIAB.LE.NIAB) THEN
          IAB=1+(IIAB-1)/NI
          II=IIAB-NI*(IAB-1)
          IIABS=II+NIES(ISYMI)
          IF(ICASE.EQ.10) THEN
            IABABS=IAB+NAGEBES(ISYMAB)
            IAABS=MAGEB(1,IABABS)
            IBABS=MAGEB(2,IABABS)
          ELSE
            IABABS=IAB+NAGTBES(ISYMAB)
            IAABS=MAGTB(1,IABABS)
            IBABS=MAGTB(2,IABABS)
          END IF
          IP=IINAIS(IIABS)
          IQ=IEXTIS(IAABS)
          IR=IEXTIS(IBABS)
          RETURN
        END IF
        IIAB=IIAB-NIAB
      END DO
      WRITE(6,*)'NSIND BJAT: Impossible situation.'
      CALL ABEND()

      END
