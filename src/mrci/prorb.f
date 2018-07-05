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
      SUBROUTINE PRORB(CNO,OCC)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION CNO(NCMO),OCC(NBAST)
      CHARACTER*(LENIN8) CLEAN_BNAME
      EXTERNAL CLEAN_BNAME
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)'NATURAL ORBITALS IN AO BASIS. IN EACH SYMMETRY,'
      CALL XFLUSH(6)
      WRITE(6,*)'THE ORBITALS PRINTED ARE THOSE UP TO AND INCLUDING'
      CALL XFLUSH(6)
      WRITE(6,*)'THE LAST ORBITAL WITH OCCUPATION NUMBER LARGER'
      CALL XFLUSH(6)
      WRITE(6,'(A,F10.7)')' THAN THRORB = ',THRORB
      CALL XFLUSH(6)
      IEB=0
      IEM=0
      NDIV=10
      DO 100 ISYM=1,NSYM
         NB=NBAS(ISYM)
         IF(NB.EQ.0) GO TO 100
         NPRT=0
         DO 10 I=1,NB
            IF(OCC(IEB+I).GE.THRORB) NPRT=I
10       CONTINUE
         IF(NPRT.EQ.0) GO TO 40
         WRITE(6,'(/28X,''SYMMETRY LABEL'',I3)') ISYM
      CALL XFLUSH(6)
         DO 30 IST=1,NPRT,NDIV
           IEND=MIN(NPRT,IST-1+NDIV)
           WRITE(6,'(/5X,''ORBITAL'',6X,10I8)') (I,I=IST,IEND)
      CALL XFLUSH(6)
           WRITE(6,'( 5X,''OCC.NO.'',8X,10F8.5)')
     *              (OCC(IEB+I),I=IST,IEND)
      CALL XFLUSH(6)
           WRITE(6,*)
      CALL XFLUSH(6)
           DO 20 I=1,NB
              JSMO=IEM+I+NB*(IST-1)
              JEMO=IEM+I+NB*(IEND-1)
              WRITE(6,'(1X,I3,2X,A,10F8.4)')
     *                   I,CLEAN_BNAME(NAME(IEB+I),LENIN),
     *                   (CNO(J),J=JSMO,JEMO,NB)
      CALL XFLUSH(6)
20         CONTINUE
30       CONTINUE
40       CONTINUE
         IEB=IEB+NB
         IEM=IEM+NB*NB
100   CONTINUE
      RETURN
      END
