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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE IIJJ_CPF(ICASE,JSY,HDIAG,FC,FIJ,FJI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),HDIAG(*),FC(*),FIJ(*),FJI(*)
      DIMENSION ICASE(*)
      DIMENSION IOC(55)
*
      JO(L)=ICUNP(ICASE,L)
      JSYM(L)=JSUNP_CPF(JSY,L)
*
      CALL QENTER('IIJJ_CPF')
      IAD27=0
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IRL=IRC(ILIM)

      DO 100 IR=1,IRL
        DO I=1,LN
          JOJ=JO(I+LN*(IR-1))
          IOC(I)=(JOJ+1)/2
        END DO
        NSS=MUL(JSYM(IR),LSYM)

        SUM=0.0D0
        DO I=1,LN
          IF(IOC(I).NE.0) THEN
            DO J=1,I-1
              IJ=(I*(I-1))/2+J
              IF(IOC(J).NE.0) THEN
                TERM=IOC(I)*(IOC(J)*FIJ(IJ)-FJI(IJ))
                SUM=SUM+TERM
              END IF
            END DO
            II=(I*(I+1))/2
            TERM=(IOC(I)-1)*FIJ(II)+IOC(I)*FC(II)
            SUM=SUM+TERM
          END IF
        END DO

        IF(IR.GT.IRC(1))GO TO 120
C IR=1..IRC(1), HDIAG(IR)=SUM
        HDIAG(IR)=SUM
        IF(IR.EQ.IRC(1))CALL dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)
        GO TO 100

120     CONTINUE
        IF(IR.GT.IRC(2))GO TO 130
C IR=IRC(1)+1 ... IRC(2)
        IND=0
        NA1=NSYS(NSS)+1
        NA2=NSYS(NSS+1)
        IF(NA2.LT.NA1)GO TO 100
        DO NA=NA1,NA2
          IND=IND+1
          IA=IROW(LN+NA)
          SUM1=SUM+FC(IA+LN+NA)
          DO I=1,LN
          IF(IOC(I).NE.0) SUM1=SUM1+
     &                  IOC(I)*FIJ(IA+I)-FJI(IA+I)
          END DO
          HDIAG(IND)=SUM1
        END DO
        CALL dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
        GO TO 100

130     CONTINUE
C IR=IRC(2)+1 ... IRC(ILIM)
        IND=0
        DO NA=1,NVIRT
          NSA=MUL(NSS,NSM(LN+NA))
          NB1=NSYS(NSA)+1
          NB2=NSYS(NSA+1)
          IF(NB2.GT.NA)NB2=NA
          IF(NB2.LT.NB1)GO TO 141
          IA=IROW(LN+NA)
          IAV=IA+LN
          DO NB=NB1,NB2
            IND=IND+1
            IB=IROW(LN+NB)
            IBV=IB+LN
            TERM=SUM+FIJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
            IF(IR.LE.IRC(3))SUM1=TERM-FJI(IAV+NB)
            IF(IR.GT.IRC(3))SUM1=TERM+FJI(IAV+NB)
            DO I=1,LN
              IF(IOC(I).NE.0) THEN
                TERM=IOC(I)*(FIJ(IA+I)+FIJ(IB+I))-FJI(IA+I)-FJI(IB+I)
                SUM1=SUM1+TERM
              END IF
            END DO
            HDIAG(IND)=SUM1
          END DO
141       CONTINUE
        END DO
        IF(IND.GT.0)CALL dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
100   CONTINUE
      CALL QEXIT('IIJJ_CPF')
      RETURN
      END
