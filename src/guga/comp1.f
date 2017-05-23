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
************************************************************************
      SUBROUTINE COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
      DIMENSION ITAI(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
*
      JO(L)=ICUNP(ICASE,L)
*
      CALL QENTER('COMP1')
      FAC=D1
      ITAIL=IX(IT2+LJ)
      DO 90 IN=1,ITAIL
      IC1=ICOUP(1)+IN
      JND1=JNDX(II+IC1)
      IF(JND1.EQ.0)GO TO 90
      IN2=ITAI(IN)
      IF(IN2.EQ.0)GO TO 90
      IC2=ICOUP1(1)+IN2
      IF(ITYP.NE.1)GO TO 91
      KK1=(JND1-1)*LN+L
      JOJ=JO(KK1)
      IF(JOJ.GT.1)JOJ=JOJ-1
      FAC=JOJ
      IF(JOJ.EQ.0)GO TO 90
91    IC1=JND1-IID
      JND2=JNDX(JJ+IC2)
      IF(JND2.EQ.0)GO TO 90
      IC2=JND2-JJD
      IOUT=IOUT+1
      KTYP=JTYP
      IF(JTYP.LE.3)GO TO 92
      ICT=IC1
      IC1=IC2
      IC2=ICT
      KTYP=JTYP-3
92    CONTINUE
*      IND=KTYP+2**6*IC2
      IND=IOR(KTYP,ISHFT(IC2,6))
*      ICOP1(IOUT)=IND+2**19*IC1
      ICOP1(IOUT)=IOR(IND,ISHFT(IC1,19))
      COP(IOUT)=FAC*COUP(1)
      IF(IOUT.LT.NBUF)GO TO 90
      ICOP1(nCOP+1)=NBUF
      IF(.FALSE.) THEN
         WRITE(6,*) 'WRITING BUFFER IN COMP1'
         WRITE(6,*) '======================='
         WRITE(6,'(A,I6)') 'ICOP1(nCOP+1) ',ICOP1(nCOP+1)
         WRITE(6,'(A,I6)') 'IADD10     ',IADD10
      END IF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
90    CONTINUE
      CALL QEXIT('COMP1')
      RETURN
      END
