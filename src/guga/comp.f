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
      SUBROUTINE COMP(I,LJ,ITYP,L,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
*
      JO(L)=ICUNP(ICASE,L)
*
      IF (IT1.NE.IT2) Then
         Write (6,*) 'Comp: IT1.NE.IT2'
         Write (6,*) 'IT1,IT2=',IT1,IT2
         Write (6,*) 'ITYP,L=',ITYP,L
         Call Abend
      End If
      FAC=D1
      KM=I
25    KM=KM-1
      IF(KM.EQ.0)GO TO 26
      IWAY(KM)=1
27    CALL PATH(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 25
      KM=KM+1
      IF(KM.EQ.I)GO TO 71
      GO TO 27
26    IVL=J2(1)
      ITAIL=IX(IT1+LJ)
      IVV=IV0-IVL
      JJ=0
      JJD=0
      IF(IVV.NE.0)JJ=IRC(IVV)
      IF(IVV.NE.0)JJD=JRC(IVV)
      DO 90 IN=1,ITAIL
      IC1=ICOUP(1)+IN
      IC2=ICOUP1(1)+IN
      JND1=JNDX(JJ+IC1)
      IF(JND1.EQ.0)GO TO 90
      IF(ITYP.NE.1)GO TO 91
      II1=(JND1-1)*LN+L
      JOJ=JO(II1)
      IF(JOJ.GT.1)JOJ=JOJ-1
      FAC=JOJ
      IF(JOJ.EQ.0)GO TO 90
91    IC1=JND1-JJD
      JND2=JNDX(JJ+IC2)
      IF(JND2.EQ.0)GO TO 90
      IC2=JND2-JJD
      IOUT=IOUT+1
      IVL0=IV0-IVL
*      IND=IVL0+2**6*IC2
*      ICOP1(IOUT)=IND+2**19*IC1
      IND=IOR(IVL0,ISHFT(IC2,6))
      ICOP1(IOUT)=IOR(IND,ISHFT(IC1,19))

      COP(IOUT)=FAC*COUP(I)
      IF(IOUT.LT.NBUF)GO TO 90
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
90    CONTINUE
      IF(I.EQ.1)GO TO 71
      KM=1
      GO TO 27
*
71      Continue
      Return
      End
