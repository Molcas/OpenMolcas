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
cpgi$g opt=1
      SUBROUTINE MKHREF(HREF,FC,FIJKL,JREFX)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION HREF(*),FC(*),FIJKL(*),JREFX(NCVAL)
*
      NHREF=(NREF*(NREF+1))/2
      CALL FZERO(HREF,NHREF)
      ICHK=0
      IK=0
      FINI=0.0D00
      IADD25=0
      CALL dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
      IADD10=IAD10(8)
100   CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0) GOTO 200
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0) THEN
        ICHK=0
        INDI=IND
*        NI=MOD(INDI,2**10)
*        NK=MOD(INDI/2**10,2**10)
        NI=IBITS(INDI, 0,10)
        NK=IBITS(INDI,10,10)
        IK=IROW(NK)+NI
        FINI=FC(IK)
        GO TO 10
      END IF
      IF(IND.EQ.0) THEN
        ICHK=1
        GO TO 10
      END IF
*      IVL=MOD(IND,2**6)
      IVL=IBITS(IND, 0, 6)
      IF(IVL.NE.IVVER)GO TO 10
*      IC2=MOD(IND/2**6,2**13)
      IC2=IBITS(IND, 6,13)
      NA=JREFX(IC2)
      IF(NA.EQ.0)GO TO 10
*      IC1=MOD(IND/2**19,2**13)
      IC1=IBITS(IND,19,13)
      NB=JREFX(IC1)
      IF(NB.EQ.0)GO TO 10
      IF(NA.LT.NB) THEN
        NAT=NA
        NA=NB
        NB=NAT
      END IF
      IVEC=(NA*(NA-1))/2+NB
      HREF(IVEC)=HREF(IVEC)+COP(IN)*FINI
10    CONTINUE
      GO TO 100
200   CONTINUE
      ICHK=0
      NIJ=IROW(LN+1)
      NIJKL=NIJ*(NIJ+1)/2
      CALL FZERO(FIJKL,NIJKL)
      FINI=0.0D00
      IADR=LASTAD(1)
201   CONTINUE
      CALL dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
      CALL iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
      LENGTH=INDSRT(NSRTMX+1)
      IADR=INDSRT(NSRTMX+2)
*      IF(LENGTH.GT.0) CALL SCATTER(LENGTH,FIJKL,INDSRT,VALSRT)
      do i=1,length
        FIJKL(INDSRT(i))=VALSRT(i)
      end do
      IF(IADR.NE.-1) GO TO 201
      IADD10=IAD10(5)
300   CONTINUE
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 300
      IF(LEN.LT.0)GOTO 400
      DO 310 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0) THEN
        ICHK=0
        INDI=IND
CPAM96        IP=IAND(INDI,255)
CPAM96        JP=IAND(ISHFT(INDI,-8),255)
CPAM96        KP=IAND(ISHFT(INDI,-16),255)
CPAM96        LP=IAND(ISHFT(INDI,-24),255)
*        IP=MOD(INDI,2**8)
*        JP=MOD(INDI/2**8,2**8)
*        KP=MOD(INDI/2**16,2**8)
*        LP=MOD(INDI/2**24,2**8)
        IP=IBITS(INDI, 0,8)
        JP=IBITS(INDI, 8,8)
        KP=IBITS(INDI,16,8)
        LP=IBITS(INDI,24,8)
        NIJ=IROW(IP)+JP
        NKL=IROW(KP)+LP
        IND=NIJ*(NIJ-1)/2+NKL
        FINI=FIJKL(IND)
        GOTO 310
      END IF
      IF(IND.EQ.0) THEN
        ICHK=1
        GOTO 310
      END IF
*      IVL=MOD(IND,2**6)
      IVL=IBITS(IND, 0, 6)
      IF(IVL.NE.0)GO TO 310
*      IC2=MOD(IND/2**6,2**13)
      IC2=IBITS(IND, 6,13)
      NA=JREFX(IC2)
      IF(NA.EQ.0)GO TO 310
*      IC1=MOD(IND/2**19,2**13)
      IC1=IBITS(IND,19,13)
      NB=JREFX(IC1)
      IF(NB.EQ.0)GO TO 310
      IF(NA.LT.NB) THEN
        NAT=NA
        NA=NB
        NB=NAT
      END IF
      IVEC=(NA*(NA-1))/2+NB
      HREF(IVEC)=HREF(IVEC)+COP(IN)*FINI
310   CONTINUE
      GO TO 300
400   CONTINUE
      IADD25=IAD25S
      IBUF=nCOP
      DO 410 I=1,IRC(1)
        IBUF=IBUF+1
        IF(IBUF.GT.nCOP) THEN
          CALL dDAFILE(Lu_25,2,COP,nCOP,IADD25)
          IBUF=1
        END IF
        IR=JREFX(I)
        IF(IR.GT.0) THEN
          IIR=(IR*(IR+1))/2
          HREF(IIR)=HREF(IIR)+COP(IBUF)+POTNUC
        END IF
410   CONTINUE
      RETURN
      END
