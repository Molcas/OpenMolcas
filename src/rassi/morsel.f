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
      INTEGER FUNCTION MorsPop(IMORS)
      IMPLICIT NONE
      DIMENSION NUM(0:15)
      INTEGER I1,I2,I3,IMORS
      INTEGER J,NUM
      DATA NUM / 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 /
      MorsPop=0 ! dummy initialize
      IF(IMORS.LT.0) GOTO 99
      I1=IMORS/16
      J=IMORS-16*I1
      MorsPop=NUM(J)
      IF(I1.EQ.0) RETURN
      I2=I1/16
      J=I1-16*I2
      MorsPop=MorsPop+NUM(J)
      IF(I2.EQ.0) RETURN
      I3=I2/16
      J=I2-16*I3
      MorsPop=MorsPop+NUM(J)
      IF(I3.EQ.0) RETURN
  99  CONTINUE
      WRITE(6,*)' MorsPop: Bad IMORS=',IMORS
      CALL ABEND()
      END
      INTEGER FUNCTION MorsParity(IMORS)
      IMPLICIT NONE
      DIMENSION ISG(0:15)
      INTEGER I1,I2,I3,IMORS
      INTEGER J,ISG
      DATA ISG / 1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1 /
      MorsParity=0 ! dummy initialize
      IF(IMORS.LT.0) GOTO 99
      I1=IMORS/16
      J=IMORS-16*I1
      MorsParity=ISG(J)
      IF(I1.EQ.0) RETURN
      I2=I1/16
      J=I1-16*I2
      MorsParity=MorsParity*ISG(J)
      IF(I2.EQ.0) RETURN
      I3=I2/16
      J=I2-16*I3
      MorsParity=MorsParity*ISG(J)
      IF(I3.EQ.0) RETURN
  99  CONTINUE
      WRITE(6,*)' MorsParity: Bad IMORS=',IMORS
      CALL ABEND()
      END
      INTEGER FUNCTION MorsSpin(IMORS,MS2ARR)
      IMPLICIT NONE
      INTEGER MORSBITS
      PARAMETER (MORSBITS=8)
      DIMENSION MS2ARR(*)
      INTEGER I,IB,IBIT,IMORS,MS2ARR
      MorsSpin=0
      IF(IMORS.LT.0) GOTO 99
      IB=IMORS
      DO I=1,MORSBITS
        IBIT=MOD(IB,2)
        IB=IB/2
        IF(IBIT.EQ.1) MorsSpin=MorsSpin+MS2ARR(I)
      END DO
      RETURN
  99  CONTINUE
      WRITE(6,*)' MorsSpin: Bad IMORS=',IMORS
      CALL ABEND()
      END
      INTEGER FUNCTION MorsSymm(IMORS,ISMARR)
      IMPLICIT NONE
      INTEGER MORSBITS
      PARAMETER (MORSBITS=8)
      DIMENSION ISMARR(*)
      INTEGER I,IB,IBIT,IMORS,ISMARR
#include "symmul.fh"
      MorsSymm=1
      IF(IMORS.LT.0) GOTO 99
      IB=IMORS
      DO I=1,MORSBITS
        IBIT=MOD(IB,2)
        IB=IB/2
        IF(IBIT.EQ.1) MorsSymm=MUL(MorsSymm,ISMARR(I))
      END DO
      RETURN
  99  CONTINUE
      WRITE(6,*)' MorsSymm: Bad IMORS=',IMORS
      CALL ABEND()
      END
      SUBROUTINE MorsWrite(IMORS,STRING)
      IMPLICIT NONE
      CHARACTER*(*) STRING
      INTEGER I,IB,IBIT,IMORS
      IF(IMORS.LT.0) GOTO 99
      IB=IMORS
      DO I=1,LEN(STRING)
        STRING(I:I)='0'
        IBIT=MOD(IB,2)
        IB=IB/2
        IF(IBIT.EQ.1) STRING(I:I)='1'
      END DO
      IF(IB.GT.0) THEN
       DO I=1,LEN(STRING)
        STRING(I:I)='*'
       END DO
      END IF
      RETURN
  99  CONTINUE
      WRITE(6,*)' MorsWrite: Bad IMORS=',IMORS
      CALL ABEND()
      END
*      SUBROUTINE MorsRead(IMORS,STRING)
*      IMPLICIT NONE
*      CHARACTER*(*) STRING
*      INTEGER I,POW2,IMORS
*      IMORS=0
*      POW2=1
*      DO I=1,LEN(STRING)
*        IF(STRING(I:I).EQ.'1') IMORS=IMORS+POW2
*        POW2=2*POW2
*      END DO
*      RETURN
*      END
      INTEGER FUNCTION Occ2Mrs(NO,IARRAY)
      IMPLICIT NONE
      INTEGER NO
      INTEGER IARRAY(NO)
      INTEGER I,POW2
      OCC2MRS=0
      POW2=1
      DO I=1,NO
        IF(IARRAY(I).NE.0) OCC2MRS=OCC2MRS+POW2
        POW2=2*POW2
      END DO
      RETURN
      END
      INTEGER FUNCTION MorsCre(IMORS,IPOS)
      IMPLICIT NONE
      INTEGER MorsParity
      EXTERNAL MorsParity
      INTEGER IMORS,IPOS,ISGN,MASK
      MorsCre=999999
      MASK=2**(IPOS-1)
      IF(IAND(MASK,IMORS).NE.0) RETURN
      ISGN=MorsParity(IMORS/MASK)
      MorsCre=ISGN*(IMORS+MASK)
      RETURN
      END
      INTEGER FUNCTION MorsAnn(IMORS,IPOS)
      IMPLICIT NONE
      INTEGER MorsParity
      EXTERNAL MorsParity
      INTEGER IMORS,IPOS,ISGN,MASK
      MorsAnn=999999
      MASK=2**(IPOS-1)
      IF(IAND(MASK,IMORS).EQ.0) RETURN
      MorsAnn=IMORS-MASK
      ISGN=MorsParity(MorsAnn/MASK)
      MorsAnn=ISGN*MorsAnn
      RETURN
      END
