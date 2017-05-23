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
      SUBROUTINE DXTYP(NDXTP,ITYP,JTYP,KTYP,LTYP,LEL1,LEL3,REL1,REL3)
*
* Obtain types of I,J,K,l so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(L).ge.type(J)
* are included
*
      INTEGER REL1,REL3
      INTEGER ITYP(6),JTYP(6),KTYP(6),LTYP(6)
*
*. To get rid of annoying and incorrect compiler warnings
      I1 = 0
      I3 = 0
      IK1 = 0
      IK3 = 0
      IKL1 = 0
      IKL3 = 0
      IKLJ1 = 0
      IKLJ3 = 0
*
      NDXTP = 0
      DO 400 ITP = 1, 3
        IF(ITP.EQ.1) THEN
          I1 = 1
          I3 = 0
        ELSE IF(ITP.EQ.2) THEN
          I1 = 0
          I3 = 0
        ELSE IF(ITP.EQ.3) THEN
          I1 = 0
          I3 = 1
        END IF
        DO 300 KTP = 1, ITP
          IF(KTP.EQ.1) THEN
            IK1 = I1+1
            IK3 = I3
          ELSE IF(KTP.EQ.2) THEN
            IK1 = I1
            IK3 = I3
          ELSE IF(KTP.EQ.3) THEN
            IK1 = I1
            IK3 = I3+1
          END IF
          IF(LEL1-IK1.LT.0) GOTO 300
          IF(LEL3-IK3.LT.0) GOTO 300
          DO 200 LTP = 1,3
            IF(LTP.EQ.1) THEN
              IKL1 = IK1-1
              IKL3 = IK3
            ELSE IF(LTP.EQ.2) THEN
              IKL1 = IK1
              IKL3 = IK3
            ELSE IF(LTP.EQ.3) THEN
              IKL1 = IK1
              IKL3 = IK3-1
            END IF
            DO 100 JTP = 1, 3
              IF(JTP.EQ.1) THEN
                IKLJ1 = IKL1-1
                IKLJ3 = IKL3
              ELSE IF(JTP.EQ.2) THEN
                IKLJ1 = IKL1
                IKLJ3 = IKL3
              ELSE IF(JTP.EQ.3) THEN
                IKLJ1 = IKL1
                IKLJ3 = IKL3-1
              END IF
              IF(IKLJ1+REL1.EQ.LEL1.AND.IKLJ3+REL3.EQ.LEL3) THEN
                NDXTP = NDXTP + 1
                ITYP(NDXTP) = ITP
                KTYP(NDXTP) = KTP
                LTYP(NDXTP) = LTP
                JTYP(NDXTP) = JTP
              END IF
  100       CONTINUE
  200     CONTINUE
  300   CONTINUE
  400 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)')
     &  ' Double excitations connecting LEL1,LEL3,LEL1,LEL3 ',
     &    LEL1,LEL3,REL1,REL3
        WRITE(6,*) '  ITYP KTYP LTYP JTYP '
        WRITE(6,*) '  ===================='
        DO 10 IDX = 1,NDXTP
          WRITE(6,'(1H ,5I5)')ITYP(IDX),KTYP(IDX),LTYP(IDX),JTYP(IDX)
   10   CONTINUE
      END IF
*
      RETURN
      END
