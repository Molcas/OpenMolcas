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
      SUBROUTINE SXTYP(NSXTP,ITP,JTP,LEL1,LEL3,REL1,REL3)
*
* Types of creation and annihilation  operators so
* <L!a+ a!R> is nonvanishing
*
* L is defined by LEL1,LEL3
* R is defined by REL1,REL3
*
      INTEGER REL1,REL3
      INTEGER ITP(3),JTP(3)
      NSXTP = 0
*
*. To get rid of annoying and incorrect compiler warnings
      I1 = 0
      I3 = 0
      IJ1 = 0
      IJ3 = 0
*
      DO 100 I123 = 1, 3
        IF(I123.EQ.1) THEN
          I1 = 1
          I3 = 0
        ELSE IF(I123.EQ.2) THEN
          I1 = 0
          I3 = 0
        ELSE IF(I123.EQ.3) THEN
          I1 = 0
          I3 = 1
        END IF
        IF(LEL1-I1.LT.0) GOTO 100
        IF(LEL3-I3.LT.0) GOTO 100
        DO 50 J123 = 1, 3
          IF(J123.EQ.1) THEN
            IJ1 = I1 - 1
            IJ3 = I3
          ELSE IF(J123.EQ.2) THEN
            IJ1 = I1
            IJ3 = I3
          ELSE IF(J123.EQ.3) THEN
            IJ1 = I1
            IJ3 = I3-1
          END IF
          IF(REL1+IJ1.EQ.LEL1.AND.REL3+IJ3.EQ.LEL3) THEN
            NSXTP = NSXTP + 1
            ITP(NSXTP) = I123
            JTP(NSXTP) = J123
          END IF
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)')
     &  ' SX  connecting LEL1,LEL3,REL1,REL3 ',LEL1,LEL3,REL1,REL3
        WRITE(6,*) ' Number of connections obtained ', NSXTP
        WRITE(6,*) ' ITYPE JTYPE '
        WRITE(6,*) ' =========== '
        DO 200 I = 1, NSXTP
         WRITE(6,'(2I5)') ITP(I),JTP(I)
  200   CONTINUE
*
      END IF
*
      RETURN
      END
