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
      SUBROUTINE MODOP(OP1,NOP2,OP2,NOP3,OP3)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2),OP3(NOP3)

C Purpose: Modify the coefficients in OP1 and OP2, using the
C input values of OP2 and OP3, so that the operators can be
C calculated using products of elementary excitation
C operators rather than normal-ordered products.


      IF(NACTEL.LE.2) GOTO 100

      DO I=1,NASHT
       DO J=1,NASHT
        IJ=I+NASHT*(J-1)
        DO K=1,NASHT
         DO L=1,NASHT
          KL=K+NASHT*(L-1)
          IF(KL.GT.IJ) GOTO 120
          DO M=1,NASHT
           DO N=1,NASHT
            MN=M+NASHT*(N-1)
            IF(MN.GT.KL) GOTO 110
            IJKLMN=((IJ+1)*IJ*(IJ-1))/6+(KL*(KL-1))/2+MN
            X=OP3(IJKLMN)
            IF(ABS(X).LT.1.0D-15) GOTO 110

            IF(K.EQ.J) THEN
              IL=I+NASHT*(L-1)
              IF(IL.GE.MN) THEN
                IND=(IL*(IL-1))/2+MN
              ELSE
                IND=(MN*(MN-1))/2+IL
              END IF
              OP2(IND)=OP2(IND)-X
              IF(M.EQ.L) THEN
                OP1(I,N)=OP1(I,N)-X
              END IF
            END IF

            IF(M.EQ.J) THEN
              IN=I+NASHT*(N-1)
              IF(IN.GE.KL) THEN
                IND=(IN*(IN-1))/2+KL
              ELSE
                IND=(KL*(KL-1))/2+IN
              END IF
              OP2(IND)=OP2(IND)-X
            END IF

            IF(M.EQ.L) THEN
              KN=K+NASHT*(N-1)
              IF(KN.GE.IJ) THEN
                IND=(KN*(KN-1))/2+IJ
              ELSE
                IND=(IJ*(IJ-1))/2+KN
              END IF
              OP2(IND)=OP2(IND)-X
            END IF

 110        CONTINUE
           END DO
          END DO
 120      CONTINUE
         END DO
        END DO
       END DO
      END DO

 100  CONTINUE
      IF(NACTEL.LE.1) GOTO 200

      DO I=1,NASHT
       DO J=1,NASHT
        IJ=I+NASHT*(J-1)
        DO K=1,NASHT
         DO L=1,NASHT
          KL=K+NASHT*(L-1)
          IF(KL.GT.IJ) GOTO 10

          IJKL=(IJ*(IJ-1))/2+KL
          IF(J.EQ.K) THEN
            OP1(I,L)=OP1(I,L)-OP2(IJKL)
          END IF

  10      CONTINUE
         END DO
        END DO
       END DO
      END DO

 200  CONTINUE
      RETURN
      END
