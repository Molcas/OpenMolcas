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
      use definitions, only: iwp, wp
      use caspt2_module, only: NASHT, NACTEL
      IMPLICIT None

      integer(kind=iwp), intent(in):: NOP2, NOP3
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT),OP2(NOP2)
      real(kind=wp), intent(in):: OP3(NOP3)

      integer(kind=iwp) I,J,K,L,IJ,KL,M,N,MN,IJKLMN,IL,IN,KN,IND,IJKL
      real(kind=wp) X
C Purpose: Modify the coefficients in OP1 and OP2, using the
C input values of OP2 and OP3, so that the operators can be
C calculated using products of elementary excitation
C operators rather than normal-ordered products.


      IF (NACTEL>2) THEN

      DO I=1,NASHT
       DO J=1,NASHT
        IJ=I+NASHT*(J-1)
        DO K=1,NASHT
         DO L=1,NASHT
          KL=K+NASHT*(L-1)
          IF(KL.GT.IJ) CYCLE

          DO M=1,NASHT
           DO N=1,NASHT
            MN=M+NASHT*(N-1)
            IF(MN.GT.KL) CYCLE
            IJKLMN=((IJ+1)*IJ*(IJ-1))/6+(KL*(KL-1))/2+MN
            X=OP3(IJKLMN)
            IF(ABS(X).LT.1.0D-15) CYCLE

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

           END DO
          END DO

         END DO
        END DO
       END DO
      END DO

      END IF

      IF(NACTEL>1) THEN

      DO I=1,NASHT
       DO J=1,NASHT
        IJ=I+NASHT*(J-1)
        DO K=1,NASHT
         DO L=1,NASHT
          KL=K+NASHT*(L-1)
          IF(KL.GT.IJ) CYCLE

          IF(J.EQ.K) THEN
            IJKL=(IJ*(IJ-1))/2+KL
            OP1(I,L)=OP1(I,L)-OP2(IJKL)
          END IF

         END DO
        END DO
       END DO
      END DO

      END IF

      END SUBROUTINE MODOP
