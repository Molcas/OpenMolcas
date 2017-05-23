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
!EB      Subroutine Get_OrbCen(nPrim,nBas,NORBI,Q_MltPl,RCHC,
      Subroutine Get_OrbCen(nPrim,NORBI,Q_MltPl,RCHC,
     &CENTX,CENTY,CENTZ,OCOF)
      Implicit Real*8 (a-h,o-z)


      Dimension RCPO(3,NORBI), RCMI(3,NORBI)
      Dimension CENTX(nPrim*(nPrim+1)/2)
      Dimension CENTY(nPrim*(nPrim+1)/2)
      Dimension CENTZ(nPrim*(nPrim+1)/2)
      Dimension Q_MltPl(nPrim*(nPrim+1)/2)
      Dimension RCHC(3,NORBI)
      Dimension oCof(NORBI,nPrim)

C
C CALCULATE A CENTER OF CHARGE FOR EACH MOLECULAR ORBITAL
C
      Do I=1,NORBI
         RCPO(1,I) = 0.0
         RCPO(2,I) = 0.0
         RCPO(3,I) = 0.0
         RCMI(1,I) = 0.0
         RCMI(2,I) = 0.0
         RCMI(3,I) = 0.0
         QPOS = 0.0
         QMIN = 0.0
         Do J=1,nPrim
            Do K=1,J
C
C THE WEIGHTED AVERAGE OF THE POSITIVE AND NEGATIVE CONTRIBUTIONS
C
               OOQ = OCOF(I,J)*OCOF(I,K)*Q_MltPl(J*(J-1)/2+K)*2.0
               IF( OOQ .GE. 0.0 ) THEN
                  QPOS = QPOS + OOQ
                  RCPO(1,I) = RCPO(1,I) + OOQ*CENTX(J*(J-1)/2+K)
                  RCPO(2,I) = RCPO(2,I) + OOQ*CENTY(J*(J-1)/2+K)
                  RCPO(3,I) = RCPO(3,I) + OOQ*CENTZ(J*(J-1)/2+K)
               ELSE
                  QMIN = QMIN + OOQ
                  RCMI(1,I) = RCMI(1,I) + OOQ*CENTX(J*(J-1)/2+K)
                  RCMI(2,I) = RCMI(2,I) + OOQ*CENTY(J*(J-1)/2+K)
                  RCMI(3,I) = RCMI(3,I) + OOQ*CENTZ(J*(J-1)/2+K)
               ENDIF
            EndDo
            OOQ = OCOF(I,J)*OCOF(I,J)*Q_MltPl(J*(J+1)/2)
            IF( OOQ .GE. 0.0 ) THEN
               QPOS = QPOS - OOQ
               RCPO(1,I) = RCPO(1,I) - OOQ*CENTX(J*(J+1)/2)
               RCPO(2,I) = RCPO(2,I) - OOQ*CENTY(J*(J+1)/2)
               RCPO(3,I) = RCPO(3,I) - OOQ*CENTZ(J*(J+1)/2)
            ELSE
               QMIN = QMIN - OOQ
               RCMI(1,I) = RCMI(1,I) - OOQ*CENTX(J*(J+1)/2)
               RCMI(2,I) = RCMI(2,I) - OOQ*CENTY(J*(J+1)/2)
               RCMI(3,I) = RCMI(3,I) - OOQ*CENTZ(J*(J+1)/2)
            ENDIF
         EndDo
         RCHC(1,I) = (RCPO(1,I) - RCMI(1,I))/(QPOS-QMIN)
         RCHC(2,I) = (RCPO(2,I) - RCMI(2,I))/(QPOS-QMIN)
         RCHC(3,I) = (RCPO(3,I) - RCMI(3,I))/(QPOS-QMIN)
      EndDo
!EB 96   FORMAT(I5,3E15.8)
      Return
      End
