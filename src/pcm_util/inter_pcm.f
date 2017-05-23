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
      Subroutine Inter_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NS,I,
     +                 IPRINT)
      IMPLICIT REAL*8 (A-H,O-Z)

      Dimension XE(*),YE(*),ZE(*),RE(*)
      DIMENSION P1(3),P2(3),P3(3),P4(3)
C
C     Trova il punto P4, sull`arco P1-P2 sotteso dal centro P3, che
C     si trova sulla superficie della sfera NS
C     P4 e' definito come combinazioe lineare di P1 e P2, con
C     il parametro ALPHA ottimizzato per tentativi.
C

      R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
      R = sqrt(R2)
      TOL = 1.D-12
      ALPHA = 0.5D0
      DELTA = 0.D0
      M = 1
  10  CONTINUE
      IF(M.GT.100) then
        IF(IPRINT.GT.0)
     +  Write(6,'(/,10X,'' INTER: too many iterations'')')
        Return
      EndIf
      ALPHA = ALPHA + DELTA
      DNORM = 0.D0
      DO 2000 JJ = 1,3
      P4(JJ)=P1(JJ)+ALPHA*(P2(JJ)-P1(JJ))-P3(JJ)
      DNORM = DNORM + P4(JJ)**2
 2000 continue
      DNORM = sqrt(DNORM)
      DO 2010 JJ = 1,3
      P4(JJ)= P4(JJ)*R/DNORM + P3(JJ)
 2010 continue
      DIFF2=(P4(1)-XE(NS))**2 + (P4(2)-YE(NS))**2 + (P4(3)-ZE(NS))**2
      DIFF = sqrt(DIFF2) - RE(NS)
      IF(abs(DIFF).LT.TOL) RETURN
      IF(I.EQ.0) THEN
      IF(DIFF.GT.0.D0) DELTA = 1.D0/(2.D0**(M+1))
      IF(DIFF.LT.0.D0) DELTA = - 1.D0/(2.D0**(M+1))
      M = M + 1
      GOTO 10
      ENDIF
      IF(I.EQ.1) THEN
      IF(DIFF.GT.0.D0) DELTA = - 1.D0/(2.D0**(M+1))
      IF(DIFF.LT.0.D0) DELTA = 1.D0/(2.D0**(M+1))
      M = M + 1
      GOTO 10
      ENDIF
      END
