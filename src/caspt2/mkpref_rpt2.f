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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKPREF_RPT2(N,G2,PREF,NPREF)
      use definitions, only: iwp, wp
      IMPLICIT NONE

      INTEGER(kind=iwp), INTENT(IN) :: N, NPREF
      REAL(kind=wp), INTENT(IN) ::  G2(N,N,N,N)
      REAL(kind=wp), INTENT(OUT) :: PREF(NPREF)

      INTEGER(kind=iwp) I,J,K,L,IJ,JI,KL,LK
      INTEGER(kind=iwp) IJKL,IJLK,JIKL,JILK
      INTEGER(kind=iwp) IJT,KLT,IJKLT

      REAL(kind=wp) P1,P2

C Compute PREF(PQRS) = <0| 0.5*Epqrs |0>
C from G2(P,Q,R,S) = <0| Epqrs |0>
C Storage differs: PREF is triangular
C in the Fortran-like indices PQ, RS.

      IJT=0
      IJKLT=0
      DO I=1,N
        DO J=1,I
          IJT=IJT+1
          IJ=I+N*(J-1)
          JI=J+N*(I-1)
          KLT=0
          DO K=1,N
            DO L=1,K
              KLT=KLT+1
              IF(KLT.GT.IJT) GOTO 130
              IJKLT=IJKLT+1
              KL=K+N*(L-1)
              LK=L+N*(K-1)

              P1=0.5D0*G2(I,J,K,L)
              P2=0.5D0*G2(I,J,L,K)
              IF(J.GE.L) THEN
                IJKL=(IJ*(IJ-1))/2+KL
              ELSE
                IJKL=(KL*(KL-1))/2+IJ
              END IF
              IF(J.GE.K) THEN
                IJLK=(IJ*(IJ-1))/2+LK
              ELSE
                IJLK=(LK*(LK-1))/2+IJ
              END IF
              JIKL=(JI*(JI-1))/2+KL
              JILK=(JI*(JI-1))/2+LK
              PREF(IJKL)=P1
              PREF(IJLK)=P2
              PREF(JIKL)=P2
              PREF(JILK)=P1
            END DO
          END DO
 130    CONTINUE
        END DO
      END DO

      END SUBROUTINE MKPREF_RPT2
