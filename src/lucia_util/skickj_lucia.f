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
* Copyright (C) 1994,1997, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE SKICKJ_LUCIA(   SKII,   CKJJ,    NKA,    NKB,  XIJKL,
     &                             NI,     NJ,     NK,     NL,   MAXK,
     &                           KBIB,  XKBIB,   KBJB,  XKBJB,  IKORD,
     &                           FACS, IROUTE)
*
*
* Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
*
*
* Jeppe Olsen, Spring of 94
*
* : Note : Route 1 has retired, March 97
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(*)
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(*)
*. Scratch
#include "mxpdim.fh"
      DIMENSION XIJILS(MXPTSOB)
*. To get rid of annoying and incorrect compiler warnings
      JKINTOF = 0
      IKINTOF = 0
*
      IF(NI.GT.MXPTSOB.OR.NJ.GT.MXPTSOB.OR.NK.GT.MXPTSOB
     &   .OR.NL.GT.MXPTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : > MXPTSOB '
         WRITE(6,*) ' N, MXPTSOB ',MAX(NI,NJ,NK,NL),MXPTSOB
*         STOP ' Redim MXPTSOB '
        CALL SYSABENDMSG('lucia_util/skickj','Redim MXPTSOB',' ')
      END IF
*
      NTEST =  000
*
      IF(IROUTE.EQ.3) THEN
* S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(NTEST.GE.100)
     &           WRITE(6,*) ' KB,K,L,IB,JB',KB,K,L,IB,JB
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
*
                   IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
                     IMAX  = K
                     JKINTOF = INTOF + (K-1)*NJ
C                    CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                     DO J = L,NL
                       XIJILS(J) = XIJKL(JKINTOF-1+J)
                     END DO
                     XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
                     DO J = L+1, NL
                      XIJKL(JKINTOF-1+J) = 0.0D0
                     END DO
                   END IF
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF),
     &                         CKJJ(ICOFF),XIJKL(INTOF),NKA,IMAX,  NKA,
     &                              NJ,     NJ,   IMAX,   FACS,FACTOR ,
     &                               0)
                   IF(IKORD.NE.0) THEN
                      DO J = L,NL
                        XIJKL(JKINTOF-1+J) =  XIJILS(J)
                      END DO
C                    CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                   END IF
*
                 END IF
               END DO
*
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
      ELSE IF(IROUTE.EQ.2) THEN
* S(I,Ka,Ib) = S(I,Ka,Ib) + sum(j) (ij!kl) C(j,Ka,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
*
                   JMAX = NJ
                   IF(IKORD.NE.0) THEN
*. Restrict so (ji) .le. (kl)
                     JMAX  = K
                     IKINTOF = INTOF + (K-1)*NI
                     CALL COPVEC(XIJKL(IKINTOF),XIJILS,NI)
                     XIJKL(IKINTOF-1+L) = 0.5D0*XIJKL(IKINTOF-1+L)
                     DO I = L+1, NL
                      XIJKL(IKINTOF-1+I) = 0.0D0
                     END DO
                   END IF
*
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF),
     &                         XIJKL(INTOF),CKJJ(ICOFF),NI,NKA, NI,
     &                              NJ,     NJ,    NKA,   FACS, FACTOR,
     &                               0)
*
                 IF(IKORD.NE.0) THEN
                   CALL COPVEC(XIJILS,XIJKL(IKINTOF),NI)
                 END IF
*
                 END IF
               END DO
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*

      ELSE IF (IROUTE.EQ.1) THEN
         WRITE(6,*) ' Sorry route 1 has retired, March 1997'
*         STOP'SKICKJ:Invalid route=1'
        CALL SYSABENDMSG('lucia_util/skickj','Internal error',' ')
C     DO 1000 KB = 1, NKB
*. Number of nonvanishing a+lb !Kb>
C       LL = 0
C       DO L = 1, NL
C         IF(KBJB(KB,L).NE.0) LL = LL + 1
C       END DO
*
C       IKEFF = 0
C       DO 900 K = 1, NK
C         IB = KBIB(KB,K)
C         IF(IB.EQ.0) GOTO 900
C         SGNK = XKBIB(KB,K)
*
C         IF(IKORD.EQ.0) THEN
C            LI = NI
C            IMIN = 1
C         ELSE
C            LI = NI-K+1
C            IMIN = K
C         END IF
*
C         DO 700 I = IMIN, NI
C           IKEFF = IKEFF + 1
C           IOFF = (IKEFF-1)*NJ*LL
*. Offset for S(1,IB,i)
C           IBOFF(IKEFF)  = (I-1)*NIB+IB
C           LEFF = 0
C           DO 800 L = 1, NL
C             JB = KBJB(KB,L)
C             IF(JB.EQ.0) GOTO 800
C             LEFF = LEFF + 1
C             SGNL = XKBJB(KB,L)
C             IF(IKORD.EQ.1.AND.I.EQ.K)THEN
C                FACTOR = 0.5D0*SGNK*SGNL
C             ELSE
C                FACTOR =       SGNK*SGNL
C             END IF
C             JL0 = (LEFF-1)*NJ
C             JLIK0 = (K-1)*NJ*NL*NI
C    &              + (I-1)*NJ*NL
C    &              + (L-1)*NJ
C             DO 600 J = 1, NJ
C               JL = JL0 + J
*. Offsets for C(1,JB,j)
C               JBOFF(JL) = (J-1)*NJB + JB
*. integral * signs in SCR(jl,ik)
*. Integrals are stored as (j l i k )
C               SCR((IKEFF-1)*NJ*LL+JL) = FACTOR*XIJKL(JLIK)
C               SCR(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
C 600         CONTINUE
C 800       CONTINUE
C 700     CONTINUE
C 900   CONTINUE
*
C       CALL GSAXPY_LUCIA(SKII,CKJJ,SCR,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)
C1000 CONTINUE
      END IF
*. End of IROUTE branchning
*
      RETURN
      END

*
