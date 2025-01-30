!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994,1997, Jeppe Olsen                                 *
!***********************************************************************
      SUBROUTINE SKICKJ_LUCIA(   SKII,   CKJJ,    NKA,    NKB,  XIJKL,  &
     &                             NI,     NJ,     NK,     NL,   MAXK,  &
     &                           KBIB,  XKBIB,   KBJB,  XKBJB,  IKORD,  &
     &                           FACS, IROUTE)
!
!
! Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
!          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
!
!
!
! Jeppe Olsen, Spring of 94
!
! : Note : Route 1 has retired, March 97
      use lucia_data, only: MXPTSOB
      IMPLICIT None
      INTEGER NKA,NKB,NI,NJ,NK,NL,MAXK,IKORD,IROUTE
      REAL*8 FACS
!. Input
      REAL*8 CKJJ(*)
      REAL*8 XIJKL(*)
      INTEGER KBIB(MAXK,*),KBJB(MAXK,*)
      REAL*8 XKBIB(MAXK,*),XKBJB(MAXK,*)
!. Input and output
      REAL*8 SKII(*)
!. Scratch
      REAL*8 XIJILS(MXPTSOB)
      INTEGER JKINTOF,IKINTOF,NTEST,KB,LL,KK,L,K,IB,JB,ISOFF,ICOFF,     &
     &        INTOF,IMAX,J,I
      REAL*8 SGNK,SGNL,FACTOR
!. To get rid of annoying and incorrect compiler warnings
      JKINTOF = 0
      IKINTOF = 0
!
      IF(NI.GT.MXPTSOB.OR.NJ.GT.MXPTSOB.OR.NK.GT.MXPTSOB                &
     &   .OR.NL.GT.MXPTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : > MXPTSOB '
         WRITE(6,*) ' N, MXPTSOB ',MAX(NI,NJ,NK,NL),MXPTSOB
!         STOP ' Redim MXPTSOB '
        CALL SYSABENDMSG('lucia_util/skickj','Redim MXPTSOB',' ')
      END IF
!
      NTEST =  000
!
      IF(IROUTE.EQ.3) THEN
! S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
        DO KB = 1, NKB
!. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
!
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(NTEST.GE.100)                                       &
     &           WRITE(6,*) ' KB,K,L,IB,JB',KB,K,L,IB,JB
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
!. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
!
                   IF(IKORD.NE.0) THEN
!. Restrict so (ij) .le. (kl)
                     IMAX  = K
                     JKINTOF = INTOF + (K-1)*NJ
!                    CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                     DO J = L,NL
                       XIJILS(J) = XIJKL(JKINTOF-1+J)
                     END DO
                     XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
                     DO J = L+1, NL
                      XIJKL(JKINTOF-1+J) = 0.0D0
                     END DO
                   END IF
!                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF),                             &
     &                         CKJJ(ICOFF),XIJKL(INTOF),NKA,IMAX,  NKA, &
     &                              NJ,     NJ,   IMAX,   FACS,FACTOR , &
     &                               0)
                   IF(IKORD.NE.0) THEN
                      DO J = L,NL
                        XIJKL(JKINTOF-1+J) =  XIJILS(J)
                      END DO
!                    CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                   END IF
!
                 END IF
               END DO
!
             END IF
           END DO
         END IF
       END DO
!. (end over loop over Kb strings )
      ELSE IF(IROUTE.EQ.2) THEN
! S(I,Ka,Ib) = S(I,Ka,Ib) + sum(j) (ij!kl) C(j,Ka,Jb)
        DO KB = 1, NKB
!. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
!
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
!. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
!
                   IF(IKORD.NE.0) THEN
!. Restrict so (ji) .le. (kl)
                     IKINTOF = INTOF + (K-1)*NI
                     CALL COPVEC(XIJKL(IKINTOF),XIJILS,NI)
                     XIJKL(IKINTOF-1+L) = 0.5D0*XIJKL(IKINTOF-1+L)
                     DO I = L+1, NL
                      XIJKL(IKINTOF-1+I) = 0.0D0
                     END DO
                   END IF
!
!                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF),                             &
     &                         XIJKL(INTOF),CKJJ(ICOFF),NI,NKA, NI,     &
     &                              NJ,     NJ,    NKA,   FACS, FACTOR, &
     &                               0)
!
                 IF(IKORD.NE.0) THEN
                   CALL COPVEC(XIJILS,XIJKL(IKINTOF),NI)
                 END IF
!
                 END IF
               END DO
             END IF
           END DO
         END IF
       END DO
!. (end over loop over Kb strings )
!

      ELSE IF (IROUTE.EQ.1) THEN
         WRITE(6,*) ' Sorry route 1 has retired, March 1997'
!         STOP'SKICKJ:Invalid route=1'
        CALL SYSABENDMSG('lucia_util/skickj','Internal error',' ')
!     DO 1000 KB = 1, NKB
!. Number of nonvanishing a+lb !Kb>
!       LL = 0
!       DO L = 1, NL
!         IF(KBJB(KB,L).NE.0) LL = LL + 1
!       END DO
!
!       IKEFF = 0
!       DO 900 K = 1, NK
!         IB = KBIB(KB,K)
!         IF(IB.EQ.0) GOTO 900
!         SGNK = XKBIB(KB,K)
!
!         IF(IKORD.EQ.0) THEN
!            LI = NI
!            IMIN = 1
!         ELSE
!            LI = NI-K+1
!            IMIN = K
!         END IF
!
!         DO 700 I = IMIN, NI
!           IKEFF = IKEFF + 1
!           IOFF = (IKEFF-1)*NJ*LL
!. Offset for S(1,IB,i)
!           IBOFF(IKEFF)  = (I-1)*NIB+IB
!           LEFF = 0
!           DO 800 L = 1, NL
!             JB = KBJB(KB,L)
!             IF(JB.EQ.0) GOTO 800
!             LEFF = LEFF + 1
!             SGNL = XKBJB(KB,L)
!             IF(IKORD.EQ.1.AND.I.EQ.K)THEN
!                FACTOR = 0.5D0*SGNK*SGNL
!             ELSE
!                FACTOR =       SGNK*SGNL
!             END IF
!             JL0 = (LEFF-1)*NJ
!             JLIK0 = (K-1)*NJ*NL*NI
!    &              + (I-1)*NJ*NL
!    &              + (L-1)*NJ
!             DO 600 J = 1, NJ
!               JL = JL0 + J
!. Offsets for C(1,JB,j)
!               JBOFF(JL) = (J-1)*NJB + JB
!. integral * signs in SCR(jl,ik)
!. Integrals are stored as (j l i k )
!               SCR((IKEFF-1)*NJ*LL+JL) = FACTOR*XIJKL(JLIK)
!               SCR(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
! 600         CONTINUE
! 800       CONTINUE
! 700     CONTINUE
! 900   CONTINUE
!
!       CALL GSAXPY_LUCIA(SKII,CKJJ,SCR,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)
!1000 CONTINUE
      END IF
!. End of IROUTE branchning
!
      END SUBROUTINE SKICKJ_LUCIA

!
