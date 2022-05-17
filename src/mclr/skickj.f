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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE SKICKJ(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,
     &                  NI,NJ,NK,NL,MAXK,
     &                  KBIB,XKBIB,KBJB,XKBJB,IKORD,
     &                  IXBOFF,JXBOFF,SXCR,IROUTE,NTEST )
*
*
* Calculate S(Ka,Ib,i) = S(Ka,Ib,i)
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
*
*
* Jeppe Olsen, Spring of 94
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "detdim.fh"
#include "stdalloc.fh"
*. Input
      DIMENSION CKJJ(NKA*NJ,*)
*. Note if Iroute = 2 the form is C(j,Ka,Jb)
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(NKA*NI,*)
*. Note if Iroute = 2 the form is S(i,Ka,Ib)
*. Scratch
      PARAMETER(MXTSOB=35)
      DIMENSION IBOFF(MXTSOB*MXTSOB),JBOFF(MXTSOB*MXTSOB)
      Real*8, Allocatable:: KSKICK(:)
*
      MAXORB = MAX(NI,NJ,NK,NL)
      LENGTH = MAXORB*MAXORB*MAXORB*MAXORB
      Call mma_allocate(KSKICK,LENGTH,Label='KSKICK')
*
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &   .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Write (6,*) ' Redim MXTSOB in SKICKJ'
         Call Abend()
      END IF
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
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
C                  ISOFF = (IB-1)*NI*NKA + 1
C                  ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
*
                   ONE = 1.0D0

                   CALL  DGEMM_('N','N',NKA,NI,NJ,
     &                         FACTOR ,CKJJ(1,jB),max(1,NKA),
     &                                 XIJKL(INTOF),max(1,NJ),
     &                         ONE,SKII(1,iB),max(1,NKA))
*
                 END IF
               END DO
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
C                  ISOFF = (IB-1)*NI*NKA + 1
C                  ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
*
                   ONE = 1.0D0
                   CALL  DGEMM_('N','N',NI,NKA,NJ,FACTOR ,
     &                         XIJKL(INTOF),max(1,NI),
     &                         CKJJ(1,JB),max(1,NJ),
     &                         ONE,SKII(1,IB),max(1,NI))
                 END IF
               END DO
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*

      ELSE IF (IROUTE.EQ.1) THEN



      DO 1000 KB = 1, NKB
*. Number of nonvanishing a+lb !Kb>
        LL = 0
        DO L = 1, NL
          IF(KBJB(KB,L).NE.0) LL = LL + 1
        END DO
*
        IKEFF = 0
        DO 900 K = 1, NK
          IB = KBIB(KB,K)
          IF(IB.EQ.0) GOTO 900
          SGNK = XKBIB(KB,K)
*
          IF(IKORD.EQ.0) THEN
             LI = NI
             IMIN = 1
          ELSE
             LI = NI-K+1
             IMIN = K
          END IF
*
          DO 700 I = IMIN, NI
            IKEFF = IKEFF + 1
            IOFF = (IKEFF-1)*NJ*LL
*. Offset for S(1,IB,i)
            IBOFF(IKEFF)  = (I-1)*NIB+IB
            LEFF = 0
            DO 800 L = 1, NL
              JB = KBJB(KB,L)
              IF(JB.EQ.0) GOTO 800
              LEFF = LEFF + 1
              SGNL = XKBJB(KB,L)
              IF(IKORD.EQ.1.AND.I.EQ.K)THEN
                 FACTOR = 0.5D0*SGNK*SGNL
              ELSE
                 FACTOR =       SGNK*SGNL
              END IF
              JL0 = (LEFF-1)*NJ
              JLIK0 = (K-1)*NJ*NL*NI
     &              + (I-1)*NJ*NL
     &              + (L-1)*NJ

              DO 600 J = 1, NJ
                JL = JL0 + J
*. Offsets for C(1,JB,j)
                JBOFF(JL) = (J-1)*NJB + JB
*. integral * signs in SCR(jl,ik)
*. Integrals are stored as (j l i k )
                KSKICK(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
  600         CONTINUE

  800       CONTINUE

  700     CONTINUE

  900   CONTINUE
*
        CALL GSAXPY(SKII,CKJJ,KSKICK,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)

 1000 CONTINUE
      END IF
*. End of IROUTE branching
*
      Call mma_deallocate(KSKICK)
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(IXBOFF)
        CALL Unused_integer(JXBOFF)
        CALL Unused_real(SXCR)
        CALL Unused_integer(NTEST)
      END IF

      END
