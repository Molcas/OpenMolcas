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
* Copyright (C) 1996, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ABTOR2(    SKII,    CKJJ,     NKA,     NIB,     NJB,
     &                       NKB,   RHO2B,      NI,      NJ,      NK,
     &                        NL,    MAXK,    KBIB,   XKBIB,    KBJB,
     &                     XKBJB,   IKORD)
*
* Obtain contributions alpha-beta contributions to two-particle
* density matrix
*
* Rho2b(ij,kl)  = RHo2b(ij,kl)
*               + sum(Ka) Skii(Ka,i,Ib)<Ib!Eb(kl)!Jb> Ckjj(Ka,j,Jb)
*
*
* Jeppe Olsen, Fall of 96
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(*),SKII(*)
      DIMENSION RHO2B(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*
      IF(IKORD.NE.0) THEN
        WRITE(6,*) ' ABTOR2 : IKORD .NE. 0 '
        WRITE(6,*) ' I am not ready for this '
*        STOP     ' ABTOR2 : IKORD .NE. 0 '
        CALL SYSABENDMSG('lucia_util/abtor2_gas',
     &                          'Internal error',' ')
      END IF
*
*. Excitations <Ib!Eb(kl)!Jb>
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
                   KLOFF= ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
*
C                  IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
C                    IMAX  = K
C                    JKINTOF = INTOF + (K-1)*NJ
C                    DO J = L,NL
C                      XIJILS(J) = XIJKL(JKINTOF-1+J)
C                    END DO
C                    XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
C                    DO J = L+1, NL
C                     XIJKL(JKINTOF-1+J) = 0.0D0
C                    END DO
C                  END IF
                   ONE = 1.0D0
                   CALL MATML7(RHO2B(KLOFF),
     &                         SKII(ISOFF),CKJJ(ICOFF),NI,  NJ, NKA,
     &                            IMAX,    NKA,     NJ,    ONE,FACTOR ,
     &                               1)
C                  IF(IKORD.NE.0) THEN
C                     DO J = L,NL
C                       XIJKL(JKINTOF-1+J) =  XIJILS(J)
C                     END DO
C                  END IF
*
                 END IF
               END DO
*
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NIB)
        CALL Unused_integer(NJB)
      END IF
      END
