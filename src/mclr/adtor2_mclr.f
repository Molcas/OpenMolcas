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
      SUBROUTINE AdToR2_MCLR(RHO2,RHO2T,ITYPE,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
*
* Add contributions to two electron density matrix RHO2
* output density matrix is in the form Rho2(ij,kl),(ij).ge.(kl)
*
*
* Jeppe Olsen, Fall of 96
*
*
* Itype = 1 => alpha-alpha or beta-beta loop
*              input is in form Rho2t(ik,jl)
* Itype = 2 => alpha-beta loop
*              input is in form Rho2t(ij,kl)
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION RHO2T(*)
*. Input and output
      DIMENSION RHO2(*)
*
      I       = 0     ! dummy initialize
      J       = 0     ! dummy initialize
      K       = 0     ! dummy initialize
      L       = 0     ! dummy initialize
      NII     = 0     ! dummy initialize
      NJJ     = 0     ! dummy initialize
      NKK     = 0     ! dummy initialize
      NLL     = 0     ! dummy initialize
      IIOFF   = 0     ! dummy initialize
      JJOFF   = 0     ! dummy initialize
      KKOFF   = 0     ! dummy initialize
      LLOFF   = 0     ! dummy initialize
      SIGN    = 0.0D0 ! dummy initialize
      IACTIVE = 0     ! dummy initialize
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to ADTOR2 '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' NI NJ NK NL = ', NI,NJ,NK,NL
        WRITE(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
        WRITE(6,*) ' ITYPE = ',ITYPE
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Initial two body density matrix '
          CALL PRSYM(RHO2,NORB**2)
        END IF
        WRITE(6,*) ' RHO2T : '
*       IF(ITYPE.EQ.1) THEN
*         IF(IOFF.EQ.KOFF) THEN
*           NROW = NI*(NI+1)/2
*         ELSE
*           NROW = NI*NK
*         END IF
*         IF(JOFF.EQ.LOFF) THEN
*           NCOL = NJ*(NJ+1)/2
*         ELSE
*           NCOL = NJ*NL
*         END IF
*       ELSE IF (ITYPE.EQ.2) THEN
*         NROW = NI*NJ
*         NCOL = NK*NL
*       END IF
*       CALL WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
      END IF
*
      IF(ITYPE.EQ.1) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
*. Four permutations
      DO IPERM = 1, 4
        IF(IPERM.EQ.1) THEN
          NII = NI
          IIOFF = IOFF
          NJJ = NJ
          JJOFF = JOFF
          NKK = NK
          KKOFF = KOFF
          NLL = NL
          LLOFF = LOFF
          SIGN = 1.0D0
          IACTIVE = 1
        ELSE IF(IPERM.EQ.2) THEN
          IF(IOFF.NE.KOFF) THEN
            NII = NK
            IIOFF = KOFF
            NKK = NI
            KKOFF = IOFF
            NJJ = NJ
            JJOFF = JOFF
            NLL = NL
            LLOFF = LOFF
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
          SIGN = -1.0D0
        ELSE IF(IPERM.EQ.3) THEN
          IF(JOFF.NE.LOFF) THEN
            NII = NI
            IIOFF = IOFF
            NKK = NK
            KKOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = -1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        ELSE IF(IPERM.EQ.4) THEN
          IF(IOFF.NE.KOFF.AND.JOFF.NE.LOFF) THEN
            NKK = NI
            KKOFF = IOFF
            NII = NK
            IIOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = 1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        END IF
*
C       IJOFF = (JJOFF-1)*NORB+IIOFF
C       KLOFF = (LLOFF-1)*NORB+KKOFF
C       IF(IACTIVE.EQ.1.AND.IJOFF.GE.KLOFF) THEN
        IF(IACTIVE.EQ.1) THEN
            DO II = 1, NII
              DO JJ = 1, NJJ
                DO KK = 1, NKK
                  DO LL = 1, NLL
                    IJ = (JJ+JJOFF-2)*NORB + II+IIOFF - 1
                    KL = (LL+LLOFF-2)*NORB + KK+KKOFF - 1
                    IF(IJ.GE.KL) THEN
                      IJKL = IJ*(IJ-1)/2+KL
                      IF(IPERM.EQ.1) THEN
                        I = II
                        K = KK
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.2) THEN
                        I = KK
                        K = II
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.3) THEN
                        I = II
                        K = KK
                        J = LL
                        L = JJ
                      ELSE IF(IPERM.EQ.4) THEN
                        I = KK
                        K = II
                        J = LL
                        L = JJ
                      END IF
                      IF(IOFF.NE.KOFF) THEN
                        IKIND = (K-1)*NI+I
                        NIK = NI*NK
                        SIGNIK = 1.0D0
                      ELSE
                        IKIND = MAX(I,K)*(MAX(I,K)-1)/2+MIN(I,K)
                        NIK = NI*(NI+1)/2
                        IF(I.EQ.MAX(I,K)) THEN
                          SIGNIK = 1.0D0
                        ELSE
                          SIGNIK = -1.0D0
                        END IF
                      END IF
                      IF(JOFF.NE.LOFF) THEN
                        JLIND = (L-1)*NJ+J
                        SIGNJL = 1.0D0
                      ELSE
                        JLIND = MAX(J,L)*(MAX(J,L)-1)/2+MIN(J,L)
                        IF(J.EQ.MAX(J,L)) THEN
                          SIGNJL = 1.0D0
                        ELSE
                          SIGNJL = -1.0D0
                        END IF
                      END IF
                      IKJLT = (JLIND-1)*NIK+IKIND
                      RHO2(IJKL) = RHO2(IJKL)
     &                           - SIGN*SIGNJL*SIGNIK*
     &                             RHO2T(IKJLT)
*. The minus : Rho2t comes as <a+i a+k aj al>, but we want
* <a+ia+k al aj>
                    END IF
                  END DO
                END DO
              END DO
            END DO
*. End of active/inactive if
        END IF
*. End of loop over permutations
      END DO
      ELSE IF(ITYPE.EQ.2) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJ = (J+JOFF-2)*NORB + I+IOFF - 1
             KL = (L+LOFF-2)*NORB + K+KOFF - 1
             IF(IJ.EQ.KL) THEN
               FACTOR = 2.0D0
             ELSE
               FACTOR= 1.0D0
             END IF
             IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
             RHO2(IJKL) = RHO2(IJKL) + FACTOR*RHO2T(IJKLT)
            END DO
          END DO
        END DO
      END DO
*
      Else If (itype.eq.3) Then
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJ = (J+JOFF-2)*NORB + I+IOFF - 1
             KL = (L+LOFF-2)*NORB + K+KOFF - 1
*            IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
             IJKL=(IJ-1)*NORB**2+KL
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
             RHO2(IJKL) = RHO2(IJKL) + RHO2T(IJKLT)
            END DO
          END DO
        END DO
      END DO

      END IF
*
*     i=1
*     j=4
*     k=5
*     l=6
*     iA1=itri((i-1)*nOrb+j,(k-1)*nOrb+l)
*     iA2=itri((j-1)*nOrb+i,(k-1)*nOrb+l)
*     iA3=itri((i-1)*nOrb+j,(l-1)*nOrb+k)
*     iA4=itri((j-1)*nOrb+i,(l-1)*nOrb+k)
*     rfel=RHO2(ia1)+RHO2(ia2)+RHO2(ia3)+RHO2(ia4)
*     IF(rfel.NE.0.0D0) THEN
*         Write(*,*) 'STOP'
*     END IF
*     IF(NTEST.GE.1000) THEN
*        WRITE(6,*) ' Updated two-body density matrix '
*        CALL PRSYM(RHO2,NORB**2)
*     END IF
*
      RETURN
      END
