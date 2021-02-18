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
*               2003, Jesper Wisborg Krogh                             *
************************************************************************
      SUBROUTINE ADTOR2(    RHO2,   RHO2S,   RHO2A,   RHO2T,   ITYPE,
     &                        NI,    IOFF,      NJ,    JOFF,      NK,
     &                      KOFF,      NL,    LOFF,    NORB,   IPACK)
*
* Add contributions to two electron density matrix RHO2
* output density matrix is in the form Rho2(ij,kl),(ij).ge.(kl)
*
*
* Jeppe Olsen, Fall of 96
* Jesper Wisborg Krogh, September 03: Can now symmetry pack on the fly
*
*
* Itype = 1 => alpha-alpha or beta-beta loop
*              input is in form Rho2t(ik,jl)
* Itype = 2 => alpha-beta loop
*              input is in form Rho2t(ij,kl)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "crun.fh"
#include "orbinp.fh"
*.Input
      DIMENSION RHO2T(*)
      LOGICAL   IPACK
*. Input and output
      DIMENSION RHO2(*),RHO2S(*),RHO2A(*)
* Program flow flags
      LOGICAL   IPROCEED,DO_IJKL,DO_IJKL_PACK,DO_JIKL_PACK
*
* Some dummy initializations
      NII = 0 ! jwk-cleanup
      NJJ = 0 ! jwk-cleanup
      NKK = 0 ! jwk-cleanup
      NLL = 0 ! jwk-cleanup
      KKOFF = 0 ! jwk-cleanup
      LLOFF = 0 ! jwk-cleanup
      SIGN = 1.0D0 ! jwk-cleanup
      IACTIVE = 1 ! jwk-cleanup
      I = 0 ! jwk-cleanup
      K = 0 ! jwk-cleanup
      J = 0 ! jwk-cleanup
      L = 0 ! jwk-cleanup
      I_PACK    = 0
      J_PACK    = 0
      K_PACK    = 0
      L_PACK    = 0
      IJ_PACK   = 0
      JI_PACK   = 0
      KL_PACK   = 0
      FACTOR    = 0.0d0
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to ADTOR2 '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' NI NJ NK NL = ', NI,NJ,NK,NL
        WRITE(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
        WRITE(6,*) ' ITYPE = ',ITYPE
        WRITE(6,*) ' NORB ', NORB
        IF(NTEST.GE.2000 .AND. .NOT. IPACK) THEN
          WRITE(6,*) ' Initial two body density matrix '
          CALL PRSYM(RHO2,NORB**2)
        END IF
*
        WRITE(6,*) ' RHO2T : '
        IF(ITYPE.EQ.1) THEN
          IF(IOFF.EQ.KOFF) THEN
            NROW = NI*(NI+1)/2
          ELSE
            NROW = NI*NK
          END IF
          IF(JOFF.EQ.LOFF) THEN
            NCOL = NJ*(NJ+1)/2
          ELSE
            NCOL = NJ*NL
          END IF
        ELSE IF (ITYPE.EQ.2) THEN
          NROW = NI*NJ
          NCOL = NK*NL
        END IF
        CALL WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
      END IF

C?    RETURN
*
      NELMNT = NORB**2*(NORB**2+1)/2
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
C         IJOFF = (JJOFF-1)*NORB+IIOFF
C         KLOFF = (LLOFF-1)*NORB+LLOFF
            DO II = 1, NII
              DO JJ = 1, NJJ
                DO KK = 1, NKK
                  DO LL = 1, NLL
                    IJ = (JJ+JJOFF-2)*NORB + II+IIOFF - 1
                    KL = (LL+LLOFF-2)*NORB + KK+KKOFF - 1

                    IPROCEED     = .FALSE.
                    DO_IJKL      = .FALSE.
                    DO_IJKL_PACK = .FALSE.
                    DO_JIKL_PACK = .FALSE.
                    IF (IPACK) THEN
                       I_PACK    = II + IIOFF - 1
                       J_PACK    = JJ + JJOFF - 1
                       K_PACK    = KK + KKOFF - 1
                       L_PACK    = LL + LLOFF - 1

                       IJ_PACK   = I_PACK*(I_PACK-1)/2 + J_PACK
                       JI_PACK   = J_PACK*(J_PACK-1)/2 + I_PACK
                       KL_PACK   = K_PACK*(K_PACK-1)/2 + L_PACK

                       IF (K_PACK .EQ. L_PACK) THEN
                          FACTOR = 0.25d0
                       ELSE
                          FACTOR = 0.5d0
                       END IF

                       IF (K_PACK .GE. L_PACK) THEN
                          IF (I_PACK .GE. J_PACK .AND.
     &                          IJ_PACK .GE. KL_PACK) THEN
                             DO_IJKL_PACK = .TRUE.
                             IPROCEED     = .TRUE.
                          END IF
                          IF (J_PACK .GE. I_PACK .AND.
     &                          JI_PACK .GE. KL_PACK) THEN
                             DO_JIKL_PACK = .TRUE.
                             IPROCEED     = .TRUE.
                          END IF
                       END IF
                    ELSE IF (IJ .GE. KL) THEN
                       DO_IJKL  = .TRUE.
                       IPROCEED = .TRUE.
                       FACTOR   = 1.0D0
                    END IF

                    IF (IPROCEED) THEN
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
                      TERM = FACTOR*SIGN*SIGNJL*SIGNIK*RHO2T(IKJLT)

                      IF (DO_IJKL_PACK) THEN
                         IJKL_PACK = IJ_PACK*(IJ_PACK-1)/2 + KL_PACK
                         RHO2S(IJKL_PACK) = RHO2S(IJKL_PACK) - TERM
                         RHO2A(IJKL_PACK) = RHO2A(IJKL_PACK) - TERM
                      END IF
                      IF (DO_JIKL_PACK) THEN
                         JIKL_PACK = JI_PACK*(JI_PACK-1)/2 + KL_PACK
                         RHO2S(JIKL_PACK) = RHO2S(JIKL_PACK) - TERM
                         RHO2A(JIKL_PACK) = RHO2A(JIKL_PACK) + TERM
                      END IF
                      IF (DO_IJKL) THEN
                         IJKL = IJ*(IJ-1)/2+KL
                         IF(IJKL.GT.NELMNT) THEN
                            WRITE(6,*) ' Problemo 1 : IJKL .gt. NELMNT'
                            WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                            WRITE(6,*) ' IJ, KL', IJ,KL
                            WRITE(6,*) ' JJ JJOFF ', JJ,JJOFF
                            WRITE(6,*) ' II IIOFF ', II,IIOFF
                            WRITE(6,*) ' IPERM = ', IPERM
                         END IF
                         RHO2(IJKL) = RHO2(IJKL) - TERM
                      END IF
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
               FACTOR = 1.0D0
             END IF
             IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
                      IF(IJKL.GT.NELMNT) THEN
                         WRITE(6,*) ' Problemo 2 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                      END IF
             TERM = FACTOR*RHO2T(IJKLT)
             IF (IPACK) THEN
                DO IPERM = 1,2
                   IF (IPERM .EQ. 1) THEN
                      I_PACK    = I + IOFF - 1
                      J_PACK    = J + JOFF - 1
                      K_PACK    = K + KOFF - 1
                      L_PACK    = L + LOFF - 1
                      IJ_PACK   = I_PACK*(I_PACK-1)/2 + J_PACK
                      JI_PACK   = J_PACK*(J_PACK-1)/2 + I_PACK
                      KL_PACK   = K_PACK*(K_PACK-1)/2 + L_PACK
                   ELSE
                      I_PACK    = K + KOFF - 1
                      J_PACK    = L + LOFF - 1
                      K_PACK    = I + IOFF - 1
                      L_PACK    = J + JOFF - 1
                      IJ_PACK   = I_PACK*(I_PACK-1)/2 + J_PACK
                      JI_PACK   = J_PACK*(J_PACK-1)/2 + I_PACK
                      KL_PACK   = K_PACK*(K_PACK-1)/2 + L_PACK
                   END IF


                IF (K_PACK .EQ. L_PACK) THEN
                   FACTOR_PACK = .25d0
                ELSE
                   FACTOR_PACK = .5d0
                END IF
                IF (I_PACK .EQ. K_PACK .AND. J_PACK .EQ. L_PACK) THEN
                   FACTOR_PACK = FACTOR_PACK * .5D0
                END IF

                IF (I_PACK .GE. J_PACK .AND. K_PACK .GE. L_PACK
     &                        .AND. IJ_PACK .GE. KL_PACK) THEN
                   IJKL_PACK = IJ_PACK*(IJ_PACK-1)/2 + KL_PACK
                   RHO2S(IJKL_PACK) = RHO2S(IJKL_PACK)
     &                              + FACTOR_PACK*TERM
                   RHO2A(IJKL_PACK) = RHO2A(IJKL_PACK)
     &                              + FACTOR_PACK*TERM
                END IF
                IF (J_PACK .GE. I_PACK .AND. K_PACK .GE. L_PACK
     &                       .AND. JI_PACK .GE. KL_PACK) THEN
                   JIKL_PACK = JI_PACK*(JI_PACK-1)/2 + KL_PACK
                   IJKL_PACK = IJ_PACK*(IJ_PACK-1)/2 + KL_PACK
                   RHO2S(JIKL_PACK) = RHO2S(JIKL_PACK)
     &                              + FACTOR_PACK*TERM
                   RHO2A(JIKL_PACK) = RHO2A(JIKL_PACK)
     &                              - FACTOR_PACK*TERM
                END IF
                END DO
             ELSE
                RHO2(IJKL) = RHO2(IJKL) + TERM
             END IF
            END DO
          END DO
        END DO
      END DO
*
      END IF
*
      RETURN
      END

