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
      SUBROUTINE MXRESCPH(     IAB,  IOCTPA,  IOCTPB,  NOCTPA,  NOCTPB,
     &                       NSMST,NSTFSMSPGP,MXPNSMST, NSMOB,MXPTOB,
     &                       NTPOB,   NTSOB,  NTESTG,   MXPKA, NEL1234,
     &                        MXCJ,  MXCIJA,  MXCIJB, MXCIJAB,  MXSXBL,
     &                    MXADKBLK,  IPHGAS,NHLFSPGP,    MNHL, IADVICE,
*
     &                    MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
* Find largest dimension of matrix S(P,Ia,I,K) for a single K-string
*
* Particle hole version : hole electrons added, particle elec removed
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups
*
* Version of Jan 98 : IPHGAS added

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
      DIMENSION IPHGAS(*)
      INTEGER NHLFSPGP(*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
      MXKAB = 0
      MXCJ = 0
      MXCJ_ALLSYM = 0
      MXADKBLK = 0
      MXADKBLK_AS = 0
      MX_NSPII = 0
      DO IAORC= 1,2
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            ITOTB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
              ITOTB = ITOTB + NSTFSMSPGP(ISM,IBTPABS)
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB,ITOTB = ', MXB,ITOTB
            DO 300 IOBTP = 1, NTPOB
*. No K strings obtained from creation in particle space
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 300
*. type of K string obtained
              CALL NEWTYP(IATPABS,IAORC,IOBTP,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
*. addi constraint to avoid calc with long columns and few rows
*. Works only in connection with active advice routine !
              IF(KATP.GT.0) THEN
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IBTPABS)+NHLFSPGP(KATP).LT.MNHL.AND.
     &          NHLFSPGP(IATPABS).GT.(NHLFSPGP(IBTPABS)+1)) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IBTPABS,KATP',
C    &            IOBTP,IBTPABS,KATP
                  KATP = 0
                END IF
              END IF
*
              IF(KATP.GT.0) THEN
C                    DIM_SPII(IASPGRP,IBSPGRP,IOBTP,IAB,IAC,NSPII)
C               CALL DIM_SPII(IATPABS,IBTPABS,IOBTP,1,IAORC,NSPII)
C               MX_NSPII = MAX(MX_NSPII,NSPII)
                MX_NSPII = 0
              END IF

              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                NSOB_AS = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                  NSOB_AS = NSOB_AS + NTSOB(IOBTP,ISMOB)
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                MXADKBLK_AS = MAX(MXADKBLK_AS,NSOB_AS*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                LCJBLK_ALLSYM = NSOB_AS*MXKA*ITOTB

                IF(LCJBLK.GT.MXCJ) THEN
                  MXCJ = LCJBLK
                  IATP_MX = IATP
                  IBTP_MX = IBTP
                  KATP_MX = KATP
                  IOBTP_MX = IOBTP
                END IF
                MXCJ_ALLSYM = MAX(MXCJ_ALLSYM,LCJBLK_ALLSYM)
                MXKAB = MAX(MXKAB,MXKA)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
      END DO
*     ^ End of anni/crea map
*
* matrix C(j,Ia,Kb)
*
      DO IAORC = 1, 2
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            ITOTA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
              ITOTA = ITOTA + NSTFSMSPGP(ISM,IATPABS)
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 2812
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IATPABS)+NHLFSPGP(KBTP).LT.MNHL.AND.
     &          NHLFSPGP(IBTPABS).GT.NHLFSPGP(IATPABS)+1) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IATPABS,KBTP',
C    &            IOBTP,IATPABS,KBTP
                  KBTP = 0
                END IF
              END IF
              IF(KBTP.GT.0) THEN
C               CALL DIM_SPII(IATPABS,IBTPABS,IOBTP,2,IAORC,NSPII)
C               MX_NSPII = MAX(MX_NSPII,NSPII)
                MX_NSPII = 0
              END IF
*
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                NSOB_AS = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                  NSOB_AS = NSOB_AS + NTSOB(IOBTP,ISMOB)
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
                MXADKBLK_AS = MAX(MXADKBLK_AS,NSOB_AS*MXKBO)
CJULY29         LCJBLK = MXSOB*MXKB*MXB
                LCJBLK = MXSOB*MXKB*MXA
                LCJBLK_ALLSYM = NSOB_AS*MXKB*ITOTA
                MXCJ = MAX(MXCJ,LCJBLK)
                MXCJ_ALLSYM = MAX(MXCJ_ALLSYM,LCJBLK_ALLSYM)
                MXKAB = MAX(MXKAB,MXKB)
*
              END IF
 2812       CONTINUE
            END DO
          END IF
        END DO
      END DO
      END DO
*     ^ End of loop over creation/annihilation
      IF(NTEST.GT.100) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
        WRITE(6,*) ' MXCJ_ALLSYM = ', MXCJ_ALLSYM
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,IAORC,IOBTP,K1ATP)
*. No N+1 mappings for particle spaces
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) K1ATP = 0
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,JAORC,JOBTP,KATP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KATP = 0
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP, JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP, IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,kb)
* both Ka and Ib blocked
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1)K1BTP = 0
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,JAORC,JOBTP,KBTP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KBTP = 0
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP,JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP,IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,kb)
* both Ka and Kb blocked
*
*. Modified : Only used if atmost two elecs in i and j
*.            No batching
*.            Used for hardwired few electron code
      MXCIJAB = 0
      MXKACTEL = 1
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
C           NEL1234(JOBTP,IATPABS)
              IF(KATP.GT.0) THEN
                IF(NEL1234(IOBTP,KATP).GT.MXKACTEL) KATP = 0
              END IF
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
*. No partitioning
C               IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA

                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,KBTP)
                  IF(KBTP.GT.0) THEN
                    IF(NEL1234(JOBTP,KBTP).GT.MXKACTEL) KBTP = 0
                  END IF
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
*. No partitioning
C                   IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IATPABS,1,JOBTP,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KATP,2,IOBTP,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IBTPABS,1,JOBTP,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KBTP,2,IOBTP,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
        WRITE(6,*) ' MXRESC_PH : MXKAB = ', MXKAB
        WRITE(6,*) ' Info on largest C(Ka,j,Jb) block'
        WRITE(6,*) ' IATP_MX, IBTP_MX, KATP_MX, IOBTP_MX ',
     &               IATP_MX, IBTP_MX, KATP_MX, IOBTP_MX
        WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
*
      RETURN
      END
