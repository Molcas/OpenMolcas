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
      SUBROUTINE MXRESC(IAB,IACLS,IBCLS,NOCTPA,NOCTPB,NSMST,
     &                  NSSOA,NSSOB,KCLS,KSSOA,KOCTPA,
     &                  KSSOB,KOCTPB,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  K2SSOA,K2OCTPA,K2SSOB,K2OCTPB,NAEL123,
     &                  NBEL123,MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXIJST,MXIJSTF)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
* Find largest dimension of a+ia+j!k> MXIJST, MXIJSTF
*
* Largest block of single excitations MXSXBL
*
*     To mess it up for the enemy: kcls
*
*     Input:
*            iab
*            iacls,ibcls
*            noctpa,noctpb
*            nsmst
*            nssoa,nssob
*            kssoa
*            koctps,koctpb
*            nsmob
*            mxtpob,ntpob,ntsob
*            ntestg,mxpka
*            k2ssoa,k2ssob
*            k2octpb
*            nael123,nbel123
*     Output
*            mxcj,mxcija,mxcijb,mxcijab,mxsxbl,MXIJST,MXIJSTF
*

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSSOA(NOCTPA,NSMST),NSSOB(NOCTPB,NSMST)
      DIMENSION KSSOA(KOCTPA,NSMST),KSSOB(KOCTPB,NSMST)
      DIMENSION K2SSOA(K2OCTPA,NSMST)
      DIMENSION K2SSOB(K2OCTPB,NSMST)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NAEL123(NTPOB,*),NBEL123(NTPOB,*)
*
*
* matrix C(j,Ka,Ib)
*
      MXCJ = 0
      DO 100 IATP = 1, NOCTPA
        DO 200 IBTP = 1, NOCTPB
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB = MAX(MXB,NSSOB(IBTP,ISM))
  210       CONTINUE
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,KSSOA(KATP,KSM))
  310           CONTINUE
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
*
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      DO 101 IATP = 1, NOCTPA
        DO 201 IBTP = 1, NOCTPB
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXA = 0
            DO 211 ISM = 1, NSMST
              MXA = MAX(MXA,NSSOA(IATP,ISM))
  211       CONTINUE
            DO 301 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,KBCLS,KBTP)
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO 311 KSM = 1, NSMST
                  MXKB = MAX(MXKB,KSSOB(KBTP,KSM))
  311           CONTINUE
                IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                DO 321 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  321           CONTINUE
*
                LCJBLK = MXSOB*MXKB*MXA
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  301       CONTINUE
          END IF
  201   CONTINUE
  101 CONTINUE
*
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      MXIJSTA = 0
      MXIJSTAF = 0
      DO  IATP = 1, NOCTPA
        DO  IBTP = 1, NOCTPB

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSSOB(IBTP,ISM))
            END DO
            IF(MXIB.GT.MXPKA.AND.MXPKA.GT.0) MXIB = MXPKA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,K1ACLS,K1ATP)
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP_MCLR(K1ACLS,K1ATP,[1],[JOBTP],1,KACLS,
     &                             KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,K2SSOA(KATP,KSM))
                    END DO
                    MXKAF = MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
*
                    MXIJSTA = MAX(MXIJSTA,MXISOB*MXJSOB*MXKA)
                    MXIJSTAF= MAX(MXIJSTAF,MXISOB*MXJSOB*MXKAF)
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
*
* matrix C(ij,Ia,kb)
* both Ka and Ib blocked
*
      MXCIJB = 0
      MXIJSTB = 0
      MXIJSTBF = 0
      MXIA     = 0 ! dummy initialize
      DO  IATP = 1, NOCTPA
        DO  IBTP = 1, NOCTPB
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSSOA(IATP,ISM))
            END DO
            IF(MXIA.GT.MXPKA.AND.MXPKA.GT.0) MXIA = MXPKA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP_MCLR(IBCLS,IBTP,[1],[IOBTP],1,K1BCLS,K1BTP)
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP_MCLR(K1BCLS,K1BTP,[1],[JOBTP],1,KBCLS,
     &                             KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,K2SSOB(KBTP,KSM))
                    END DO
                    MXKBF = MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
*
                    MXIJSTB = MAX(MXIJSTB,MXISOB*MXJSOB*MXKB)
                    MXIJSTBF= MAX(MXIJSTBF,MXISOB*MXJSOB*MXKBF)
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
*
* matrix C(ij,Ka,kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      DO  IATP = 1, NOCTPA
        DO  IBTP = 1, NOCTPB

          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP_MCLR(IACLS,IATP,[1],[IOBTP],1,KACLS,KATP)
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,KSSOA(KATP,KSM))
                END DO
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,KSSOB(KBTP,KSM))
                    END DO
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
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
      MXSXBL = 0
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSSOA(IATP,ISM))
        END DO
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_MCLR(IACLS,IATP,[1],[JOBTP],1,KACLS,KATP)
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NAEL123(JOBTP,IATP))
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_MCLR(KACLS,KATP,[2],[IOBTP],1,JACLS,JATP)
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
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSSOB(IBTP,ISM))
        END DO
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_MCLR(IBCLS,IBTP,[1],[JOBTP],1,KBCLS,KBTP)
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NBEL123(JOBTP,IBTP))
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_MCLR(KBCLS,KBTP,[2],[IOBTP],1,JBCLS,JBTP)
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
      MXIJST = MAX(MXIJSTA,MXIJSTB)
      MXIJSTF = MAX(MXIJSTAF,MXIJSTBF)
*
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(KCLS)
        CALL Unused_integer(NTESTG)
      END IF
      END
