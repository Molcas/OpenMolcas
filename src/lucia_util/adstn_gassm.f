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
      SUBROUTINE ADSTN_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,NACGSOB,
     *              ISTMAP,SGNMAP,SCLFAC)
*
* Creation mappings from K-strings of given sym in each gasspace
*
* Input
* NSTB : Number of strings before active gasspace
* NSTA : Number of strings after accive gasspace
* IOFFK : Offset for K group of strings in active gasspace, i.e. start of
*         this symmetry of active K group strings
* IOFFI : Offset for I group of strings in active gasspace, i.e. start of
*         this symmetry of active I group strings
* IOFFISP: Offset for this symmetrydistribution of active I supergroup strings
* IOFFKSP: Offset for this symmetrydistribution of active K supergroup strings
* ICREORB : Orbital part of creation map for active K groupstrings
* ICRESTR : String  part of creation map for active K groupstrings
* IORBTSF   : First active orbital ( first orbital in in active GASspace
*           with required sym)
* IORBTF   : First orbital in active gas space, (can have any sym)
* NORBTS  : Number of orbitals of given symmetry and type
* NSTAK : Number of K groupstrings with given correct symmetry
* NSTAKT: Total Number of K groupstrings in active group (all symmetries)
* NSTAKTS: Total Number of K supergroup strings with correct symmetry
* ISTAKTS: Offset for K supergroup strings with hiven symmetrydistribution
* NSTAI : Number of I groupstrings in active gasspace
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION ICREORB(NACGSOB,*), ICRESTR(NACGSOB,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
*
C?    WRITE(6,*) ' ADSTN_GASSM : NSTA, NSTB, NSTAK',NSTA,NSTB,NSTAK
C?    WRITE(6,*) ' IOFFISP,IOFFKSP', IOFFISP, IOFFKSP
C?    WRITE(6,*) ' IORBTSF IORBTF ', IORBTSF,IORBTF
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
C?    WRITE(6,*) ' NSTAKT ', NSTAKT
*
*      SIGN0 = DBLE((-1)**NELB)*SCLFAC
      IF( MOD(NELB,2).EQ.0 ) THEN
        SIGN0 = SCLFAC
      ELSE
        SIGN0 = -SCLFAC
      END IF
C?    WRITE(6,*) ' NELB sign0 = ', NELB, SIGN0
      DO KSTR = IOFFK, NSTAK+IOFFK-1
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
C?         write(6,*) 'IORB IORBRT KSTR ', IORB,IORBRT, KSTR
C?         WRITE(6,*) 'ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)',
C?   &                 ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)
          IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Excitation is open, corresponding active I string
            IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
              SIGN = SIGN0
              ISTR = ICRESTR(IORBRT,KSTR)
            ELSE
              SIGN = -SIGN0
              ISTR = -ICRESTR(IORBRT,KSTR)
            END IF
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          write(6,*) ' ISTR IADRK0 IADRI0 = ', ISTR, IADRK0,IADRI0
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
            DO IB = 1, NSTB
              DO IA = 1, NSTA
C               IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
C         ELSE
C            SIGN = 0.0D0
C            ISTR = 0
*. This Creation is inactive for all choices of strings in supergroup
*. before and after the active type.
C           IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
C?          write(6,*) ' ISTR IADRK0 = ', ISTR, IADRK0
*
C           DO IB = 1, NSTB
C             DO IA = 1, NSTA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
C               IF(ISTMAP(KBKA,IORBRTS).NE.999) THEN
C                 WRITE(6,*) ' overwriting ??? '
C                 WRITE(6,*) ' Element ', (IORBRTS-1)*NSTAKTS+KBKA
C                 STOP
C               END IF
C               ISTMAP(KBKA,IORBRTS) = ISTR
C               SGNMAP(KBKA,IORBRTS) = SIGN
C
C             END DO
C           END DO
          END IF
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
COLD        IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
COLD        IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          write(6,*) ' ISTR IADRK0 IADRI0 = ', ISTR, IADRK0,IADRI0
*
COLD        DO IB = 1, NSTB
COLD          DO IA = 1, NSTA
COLD            IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
COLD            KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
COLD            ISTMAP(KBKA,IORBRTS) = IBKA
COLD            SGNMAP(KBKA,IORBRTS) = SIGN
COLD          END DO
COLD        END DO
C         END IF
*
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADSTN_GASSM '
        WRITE(6,*) ' ======================== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
            CALL WRTMAT(SGNMAP(1,IORBR),1,NK,1,NK)
          END DO
        END IF
      END IF

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NSTAKT)
        CALL Unused_integer(ISTAKTS)
      END IF
      END
