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
      SUBROUTINE ADAST_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,NACGSOB,
     *              ISTMAP,SGNMAP,SCLFAC,IAC,LROW_IN,IEC)
*
* Annihilation/Creation mappings from K-strings of given sym in each gasspace
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
      DIMENSION ICREORB(LROW_IN,*), ICRESTR(LROW_IN,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
* Some dummy initializations
      SIGN = 0.0D0 ! jwk-cleanup
      ISTR = 0 ! jwk-cleanup
*
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
*
C?    WRITE(6,*) ' ICRESTR '
C?    CALL IWRTMA(ICRESTR,LROW_IN,NSTAK,LROW_IN,NSTAK)
C?    WRITE(6,*) ' IOFFI = ', IOFFI
*PAM2009      SIGN0 = DBLE((-1)**NELB) * SCLFAC
      IF(MOD(NELB,2).EQ.0) THEN
        SIGN0=SCLFAC
      ELSE
        SIGN0=-SCLFAC
      END IF
      DO KSTR = IOFFK, NSTAK+IOFFK-1
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
C?        write(6,*) ' IORBRTS IORBRT', IORBRTS,IORBRT
*. Change of active group
          I_AM_ACTIVE = 0
          IF(IAC.EQ.2) THEN
C?          WRITE(6,*) ' ICREORB = ', ICREORB(IORBRT,KSTR)
C?          WRITE(6,*) ' ICRESTR = ', ICRESTR(IORBRT,KSTR)
            IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Creation is nonvanishing
              I_AM_ACTIVE = 1
              IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                SIGN = SIGN0
                ISTR = ICRESTR(IORBRT,KSTR)
              ELSE
                SIGN = -SIGN0
                ISTR = -ICRESTR(IORBRT,KSTR)
              END IF
            END IF
          ELSE IF(IAC.EQ.1) THEN
             IF(IEC.EQ.1) THEN
*. Expanded map
               IF(ICREORB(IORBRT,KSTR) .LT. 0 ) THEN
*. Annihilation is non-vanishing
                 I_AM_ACTIVE = 1
                 IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                   SIGN = SIGN0
                   ISTR = ICRESTR(IORBRT,KSTR)
                 ELSE
                   SIGN = -SIGN0
                   ISTR = -ICRESTR(IORBRT,KSTR)
                 END IF
               END IF
             ELSE
*. Compressed map
               DO IROW = 1, LROW_IN
C?              write(6,*) ' IROW, ICREORB(IROW,KSTR)' ,
C?   &                       IROW, ICREORB(IROW,KSTR)
COLD             IF(ICREORB(IROW,KSTR) .EQ. -IORBRT ) THEN
                 IF(ICREORB(IROW,KSTR) .EQ. -IORB   ) THEN
*. Annihilation is non-vanishing
                   I_AM_ACTIVE = 1
                   IF(ICRESTR(IROW,KSTR) .GT. 0 ) THEN
                     SIGN = SIGN0
                     ISTR = ICRESTR(IROW,KSTR)
                   ELSE
                     SIGN = -SIGN0
                     ISTR = -ICRESTR(IROW,KSTR)
                   END IF
                 END IF
               END DO
             END IF
*            ^ End of expanded/compact switch
           END IF
*          ^ End of Creation/annihilation switch

          IF(I_AM_ACTIVE .EQ. 1  ) THEN
*. Excitation is open, corresponding active I string
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
C?          WRITE(6,*) ' ISTR, relative = ', ISTR
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          WRITE(6,*) ' IADRK0 IOFFK IOFFKSP ',
C?   &                   IADRK0,IOFFK,IOFFKSP
C?          WRITE(6,*) ' IADRI0, IOFFISP ',
C?   &                   IADRI0, IOFFISP
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
*
C?          WRITE(6,*) ' ISTR NSTA NSTB ',ISTR,NSTA,NSTB
C?          WRITE(6,*) ' NSTAI,NSTAK',NSTAI,NSTAK
            DO IB = 1, NSTB
              DO IA = 1, NSTA
C               IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
          END IF
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADAST_GASSM '
        WRITE(6,*) ' ======================== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Mapped strings and sign '
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
        CALL Unused_integer(NACGSOB)
      END IF
      END
