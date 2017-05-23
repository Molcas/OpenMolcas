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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ORBORD(NSMOB,MXPOBS,NR4TP,NDEOBS,NINOBS,NR0OBS,NACOBS,
     &                  NRSOBS,NR4OBS,NOCOBS,NTOOBS,
     &                  IREOST,IREOTS,ISFTO,ITFSO,IPRNT,IBSO,
     &                  NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,MXPR4T,
     &                  ISMFSO,ITPFTO,NOBPT)
*
* Obtain Reordering arrays for orbitals
* ( See note below for assumed ordering )
*
* =====
* Input
* =====
*  NSMOB  : Number of orbital symmetries
*  MXPOBS : MAx number of orbital symmetries
*  NR4TP  : Number of RAS4 types
*  NDEOBS : Number of deleted orbitals per symmetry
*  NINOBS : Number of inactive  orbitals per symmetry
*  NR0OBS : Number of RAS 0 (core ) orbitals per symmetry
*  NACOBS : Number of Active orbitals per symmetry
*  NRSOBS : Number of orbitals per symmetry in RAS1,RAS2,RAS3
*  NR4OBS : Number of RAS 4 orbitals per symmetry and type
*  NOCOBS : Number of occupied orbitals per symmetry
*  NTOOBS : Number of orbitals per symmetry,all types
*
* ======
* Output
* ======
*  IREOST : Reordering array symmetry => type
*  IREOTS : Reordering array type     => symmetry
*  ISFTO  : Symmetry array for type ordered orbitals
*  ITFSO  : Type array for symmetry ordered orbitals( not activated )
*  IBSO   : First orbital of given symmetry ( symmetry ordered )
*  NTSOB  : Number of active orbitals of give RAS type and SYM
*  IBTSOB : Off set for active orbitals of given RAS type and SYM
*  ITSOB  : Orbitals of given RAS type and sym
*
*  NOBPTS : Number of orbitals per subtype and symmetry
*  IOBPTS : Off sets for orbitals of given subtype and symmetry
*           ordered according to input integrals
*
* ISMFSO  : Symmetry of orbitals, symmetry ordereing
* ITPFTO  : Type of orbital, type ordering
* Jeppe Olsen, Winter 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NDEOBS(*),NINOBS(*),NR0OBS(*),NACOBS(*),
     &          NRSOBS(MXPOBS,3),NR4OBS(MXPOBS,*),NOCOBS(*),NTOOBS(*)
*. Output
      DIMENSION IREOST(*),IREOTS(*),ISFTO(*),ITFSO(*),IBSO(*)
      DIMENSION ISMFSO(*),ITPFTO(*)
      DIMENSION NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
      DIMENSION NOBPTS(6+MXPR4T,*),IOBPTS(6+MXPR4T,*)
      DIMENSION NOBPT(6+MXPR4T)

* ==========================
* Note on order of orbitals
* ==========================
*
* The orbitals are supposed to be imported ordered symmetry-type
* ordered as
* Loop over symmetries of orbitals
*   Inactive  of this symmetry
*   Core      of this symmetry
*   RAS1      of this symmetry
*   RAS2      of this symmetry
*   RAS3      of this symmetry
*   Secondary of this symmetry
*   Deleted   of this symmetry
* End of loop over symmetries
*
* Internally the orbitals are reordered to type symmetry order
* where the outer loop os over types and the inner loop is
* over symmetries.The types are arranged as
*  Ras1
*  Ras2
*  Ras3
*  Core
*  Secondary
*  Inactive
*  Deleted orbitals
*
*
* Active orbitals
*
      IAC = 0
      IBSM  =0 ! dummy intitialize
      NPREVS=0 ! dummy initialize
      IORB  =0 ! dummy initialize
      DO 11 IRS = 1, 3
      DO 10 ISM = 1,NSMOB
        IF(ISM.EQ.1) THEN
          IBSM = 1
        ELSE
          IBSM = IBSM + NTOOBS(ISM-1)
        END IF
        IF(IRS.EQ.1) THEN
          NPREVS = NINOBS(ISM)+NR0OBS(ISM)
          IORB = NRSOBS(ISM,1)
        ELSE IF(IRS.EQ.2) THEN
          NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)
          IORB = NRSOBS(ISM,2)
        ELSE IF(IRS.EQ.3) THEN
          NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NRSOBS(ISM,1)+NRSOBS(ISM,2)
          IORB = NRSOBS(ISM,3)
        END IF
        DO 9 IIAC = 1, IORB
*. Type ordered index
          IAC = IAC + 1
*. Symmetry ordered index
          IACS = IBSM + NPREVS - 1 + IIAC
          ISFTO(IAC) = ISM
C         ISMFSO(IACS) = ISM
          ITPFTO(IAC) = IRS
*         ITFSO(IACS) = IRS
          IREOST(IACS) = IAC
          IREOTS(IAC) = IACS

    9   CONTINUE
   10 CONTINUE
   11 CONTINUE
      NACOB = IAC
C?    write(6,*) ' IAC ',IAC
*
*. RAS 0 orbitals
*
      IR0 = NACOB
      DO 20 ISM = 1,NSMOB
        IF(ISM.EQ.1) THEN
          IBSM = 1
        ELSE
          IBSM = IBSM + NTOOBS(ISM-1)
        END IF
        NPREVS = NINOBS(ISM)
        DO 19 IIR0 = 1, NR0OBS(ISM)
*. Type ordered index
          IR0 = IR0 + 1
*. Symmetry ordered index
          IR0S = IBSM + NPREVS - 1 + IIR0
          ISFTO(IR0) = ISM
          ITFSO(IR0S) = 4
          IREOST(IR0S) = IR0
          IREOTS(IR0) = IR0S
   19   CONTINUE
   20 CONTINUE
      NR0OB = IR0 - NACOB
C?    write(6,*) ' IR0 ',IR0
*
*. RAS 4 orbitals
*
      IR4 = NACOB+NR0OB
      DO 30 ISM = 1,NSMOB
        IF(ISM.EQ.1) THEN
          IBSM = 1
        ELSE
          IBSM = IBSM + NTOOBS(ISM-1)
        END IF
        NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)
        DO 29 ITP = 1, NR4TP
        DO 29 IIR4 = 1, NR4OBS(ISM,ITP)
*. Type ordered index
          IR4 = IR4 + 1
*. Symmetry ordered index
          IR4S = IBSM + NPREVS - 1 + IIR4
          ISFTO(IR4) = ISM
          ITFSO(IR4S) = 5
          IREOST(IR4S) = IR4
          IREOTS(IR4) = IR4S
   29   CONTINUE
   30 CONTINUE
      NR4OB = IR4 - NACOB - NR0OB
C?    write(6,*) ' IR4 ',IR4
*
*. Inactive orbitals
*
      IIN = NACOB+NR0OB+NR4OB
      DO 40 ISM = 1,NSMOB
        IF(ISM.EQ.1) THEN
          IBSM = 1
        ELSE
          IBSM = IBSM + NTOOBS(ISM-1)
        END IF
        NPREVS = 0
        DO 39 IIIN = 1, NINOBS(ISM)
*. Type ordered index
          IIN = IIN + 1
*. Symmetry ordered index
          IINS = IBSM + NPREVS - 1 + IIIN
          ISFTO(IIN) = ISM
          ITFSO(IINS) = 6
          IREOST(IINS) = IIN
          IREOTS(IIN) = IINS
   39   CONTINUE
   40 CONTINUE
      NINOB = IIN - NACOB - NR0OB - NR4OB
C?    write(6,*) ' IIN   ',IIN
*
*. Deleted orbitals
*
      IDE = NACOB+NR0OB+NR4OB+NINOB
      DO 50 ISM = 1,NSMOB
        IF(ISM.EQ.1) THEN
          IBSM = 1
        ELSE
          IBSM = IBSM + NTOOBS(ISM-1)
        END IF
        IR4 = 0
        DO 45 ITP = 1, NR4TP
          IR4 = IR4 + NR4OBS(ISM,ITP)
   45   CONTINUE

        NPREVS = NINOBS(ISM)+NR0OBS(ISM)+NACOBS(ISM)+IR4
        DO 49 IIDE = 1, NDEOBS(ISM)
*. Type ordered index
          IDE = IDE + 1
*. Symmetry ordered index
          IDES = IBSM + NPREVS - 1 + IIDE
          ISFTO(IDE) = ISM
          ITFSO(IDES) = 7
          IREOST(IDES) = IDE
          IREOTS(IDE) = IDES
   49   CONTINUE
   50 CONTINUE
      NTOOB = IDE
C?    write(6,*) ' IDE ', IDE
*
      IOFF = 1
      DO 100 ISM = 1, NSMOB
       IBSO(ISM) = IOFF
       IOFF = IOFF + NTOOBS(ISM)
  100 CONTINUE
*
* ==================
* NTSOB,IBTSOB,ITSOB
* ==================
*
      IOFF = 1
      DO 300 I123 = 1, 3
        DO 200 ISM = 1, NSMOB
          NTSOB(I123,ISM) = NRSOBS(ISM,I123)
          IBTSOB(I123,ISM) = IOFF
          CALL ISTVC2(ITSOB(IOFF),IOFF-1,1,NRSOBS(ISM,I123))
          IOFF = IOFF + NRSOBS(ISM,I123)
  200   CONTINUE
  300 CONTINUE
*
* =====================
* NOBPTS NOBPT IOBPTS
* =====================
*
*. Loop over types in input order
      Call iCopy(NR4TP+6,0,0,NOBPT,1)
      LORB  = 0 ! dummy initialize
      IOTYPE= 0 ! dummy initialize
      DO 2000 ISMOB = 1, NSMOB
        LSMOB = 0
        DO 1000 ITYPE = 1, NR4TP + 6
          IF(ITYPE.EQ.1) THEN
*.Inactive ( frozen in normal notation )
            LORB = NINOBS(ISMOB)
            IOTYPE = 5+NR4TP
          ELSE IF(ITYPE.EQ.2) THEN
*.RAS0 ( inactive in normal notation
            LORB = NR0OBS(ISMOB)
            IOTYPE = 4
          ELSE IF (ITYPE.EQ.3) THEN
*.RAS1
            LORB = NRSOBS(ISMOB,1)
            IOTYPE = 1
          ELSE IF (ITYPE.EQ.4) THEN
*.RAS2
            LORB = NRSOBS(ISMOB,2)
            IOTYPE = 2
          ELSE IF (ITYPE.EQ.5) THEN
*.RAS3
            LORB = NRSOBS(ISMOB,3)
            IOTYPE = 3
          ELSE IF (ITYPE.GE.6.AND.ITYPE.LE.6+NR4TP-1) THEN
*.RAS4
            LORB = NR4OBS(ISMOB,ITYPE-5)
            IOTYPE = ITYPE -1
          ELSE IF (ITYPE.EQ.6+NR4TP) THEN
*. deleted orbitals
            LORB = NDEOBS(ISMOB)
            IOTYPE = ITYPE
          END IF
          IOBPTS(IOTYPE,ISMOB) = LSMOB+1
          NOBPTS(IOTYPE,ISMOB) = LORB
          NOBPT(IOTYPE) = NOBPT(IOTYPE)+LORB
          LSMOB = LSMOB + LORB
 1000   CONTINUE
 2000 CONTINUE
*
* =======
* ISMFSO
* =======
*
      IORB = 0
      DO ISM = 1, NSMOB
        DO IOB = 1, NTOOBS(ISM)
          IORB = IORB + 1
          ISMFSO(IORB) = ISM
        END DO
      END DO
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Output from ORBORD '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Symmetry of orbitals , type ordered '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Symmetry => type reordering array '
        CALL IWRTMA(IREOST,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Type => symmetry reordering array '
        CALL IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' IBSO array '
        CALL IWRTMA(IBSO,1,NSMOB,1,NSMOB)
*
        WRITE(6,*) ' NTSOB array : '
        CALL IWRTMA(NTSOB,3,NSMOB,3,NSMOB)
        WRITE(6,*) ' IBTSOB array '
        CALL IWRTMA(IBTSOB,3,NSMOB,3,NSMOB)
        WRITE(6,*) ' ITSOB '
        CALL IWRTMA(ITSOB,1,NACOB,1,NACOB)
*
        WRITE(6,*) ' NOBPTS '
        CALL IWRTMA(NOBPTS,6+NR4TP,NSMOB,6+MXPR4T,MXPOBS)
        WRITE(6,*) ' NOBPT '
        CALL IWRTMA(NOBPTS,6+NR4TP,1,6+MXPR4T,1)
        WRITE(6,*) ' IOBPTS '
        CALL IWRTMA(IOBPTS,6+NR4TP,NSMOB,6+MXPR4T,MXPOBS)
*
        WRITE(6,*) ' ISFTO array : '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFSO array : '
        CALL IWRTMA(ITFSO,1,NTOOB,1,NTOOB)
*
        WRITE(6,*) ' ISMFSO array : '
        CALL IWRTMA(ISMFSO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITPFTO array : '
        CALL IWRTMA(ITPFTO,1,NTOOB,1,NTOOB)
      END IF
*

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(NOCOBS)
      END
