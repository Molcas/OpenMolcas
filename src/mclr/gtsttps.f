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
      SUBROUTINE GTSTTPS(IEL1,IEL3,NEL1,NEL3,NTYP,ITYP,IWAY)
*
* ITYP : type of strings with IEL1,IEL3 electrons
*
* IWAY = 1 : IEL1, IEL3 known, find ITYP
* IWAY = 2 : ITYP known, find IEL1, IEL3
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      DIMENSION NEL1(*),NEL3(*)
*
      IF(IWAY.EQ.1) THEN
        ITYP = -1
        DO 10 IITYP = 1, NTYP
          IF(IEL1.EQ.NEL1(IITYP).AND.IEL3.EQ.NEL3(IITYP))
     &    ITYP = IITYP
   10   CONTINUE
*
C       IF(ITYP .EQ. -1 ) THEN
C         WRITE(6,*) ' Error in GSTTPS '
C         WRITE(6,*) ' Error : Type could not be identified'
C         WRITE(6,*) ' Error : IEL1 IEL3 : , IEL1,IEL3 '
C         WRITE(6,*) ' I am going to STOP '
C         STOP'GSTTPS'
C       END IF
      ELSE IF (IWAY .EQ. 2 ) THEN
        IEL1 = NEL1(ITYP)
        IEL3 = NEL3(ITYP)
      END IF
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,5I4)') ' GSTTPS : IWAY IEL1 IEL3 ITYP ',
     &  IWAY,IEL1,IEL3,ITYP
      END IF
*
      RETURN
      END
