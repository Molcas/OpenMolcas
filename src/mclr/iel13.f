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
* Copyright (C) 1991,1994, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE IEL13(MNRS1,MXRS1,MNRS3,MXRS3,NELEC,NOCTYP,
     &           IEL1,IEL3,IEL123,IACTIV)
*
* A type of strings contains NEL electrons and
*
* The number of electrons in RAS 1 is between MNRS1 AND MXRS1
* The number of electrons in RAS 3 is between MNRS3 AND MXRS3
*
* Find the number of electrons in RAS1 and RAS3 in each subtype
* and determine if subtype  is active
*
* Jeppe Olsen , Spring 1991
*               IEL123 added Jan 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IEL1(*),IEL3(*),IEL123(3,*),IACTIV(*)
*
      Call iCopy(NOCTYP,0,0,IACTIV,1)
      Call iCopy(NOCTYP,0,0,IEL1,1)
      Call iCopy(NOCTYP,0,0,IEL3,1)
      DO 300 KEL3 = MNRS3,MXRS3
        DO 100 KEL1 = MNRS1,MXRS1
          ITYP = (MXRS1-KEL1) * (MXRS3-MNRS3+1 )
     &         + KEL3-MNRS3+1
            IEL1(ITYP) = KEL1
            IEL3(ITYP) = KEL3
*. Added
            IEL123(1,ITYP) = KEL1
            IEL123(2,ITYP) = NELEC - KEL1 - KEL3
            IEL123(3,ITYP) = KEL3
*. Added End
          IF(KEL1+KEL3.LE.NELEC) THEN
            IACTIV(ITYP) = 1
          END IF
  100   CONTINUE
  300 CONTINUE
*
      NTEST = 00
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ============== '
        WRITE(6,*) ' IEL13 speaking '
        WRITE(6,*) ' ============== '
        WRITE(6,'(A,4I3)')
     &  ' IEL1 IEL3 IACTIV for MNRS1 MXRS1 MNRS3 MXRS3 ',
     &    MNRS1,MXRS1,MNRS3,MXRS3
        DO 400 ITYP = 1, NOCTYP
          WRITE(6,'(3I3)') IEL1(ITYP),IEL3(ITYP),IACTIV(ITYP)
  400   CONTINUE
*
        WRITE(6,*) ' IEL123 matrix '
        CALL IWRTMA(IEL123,3,NOCTYP,3,NOCTYP)
      END IF
*
      RETURN
      END
