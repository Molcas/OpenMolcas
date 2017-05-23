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
* Copyright (C) 1990, Jeppe Olsen                                      *
************************************************************************
      FUNCTION ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
*
* Adress of string IOCC
*
* version of Winter 1990 , Jeppe Olsen
*
      INTEGER Z
      DIMENSION IOCC(*),NEWORD(*),Z(NORB,*)
*
      IZ = 1
      DO 100 I = 1,NEL
        IZ = IZ + Z(IOCC(I),I)
  100 CONTINUE
*
      IF(IREORD.EQ.0) THEN
        ISTRNM = IZ
      ELSE
        ISTRNM = NEWORD(IZ)
      END IF
*
      NTEST = 0
      IF ( NTEST .GT. 1 ) THEN
        WRITE(6,*) ' STRING'
        CALL IWRTMA(IOCC,1,NEL,1,NEL)
        WRITE(6,*) ' Z matrix '
        CALL IWRTMA(Z,NORB,NEL,NORB,NEL)
C       WRITE(6,*) ' First two elements of reorder array'
C       CALL IWRTMA(NEWORD,1,2,1,2)
        WRITE(6,*) ' ADRESS OF STRING ',ISTRNM
        WRITE(6,*) ' REV LEX number : ', IZ
      END IF
*
      RETURN
      END
