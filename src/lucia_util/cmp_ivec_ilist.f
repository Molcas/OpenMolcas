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
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE CMP_IVEC_ILIST(IVEC,ILIST,LLIST,NLIST,INUM)
* An integer IVEC of LLIST entries are given.
* compare with list of vectors in ILIST and find first
* vector in LLIST that is identical to IVEC.
*
* If INUM = 0, the list was not found
*
*  Jeppe Olsen, December 2001
*
#include "implicit.fh"
*. General input
      INTEGER ILIST(LLIST,NLIST)
*. Specific input
      INTEGER IVEC(LLIST)
*
      INUM = 0
      DO JLIST = 1, NLIST
        IFOUND = 1
        DO IELMNT = 1, LLIST
          IF(IVEC(IELMNT).NE.ILIST(IELMNT,JLIST))IFOUND = 0
        END DO
        IF(IFOUND.EQ.1) INUM = JLIST
        IF(INUM.NE.0) GOTO 1001
      END DO
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input list : '
        CALL IWRTMA(IVEC,1,LLIST,1,LLIST)
        WRITE(6,*) ' Address of list : ', INUM
      END IF
*
      RETURN
      END
