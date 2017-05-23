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
      SUBROUTINE KILLSCTAB(LSCTAB)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      NSIZE=IWORK(LSCTAB+0)
      ITYPE=IWORK(LSCTAB+1)
      IF(ITYPE.NE.47) THEN
        WRITE(6,*)' KILLSCTAB Error: Asked to kill innocent bystander!'
        WRITE(6,*)' Argument LSCTAB=',LSCTAB
        WRITE(6,*)' Table length   =',IWORK(LSCTAB)
        WRITE(6,*)' Table type ID  =',IWORK(LSCTAB+1)
        CALL ABEND()
      END IF
      NSIZE=IWORK(LSCTAB+0)
      ITYPE=IWORK(LSCTAB+1)
      LTRANS=IWORK(LSCTAB+6)
      NTRANS=IWORK(LSCTAB+7)
      CALL GETMEM('SpnCplTb','FREE','Inte',LSCTAB,NSIZE)
      CALL GETMEM('SpnCplCf','FREE','Real',LTRANS,NTRANS)
      RETURN
      END
