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
      INTEGER FUNCTION ISYMST(STRING,NEL)
*
* Master routine for symmetry of string
*
      use lucia_data, only: PNTGRP
      IMPLICIT NONE
      INTEGER NEL
      INTEGER STRING(*)

      INTEGER, EXTERNAL :: ISYMS1
      IF(PNTGRP.EQ.1) THEN
*.D2h
        ISYMST = ISYMS1(STRING,NEL)
      ELSE
        WRITE(6,*) ' Sorry PNTGRP option not programmed ', PNTGRP
        WRITE(6,*) ' Enforced stop in ISYMST '
*        STOP 5
         CALL SYSABENDMSG('lucia_util/isymst','Internal error',' ')
         ISYMST=-9999
      END IF
*
      END FUNCTION ISYMST
