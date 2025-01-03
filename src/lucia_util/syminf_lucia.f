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
      SUBROUTINE SYMINF_LUCIA(IPRNT)
*
* Information about number of symmetries
*
      use lucia_data, only: PNTGRP,NIRREP
      IMPLICIT None
      Integer IPRNT
*
      IF(PNTGRP.EQ.1) THEN
* =====
* D 2 h
* =====
        CALL ZSYM1(NIRREP,IPRNT)
      ELSE
        WRITE(6,*) ' You are too early , sorry '
        WRITE(6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
*        STOP 11
        CALL SYSABENDMSG('lucia_util/syminf','Internal error',' ')
      END IF
*
      END SUBROUTINE SYMINF_LUCIA
