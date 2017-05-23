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
* Copyright (C) 1997, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GET_SPGP_INF(ISPGP,ITP,IGRP)
*
* Obtain groups defining supergroup ISPGP,ITP
*
* Jeppe Olsen, November 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.fh"
#include "cgas.fh"
#include "gasstr.fh"
*. Output
      DIMENSION IGRP(*)
*
      NTEST = 00
*. Absolute group number
C?    WRITE(6,*) ' GET_SPGP_INF : ISPGP, ITP', ISPGP, ITP
      ISPGPABS = ISPGP + IBSPGPFTP(ITP) -1
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),IGRP,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_SPGP_INF : ISPGP ITP ISPGPABS',
     &              ISPGP, ITP, ISPGPABS
        WRITE(6,*) ' Groups of supergroups'
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
