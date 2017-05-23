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
* Copyright (C) 2011, Giovanni Li Manni                                *
************************************************************************
      FUNCTION ICHECK_OCC_IN_ACCSPC(IOCC,IMINMAX,NGAS,MXPNGAS)
* Check if Occupation of GAS Spaces defined by IOCC are
* within the constraints of  IMINMAX chosen by the user
*
* Giovanni Li Manni 7Nov2011, for BK implementation
*
#include "implicit.fh"
*. Input
      INTEGER IOCC(NGAS)
      INTEGER IMINMAX(MXPNGAS,2)
*
      I_AM_IN = 1
      DO IGAS = 1, NGAS
       IF(IOCC(IGAS).LT.IMINMAX(IGAS,1) .OR.
     &    IOCC(IGAS).GT.IMINMAX(IGAS,2)) I_AM_IN=0
      END DO
      ICHECK_OCC_IN_ACCSPC = I_AM_IN
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IMINMAX'
       CALL IWRTMA(IMINMAX,NGAS,2,MXPNGAS,2)
      END IF
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IOCC'
       CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
       WRITE(6,*) ' And the verdict is ', I_AM_IN
      END IF
      RETURN
      END
