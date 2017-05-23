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
      SUBROUTINE SYMCOM(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* IOBJ = 1 : I1,I2 are strings I12 determinant
* ( Other things can follow )
* IOBJ = 2 : I1,I2,I3 are externals
* IOBJ = 3 : I1 is an external, I2,I3 are dets
* IOBJ = 4 : I1 is orbital, I2 is string,l, I12 is string
* IOBJ = 5 : I1 is single excitation, I2 is string,l, I12 is string
* IOBJ = 6 : I1 is orbital, I2 is Orbital I12 is single excitation
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
* ================
*. Driver routine
* ================
#include "mxpdim.fh"
#include "lucinp.fh"
*
      IF(PNTGRP.EQ.1) THEN
        CALL SYMCM1(ITASK,IOBJ,I1,I2,I12)
      ELSE
        WRITE(6,*) ' PNTGRP parameter out of bounds ', PNTGRP
        WRITE(6,*) ' Enforced stop in SYMCOM '
        CALL SYSABENDMSG('lucia_util/symcom','Internal error',' ')
      END IF
*
      RETURN
      END
