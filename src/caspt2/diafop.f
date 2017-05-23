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
      SUBROUTINE DIAFOP(NGRP,FOPXMS,SCR,EVEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FOPXMS(NGRP,NGRP)
      DIMENSION SCR((NGRP*(NGRP+1))/2)
      DIMENSION EVEC(NGRP,NGRP)

      IJ=0
      DO I=1,NGRP
       DO J=1,I
        IJ=IJ+1
        SCR(IJ)=FOPXMS(I,J)
       END DO
      END DO
      CALL DCOPY_(NGRP**2,0.0D0,0,EVEC,1)
      CALL DCOPY_(NGRP,1.0D0,0,EVEC,NGRP+1)
      CALL JACOB(SCR,EVEC,NGRP,NGRP)

      RETURN
      END
