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
      SUBROUTINE NGETH1(H,ISM,ITP,JSM,JTP)
*
* One-electron integrals over orbitals belonging to
* given OS class
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "detdim.fh"
#include "Input.fh"
#include "orbinp_mclr.fh"
*.Output
      DIMENSION H(*)
*
      NI = NTSOB(ITP,ISM)
      NJ = NTSOB(JTP,JSM)
      IJ = 0
      DO 100 J = 1, NJ
        DO 50 I = 1, NI
          IJ = IJ+1
          H(IJ) = GTH1EN(I,ITP,ISM,J,JTP,JSM)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
        CALL WRTMAT(H,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
