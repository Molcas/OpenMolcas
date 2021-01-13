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
      SUBROUTINE ZSMCL(NSMST,NOCTP,NSSO,ISTSM,ISTCL)
*
* set symmetry and class arrays for strings
*
      INTEGER ISTSM(*),ISTCL(*),NSSO(NOCTP,NSMST)
*
      IOFF = 1
      DO 100 ISM = 1, NSMST
      DO 101 ICL = 1, NOCTP
        CALL  ISETVC(ISTSM(IOFF),ISM,NSSO(ICL,ISM))
        CALL  ISETVC(ISTCL(IOFF),ICL,NSSO(ICL,ISM))
        IOFF = IOFF + NSSO(ICL,ISM)
  101 CONTINUE
  100 CONTINUE
*
      RETURN
      END
