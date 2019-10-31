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
      SUBROUTINE mk_IOFF(IOFF,mSYM,NBASF,ISY12)
#include "symmul.fh"
      INTEGER IOFF(mSYM), NBASF(mSym)

C FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
      IOF=0
      Call IZERO(IOFF,8)
      DO ISY1=1,mSYM
        ISY2=MUL(ISY1,ISY12)
        IF(ISY1.LT.ISY2) GOTO 10
        IOFF(ISY1)=IOF
        IOFF(ISY2)=IOF
        NB1=NBASF(ISY1)
        NB2=NBASF(ISY2)
        NB12=NB1*NB2
        IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
        IOF=IOF+NB12
  10    CONTINUE
      END DO
*
      RETURN
      END
