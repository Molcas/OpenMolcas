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
      FUNCTION IFRMR(WORK,IROFF,IELMNT)
*
* An integer array is stored in real array WORK,
* starting from WORK(IROFF). Obtain element
* IELMNT of this array
*
      INTEGER WORK(*)
*
#include "irat.fh"
*. offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      IFRMR = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
