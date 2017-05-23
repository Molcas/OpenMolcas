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
      SubRoutine Cho_PreScr(Thr1,Thr2)
C
C     Purpose: read integral prescreening thresholds from common block.
C
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"

      Thr1 = CutInt
      Thr2 = ThrInt

      End
