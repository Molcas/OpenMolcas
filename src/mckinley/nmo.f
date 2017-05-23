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
      Function nMo(mIrr)
      Implicit Integer (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
      nInt=0
      nA=0
      Do iS=0,nIrrep-1
       nA=nA+nAsh(is)
      End Do
      NMM=nA*(nA+1)/2
      nInt=nMM*(nMM+1)/2
      nMo=nInt
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(mIrr)
      End
