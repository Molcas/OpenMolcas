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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine CtrlMO(moip,nAcO)
*
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "etwas.fh"
      Integer moip(0:nIrrep-1)
*
      jAsh=0
      ii=0
      iTot=0
      Do iIrrep=0,nIrrep-1
            moip(iIrrep)=iTot
            iTot=iTot+nAsh(iIrrep)
      End Do
      nACO=iTot
      Return
      End
