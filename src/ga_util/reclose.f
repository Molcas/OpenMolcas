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
      Subroutine ReClose()
      Implicit None
#include "WrkSpc.fh"
#include "timtra.fh"
#include "para_info.fh"
*                                                                      *
************************************************************************
* Deallocate the arrays allocated by ReStart                           *
************************************************************************
*                                                                      *
      If (nfld_tim.ne.0)
     &  Call GetMem('iGATim','Free','Real',iGATim,nProcs*nfld_tim*2)
      If (nfld_stat.ne.0)
     &  Call GetMem('iGAStat','Free','Real',iGAStat,nProcs*nfld_stat)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
