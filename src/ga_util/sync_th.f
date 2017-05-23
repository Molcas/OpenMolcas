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
      Subroutine Sync_TH(TwoHam,nDens)
      Implicit Real*8 (a-h,o-z)
      Real*8 TwoHam(nDens)
*
#include "para_info.fh"
*
#ifdef _MOLCAS_MPP_
      If (.Not. Is_Real_Par()) Return
      If (nProcs.eq.1) Return
      Call BCTwoHam(TwoHam,nDens,TCPU,TWall)
      Call SavTim(10,TCPU,TWall)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(TwoHam)
#endif
*
      Return
      End
