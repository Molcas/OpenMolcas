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
      Subroutine GR_DArray(Array,nArray)
      Implicit Real*8 (a-h,o-z)
#include "para_info.fh"
*
      Real*8 Array(nArray)
*
#ifdef _MOLCAS_MPP_
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
      Call CWTime(TCpu1,TWall1)
      Call GADGOP(Array,nArray,'+')
      Call CWTime(TCpu2,TWall2)
      Call SavTim(10,TCpu2-TCpu1,TWall2-TWall1)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Array)
#endif
*
      Return
      End
