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
* Copyright (C) Roland Lindh                                           *
************************************************************************
*  ip_of_iWork
*
*> @brief
*>   Return pointer to be used in \c iWork which corresponds to the first element of array \p A
*> @author R. Lindh
*>
*> @details
*> Function returns pointer such that \p A and \c iWork(ip_of_iWork(A)) refer to
*> the same memory location.
*>
*> @param[in] A first element of integer array
*>
*> @return pointer to \p A in \c iWork
************************************************************************
      Integer Function ip_of_iWork(A)
      Implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Integer A
*
      loc1=(iiLoc(A)-iiLoc(iWork(ip_iDummy)))
      loc2=(iiLoc(iWork(ip_iDummy+1))-iiLoc(iWork(ip_iDummy)))
      ip_of_iWork = ip_iDummy + loc1/loc2
*
      Return
      End
