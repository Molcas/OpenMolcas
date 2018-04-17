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
*  ip_of_Work
*
*> @brief
*>   Return pointer to be used in \c Work which corresponds to the first element of array \p A
*> @author R. Lindh
*>
*> @details
*> Function returns pointer such that \p A and \c Work(ip_of_Work(A)) refer to
*> the same memory location.
*>
*> @param[in] A first element of real array
*>
*> @return pointer to \p A in \c Work
************************************************************************
      Integer Function ip_of_work(A)
      Implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 A
*
      loc1=(idLoc(A)-idLoc(Work(ip_Dummy)))
      loc2=(idLoc(Work(ip_Dummy+1))-idLoc(Work(ip_Dummy)))
      ip_of_Work = ip_Dummy + loc1/loc2
*
      Return
      End
