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
      Logical Function Is_Real_Par()
CSVC: determine if multiple processes are working together
      Implicit None
#ifdef _MOLCAS_MPP_
#  include "mpp_info.fh"
      Is_Real_Par = mpp_workshare
#else
      Is_Real_Par = .FALSE.
#endif
      End
