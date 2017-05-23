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
      Function iPL_espf()
      Implicit Real*8(A-H,O-Z)
*
*     Returns the print level
*
#include "espf.fh"
*
      Logical Reduce_Prt
      External Reduce_Prt
*
      iPL = iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
      iPL_espf = iPL
      Return
      End
