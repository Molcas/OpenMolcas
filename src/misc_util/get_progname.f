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
      Function Get_ProgName()
      Character*100 Get_ProgName
#include "unixinfo.fh"
*
      Get_ProgName=ProgName(1:100)
*
      Return
      End
      Function Get_SuperName()
      Character*100 Get_SuperName
#include "unixinfo.fh"
*
      Get_SuperName=SuperName(1:100)
*
      Return
      End
