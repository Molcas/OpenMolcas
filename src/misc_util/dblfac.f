************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
* Copyright (C) 1990,2020, Roland Lindh                                *
************************************************************************
      Function DblFac(n)
************************************************************************
*                                                                      *
* Object: to compute the double factorial of n.                        *
*                                                                      *
* Called from: NrmSph                                                  *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 DblFac
*
      DblFac = One
      Do 20 i = n , 1, -2
         DblFac = DblFac * DBLE(i)
 20   Continue
      Return
      End
