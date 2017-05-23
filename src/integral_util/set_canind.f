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
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine Set_CanInd
************************************************************************
*                                                                      *
*     Object: to set up a table which resolves the canonical           *
*             angular index, ip, into ix, iy, and iz.                  *
*                                                                      *
*     Author: Roland Lindh                                             *
*             Dept of Chem. Phys.                                      *
*             Univ. of Lund, Sweden                                    *
*             February 1999                                            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "ican.fh"
*
      ip=0
      Do la = 0, lab-1
         Do ix = la , 0, -1
            Do iy = la-ix, 0, -1
               iz = la-ix-iy
               ip = ip + 1
               iCan(1,ip)=ix
               iCan(2,ip)=iy
               iCan(3,ip)=iz
            End Do
         End Do
      End Do
*
      Return
      End
