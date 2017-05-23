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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      Function Gamma2(m,T)
************************************************************************
*                                                                      *
* Object: to compute the auxiliary function in the high argument       *
*         approximation.                                               *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Gamma2
*
      Gamma2= Sqrt(Two*ACos(Zero)/T)/Two
      Do 10 i = 1, m
         Gamma2 = ((Two*DBLE(i)-One)/(Two*T))*Gamma2
 10   Continue
      Return
      End
