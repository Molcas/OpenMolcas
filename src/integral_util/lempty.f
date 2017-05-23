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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Logical Function lEmpty(Coeff,n2,ld2,m2)
************************************************************************
*                                                                      *
* Object: to set if partial or whole contraction matrix is empty.      *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             June '91                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Coeff(ld2,m2)
#include "real.fh"
*
      lEmpty = .True.
      Temp = Zero
      Do 10 i = 1, n2
         Do 20 j = 1, m2
            Temp = Temp + Abs(Coeff(i,j))
 20      Continue
 10   Continue
      lEmpty = Temp.eq.Zero
*
      Return
      End
