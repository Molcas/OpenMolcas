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
      Subroutine Set_binom
************************************************************************
*                                                                      *
*     Object: to set up a table of binomial coefficients.              *
*                                                                      *
*     Author: Roland Lindh                                             *
*             Dept of Chem. Phys.                                      *
*             Univ. of Lund, Sweden                                    *
*             February 1999                                            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "binom.fh"
#include "real.fh"
*
      Do i = 0, 2*iTabMx
         Do j = -1, 2*iTabMx
            binom(i,j)=Zero
         End Do
      End Do
      binom(0,0)=One
*
*---- Do recurrence according to Pascal's triangle!
*
      Do n = 1, 2*iTabMx
         Do i = 0, n
            binom(n,i)=binom(n-1,i-1)+binom(n-1,i)
         End Do
      End Do
*
      Return
      End
