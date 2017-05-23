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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine Binom(n,i,iBin)
************************************************************************
*                                                                      *
* Object: to compute the binomial factor                               *
*                                                                      *
* Called from: Azmthl                                                  *
*              R2N                                                     *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit real*8 (a-h,o-z)
      Num = 1
      iDen = 1
      Do 10 j = 1, i
         Num  = Num * (n-j+1)
         iDen = iDen * j
 10   Continue
      iBin = Num/iDen
      Return
      End
