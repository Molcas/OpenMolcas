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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine NOphase(a,n)
      Implicit Real*8 (a-h,o-z)
      Dimension a(n,n)
      sign=1.0d0
      Do 100 i=1,n
         If(sign*a(1,i).lt.0.0d0) Then
            Do 110 j=1,n
               a(j,i)=-a(j,i)
110         Continue
         End If
         sign=-sign
100   Continue
      Return
      End
