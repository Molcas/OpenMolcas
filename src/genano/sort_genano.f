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
      Subroutine Sort_genano(eval,evec,n,nb)
      Implicit Real*8 (a-h,o-z)
      Dimension eval(n),evec(nb,n)
      Do 100 i=1,n-1
         k=i
         Do 110 j=i+1,n
            If(eval(j).gt.eval(k)) k=j
110      Continue
         If(k.ne.i) Then
            swap    = eval(k)
            eval(k) = eval(i)
            eval(i) = swap
            Do 120 l=1,nb
               swap      = evec(l,k)
               evec(l,k) =-evec(l,i)
               evec(l,i) = swap
120         Continue
         End If
100   Continue
      Return
      End
