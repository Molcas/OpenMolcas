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
      Subroutine PrHead(str,iOpt)
      Character*(*) str
      Character*72 line
      Character*1 fill
      n=Len(str)
      line(1:1)='*'
      line(72:72)='*'
      fill=' '
*----------------------------------------------------------------------*
      If(iAnd(iOpt,1).ne.0) Then
         fill='*'
         n=0
      End If
*----------------------------------------------------------------------*
      k=(73-n)/2
      m=72-k-n
      Do 100 i=2,k
         line(i:i)=fill
100   Continue
      Do 200 i=k+n+1,71
         line(i:i)=fill
200   Continue
      If(n.gt.0) line(k+1:k+n)=str
*----------------------------------------------------------------------*
      Write(6,*) line
      Return
      End
