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
      Subroutine VpMem(nRys,MemVp,la,lb,lr)
      Integer iAngV(4)
*
      nElem(i)=(i+1)*(i+2)/2
*
      Call mHrr(la,lb+1,nFlop,nMem)
*
      nRys=(la+lb+1+lr+2)/2
      iAngV(1) = la
      iAngV(2) = lb+1
      iAngV(3) = 0
      iAngV(4) = 0
      Call MemRys(iAngV,MemNA1)
*
      MemNA1 = Max(nMem,MemNA1)
*
      If (lb.ne.0) Then
         Call mHrr(la,lb-1,nFlop,nMem)
*
         nRys=(la+lb-1+lr+2)/2
         iAngV(1) = la
         iAngV(2) = lb-1
         iAngV(3) = 0
         iAngV(4) = 0
         Call MemRys(iAngV,MemNA2)
*
         MemNA2 = Max(nMem,MemNA2)
      Else
         MemNA2=0
      End If
*
      MemVp=Max(MemNA1,MemNA2)
*
      MemVp=MemVp+1
*
      MemVp = MemVp+nElem(la)*nElem(lb+1)
      If (lb.ne.0) MemVp=MemVp+nElem(la)*nElem(lb-1)
*
      Return
      End
