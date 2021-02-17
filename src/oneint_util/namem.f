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
      Subroutine NAMem(
#define _CALLING_
#include "mem_interface.fh"
     &)
      use Basis_Info
#include "mem_interface.fh"
      Integer iAngV(4)
*
*                                                                      *
************************************************************************
*                                                                      *
      Call mHrr(la,lb,nFlop,nMem)
*
      iAngV(1) = la
      iAngV(2) = lb
      iAngV(3) = lr
      iAngV(4) = 0
      Call MemRys(iAngV,Mem)
      nHer=(la+lb+lr+2)/2
      If (Nuclear_Model.eq.mGaussian_Type) Then
*
         labcd = (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
*
         iAngV(3)=lr+2
         Call MemRys(iAngV,MemNA2)
         Mem=Max(Mem,MemNA2)
         nHer=(la+lb+lr+4)/2
         Mem=Mem+labcd
      End If
*
      Mem = Max(nMem,Mem)
*
      Return
      End
