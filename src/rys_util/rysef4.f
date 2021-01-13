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
* Copyright (C) 1990,1994, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RysEF4(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,
     &                  nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                  PreFct,ixe,ixf,ixye,ixyf,
     &                  nzeMin,nzeMax,nzfMin,nzfMax)
************************************************************************
*                                                                      *
* Object: kernel routine to assemble the integrals from the Ixy        *
*         and Iz integrals.                                            *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             August '90                                               *
*                                                                      *
*             Modified for decreased memory access January '94.        *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 xyz2D(nRys,mArg,3,0:neMax,0:nfMax),
     &       PreFct(mArg), EFInt(nArg,meMin:meMax,mfMin:mfMax)
*
*     Statement function to compute canonical index
*
      iCan(ixyz,ix,iz) = ixyz*(ixyz+1)*(ixyz+2)/6 +
     &   (ixyz-ix)*(ixyz-ix+1)/2 + iz
*
      iyf=ixyf-ixf
      iye=ixye-ixe
      izf=nzfMax
      ize=nzeMax
      Indf=iCan(ixyf+izf,ixf,izf)
      Inde=iCan(ixye+ize,ixe,ize)
      If (nRys.eq.1) Then
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) = PreFct(iArg) *
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
         End Do
       Else If (nRys.eq.2) Then
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) = PreFct(iArg) * (
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
     &                       + xyz2D(2,iArg,1,ixe,ixf) *
     &                         xyz2D(2,iArg,2,iye,iyf) *
     &                         xyz2D(2,iArg,3,ize,izf))
         End Do
       Else If (nRys.eq.3) Then
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) = PreFct(iArg) * (
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
     &                       + xyz2D(2,iArg,1,ixe,ixf) *
     &                         xyz2D(2,iArg,2,iye,iyf) *
     &                         xyz2D(2,iArg,3,ize,izf)
     &                       + xyz2D(3,iArg,1,ixe,ixf) *
     &                         xyz2D(3,iArg,2,iye,iyf) *
     &                         xyz2D(3,iArg,3,ize,izf))
         End Do
      Else If (nRys.eq.4) Then
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) = PreFct(iArg) * (
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
     &                       + xyz2D(2,iArg,1,ixe,ixf) *
     &                         xyz2D(2,iArg,2,iye,iyf) *
     &                         xyz2D(2,iArg,3,ize,izf)
     &                       + xyz2D(3,iArg,1,ixe,ixf) *
     &                         xyz2D(3,iArg,2,iye,iyf) *
     &                         xyz2D(3,iArg,3,ize,izf)
     &                       + xyz2D(4,iArg,1,ixe,ixf) *
     &                         xyz2D(4,iArg,2,iye,iyf) *
     &                         xyz2D(4,iArg,3,ize,izf))
         End Do
      Else If (nRys.eq.5) Then
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) = PreFct(iArg) * (
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
     &                       + xyz2D(2,iArg,1,ixe,ixf) *
     &                         xyz2D(2,iArg,2,iye,iyf) *
     &                         xyz2D(2,iArg,3,ize,izf)
     &                       + xyz2D(3,iArg,1,ixe,ixf) *
     &                         xyz2D(3,iArg,2,iye,iyf) *
     &                         xyz2D(3,iArg,3,ize,izf)
     &                       + xyz2D(4,iArg,1,ixe,ixf) *
     &                         xyz2D(4,iArg,2,iye,iyf) *
     &                         xyz2D(4,iArg,3,ize,izf)
     &                       + xyz2D(5,iArg,1,ixe,ixf) *
     &                         xyz2D(5,iArg,2,iye,iyf) *
     &                         xyz2D(5,iArg,3,ize,izf))
         End Do
      Else
*
*--------------General code
*
         Do iArg = 1, mArg
            EFInt(iArg,Inde,Indf) =
     &                         xyz2D(1,iArg,1,ixe,ixf) *
     &                         xyz2D(1,iArg,2,iye,iyf) *
     &                         xyz2D(1,iArg,3,ize,izf)
               Do iRys = 2, nRys
               EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf) +
     &                            xyz2D(iRys,iArg,1,ixe,ixf) *
     &                            xyz2D(iRys,iArg,2,iye,iyf) *
     &                            xyz2D(iRys,iArg,3,ize,izf)
               End Do
               EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf) *
     &           PreFct(iArg)
         End Do
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(neMin)
         Call Unused_integer(nfMin)
         Call Unused_integer(nzeMin)
         Call Unused_integer(nzfMin)
      End If
      End
