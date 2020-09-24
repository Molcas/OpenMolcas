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
* Copyright (C) 1990,1991,1994, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RysEF(xyz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,
     &                 EFInt,meMin,meMax,mfMin,mfMax,Scrtch,PreFct,
     &                 AeqB, CeqD)
************************************************************************
*                                                                      *
*     Object: to compute integrals corresponding to the primitive set  *
*             used for the HRR. The primitive integrals are generated  *
*             from the 2D-integrals according to the Rys quadrature.   *
*                                                                      *
* Called from: Rys                                                     *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              RysEF0                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90.                                             *
*                                                                      *
*             Modified for kernal routine RysEF0 August '90.           *
*             Modified for kernal routines RysS1, RysS2, and RysS3     *
*             September '90.                                           *
*             Modified for improved vectorization August '91.          *
*             Modified for decreased memory access January '94.        *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "TriInd.fh"
#include "real.fh"
#include "print.fh"
      Real*8 xyz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg),
     &       Scrtch(nRys,mArg), EFInt(nArg,meMin:meMax,mfMin:mfMax)
      Logical AeqB, CeqD
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*80 Label
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 16
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      ne = (neMax+1)*(neMax+2)/2
      nf = (nfMax+1)*(nfMax+2)/2
*
      If (ne.gt.IJ_Max .or. nf.gt.IJ_Max) Then
         Write (6,*) 'ne,nf=',ne,nf
         Call WarningMessage(2,
     &               'Increase IJ_Max to the larger of the above!')
         Call Abend()
      End If
*
      Do ief = 1, ne*nf
         if = (ief-1)/ne + 1
         ie = ief - (if-1)*ne
*
         ixe = iTriInd(1,ie)
         iye = iTriInd(2,ie)
         ixye = ixe + iye
*
         ixf = iTriInd(1,if)
         iyf = iTriInd(2,if)
         ixyf = ixf + iyf
*                                                                      *
************************************************************************
*                                                                      *
         nzeMax=Max(0,neMax-ixe-iye)
         nzfMax=Max(0,nfMax-ixf-iyf)
         nzeMin=Max(0,neMin-ixe-iye)
         If (AeqB) nzeMin = nzeMax
         nzfMin=Max(0,nfMin-ixf-iyf)
         If (CeqD) nzfMin = nzfMax
*
         nItem=(nzeMax-nzeMin+1)*(nzfMax-nzfMin+1)
         If (nItem.gt.1) Then
*
*           Precompute for all arguments Ix*Iy, avoid multiplying
*           with ones.
*
*
*           Combine with possible Iz
*
            If (ixe+ixf+iye+iyf.eq.0) Then
*
               Call RysEF1(                     xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else If(ixe+ixf.eq.0) Then
*
               Call RysEF0(xyz2D(1,1,2,iye,iyf),xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else If(iye+iyf.eq.0) Then
*
               Call RysEF0(xyz2D(1,1,1,ixe,ixf),xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else
*
               Do iArg = 1, mArg
                  Do iRys=1,nRys
                     Scrtch(iRys,iArg) = xyz2D(iRys,iArg,1,ixe,ixf)
     &                                 * xyz2D(iRys,iArg,2,iye,iyf)
                  End Do
               End Do
               Call RysEF0(Scrtch,              xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            End If
*
         Else
*
*           Here if only one triplet of 2D-integrals
*
*           Contract over roots
*
            If (ixe+ixf+iye+iyf.eq.0) Then
*
               Call RysEF2(                     xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else If (ixe+ixf.eq.0) Then
*
               Call RysEF3(xyz2D(1,1,2,iye,iyf),xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else If (iye+iyf.eq.0) Then
*
               Call RysEF3(xyz2D(1,1,1,ixe,ixf),xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            Else
*
               Call RysEF4(                     xyz2D,
     &                     nArg,mArg,nRys,neMin,neMax,nfMin,
     &                     nfMax,EFInt,meMin,meMax,mfMin,mfMax,
     &                     PreFct,ixe,ixf,ixye,ixyf,
     &                     nzeMin,nzeMax,nzfMin,nzfMax)
*
            End If
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      Do iab = meMin, meMax
         Do icd = mfMin, mfMax
            Write (Label,'(A,I3,A,I3,A)') ' In RysEF: [', iab, ',0|',
     &                                      icd, ',0]'
            Call RecPrt(Label,' ',EFInt(1,iab,icd),1,nArg)
         End Do
      End Do
#endif
      Return
      End
