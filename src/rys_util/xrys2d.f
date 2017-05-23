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
* Copyright (C) 1990,1991,1995, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine XRys2D(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,
     &                 B10,laa,B00,lac,B01,lcc)
************************************************************************
*                                                                      *
*     Object: to compute the 2-dimensional integrals of the Rys        *
*             quadrature. The z components are assumed to be pre-      *
*             conditioned with the weights of the roots of the         *
*             Rys polynomial.                                          *
*                                                                      *
* Called from: Rys                                                     *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
* Modified loop structure for RISC 1991 R. Lindh Dept. of Theoretical  *
* Chemistry, University of Lund, Sweden.                               *
* Modified for external field version, Feb '95.                        *
* VV: improve loop structure                                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
      Real*8 xyz2D(nArg*lRys*3,0:nabMax,0:ncdMax),
     &       PAWP(nArg*lRys*3), QCWQ(nArg*lRys*3),
     &       B10(nArg*lRys*3), B00(nArg*lRys*3), B01(nArg*lRys*3)
#ifdef _DEBUG_
      Character*30 Label
#endif
*
      iRout = 15
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      If (iPrint.ge.59) Then
         If (nabMax.gt.0) Call RecPrt('PAWP',' ',PAWP,nArg,lRys*3)
         If (ncdMax.gt.0) Call RecPrt('QCWQ',' ',QCWQ,nArg,lRys*3)
         If (laa.ne.0) Call RecPrt(' B10',' ',B10,nArg*lRys,3)
         If (lac.ne.0) Call RecPrt(' B00',' ',B00,nArg*lRys,3)
         If (lcc.ne.0) Call RecPrt(' B01',' ',B01,nArg*lRys,3)
      End If
#endif
*
*     Compute 2D integrals with index (0,0). Observe that the z
*     component already contains the weight factor.
*
      call dcopy_(2*nArg*lRys,One,0,xyz2D(1,0,0),1)
*
*---- Span first I(i,0)
*
      If (nabMax.ge.1) Then
         Do i = 1, nArg*lRys*3
            xyz2D(i,1,0) = PAWP(i)*xyz2D(i,0,0)
         End Do
      End If
      Do iab = 1, nabMax-1
         Do i = 1, nArg*lRys*3
            xyz2D(i,iab+1,0) = PAWP(i)*xyz2D(i,iab,0)
     &                       + Dble(iab)*B10(i)*xyz2D(i,iab-1,0)
         End Do
      End Do
*
*---- Now do the rest!
*
      If (ncdMax.ge.1) Then
         Do i = 1, nArg*lRys*3
            xyz2D(i,0,1)=QCWQ(i)*xyz2D(i,0,0)
         End Do
         Do iab = 1, nabMax
            Do i = 1, nArg*lRys*3
               xyz2D(i,iab,1)=QCWQ(i)*xyz2D(i,iab,0)
     &                       +Dble(iab)*B00(i)*xyz2D(i,iab-1,0)
            End Do
         End Do
      End If
      Do in = 1, ncdMax-1
         Do i = 1, nArg*lRys*3
            xyz2D(i,0,in+1)=QCWQ(i)*xyz2D(i,0,in)
     &                     -Dble(in)*B01(i)*xyz2D(i,0,in-1)
         End Do
         Do iab = 1, nabMax
            Do i = 1, nArg*lRys*3
               xyz2D(i,iab,in+1)=QCWQ(i)*xyz2D(i,iab,in)
     &                          +Dble(iab)*B00(i)*xyz2D(i,iab-1,in)
     &                          -Dble( in)*B01(i)*xyz2D(i,iab,in-1)
            End Do
         End Do
      End Do
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Write (6,*) ' 2D-integral computed in XRys2D'
         Do 600 iab = 0, nabMax
            Do 610 icd = 0, ncdMax
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
               Call RecPrt(Label,' ',
     &                     xyz2D(1,iab,icd),nArg,lRys)
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
               Call RecPrt(Label,' ',
     &                     xyz2D(1+nArg*lRys,iab,icd),nArg,lRys)
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
               Call RecPrt(Label,' ',
     &                     xyz2D(1+2*nArg*lRys,iab,icd),nArg,lRys)
 610        Continue
 600     Continue
      End If
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(laa)
         Call Unused_integer(lac)
         Call Unused_integer(lcc)
      End If
#endif
      Return
      End
