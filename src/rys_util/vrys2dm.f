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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine vRys2Dm(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,
     &                 B10,laa,B00,lac,B01,lcc,
     &                 la,lb,lc,ld,IfGrad)
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
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
      Real*8 xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax),
     &       PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3),
     &       B10(nArg*lRys), B00(nArg*lRys),
     &       B01(nArg*lRys)
      Logical IfGrad(3,4)
#ifdef _DEBUGPRINT_
      Character*30 Label
#endif
*
      iRout = 15
      iPrint = nPrint(iRout)
*     Call QEnter('Rys2Dm')
#ifdef _DEBUGPRINT_
      iPrint=99
      If (iPrint.ge.99) Then
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
      call dcopy_(2*nArg*lRys,[One],0,xyz2D(1,1,0,0),1)
*
*     Compute 2D integrals with index (i,0)
*
      Do 200 iCar = 1, 3
      llab = 0
      If (IfGrad(iCar,1).or.IfGrad(iCar,2)) llab = 1
      llcd = 0
      If (IfGrad(iCar,3).or.IfGrad(iCar,4)) llcd = 1
      mabMax = la + lb + llab
      mcdMax = lc + ld + llcd
*
      If (mabMax.ne.0) Then
         Do 201 i = 1, nArg*lRys
            xyz2D(i,iCar,1,0) = PAWP(i,iCar) * xyz2D(i,iCar,0,0)
 201     Continue
      End If
      If (mabMax-1.eq.1) Then
         Do 210 i = 1, nArg*lRys
            xyz2D(i,iCar,2,0) = PAWP(i,iCar) * xyz2D(i,iCar,1,0)
     &                       + B10(i) * xyz2D(i,iCar,0,0)
 210     Continue
      Else If (mabMax-1.gt.1) Then
         Fact = One
         Do 250 iab = 1, mabMax-1
            Do 260 i = 1, nArg*lRys
               temp1 = PAWP(i,iCar) * xyz2D(i,iCar,iab,0)
               temp2 = Fact * B10(i) * xyz2D(i,iCar,iab-1,0)
               xyz2D(i,iCar,iab+1,0) = temp1 + temp2
 260        Continue
            Fact = Fact + One
 250     Continue
      End If
*
*     Compute 2D integrals with index (0,i)
*
      If (mcdMax.ne.0) Then
         Do 301 i = 1, nArg*lRys
            xyz2D(i,iCar,0,1) = QCWQ(i,iCar) * xyz2D(i,iCar,0,0)
 301     Continue
      End If
      If (mcdMax-1.eq.1) Then
         Do 310 i = 1, nArg*lRys
            xyz2D(i,iCar,0,2) = QCWQ(i,iCar) * xyz2D(i,iCar,0,1)
     &                       + B01(i     ) * xyz2D(i,iCar,0,0)
 310     Continue
      Else If (mcdMax-1.gt.1) Then
         Fact = One
         Do 350 icd = 1, mcdMax-1
            Do 360 i = 1, nArg*lRys
               temp1 = QCWQ(i,iCar) * xyz2D(i,iCar,0,icd)
               temp2 = Fact * B01(i     ) * xyz2D(i,iCar,0,icd-1)
               xyz2D(i,iCar,0,icd+1) = temp1 + temp2
 360        Continue
            Fact = Fact + One
 350     Continue
      End If
*
*     Compute 2D integrals with index (i,iCar,j)
*
      If (mcdMax.le.mabMax) Then
         Fac1 = One
         Do 400 icd = 1, mcdMax
            Do 425 i = 1, nArg*lRys
               xyz2D(i,iCar,1,icd)=PAWP(i,iCar)*xyz2D(i,iCar,0,icd)
     &                     + Fac1 * B00(i     )*xyz2D(i,iCar,0,icd-1)
 425        Continue
            If (mabMax-1.eq.1) Then
               Do 420 i = 1, nArg*lRys
                  xyz2D(i,iCar,2,icd)=PAWP(i,iCar)*xyz2D(i,iCar,1,icd)
     &                              +  B10(i     )*xyz2D(i,iCar,0,icd)
     &                         + Fac1 *B00(i     )*xyz2D(i,iCar,1,icd-1)
 420           Continue
            Else If (mabMax-1.gt.1) Then
               Fac2 = One
               Do 450 iab = 1, mabMax-1
                  Do 460 i = 1, nArg*lRys
                     temp1 = PAWP(i,iCar) * xyz2D(i,iCar,iab,icd)
                     temp2 = Fac2 *B10(i     ) *xyz2D(i,iCar,iab-1,icd)
                     temp3 = Fac1 *B00(i     ) *xyz2D(i,iCar,iab,icd-1)
                     xyz2D(i,iCar,iab+1,icd) = temp1 + temp2 + temp3
 460              Continue
                  Fac2 = Fac2 + One
 450           Continue
            End If
            Fac1 = Fac1 + One
 400     Continue
      Else
         Fac1 = One
         Do 500 iab = 1, mabMax
            Do 525 i = 1, nArg*lRys
               xyz2D(i,iCar,iab,1)=QCWQ(i,iCar)*xyz2D(i,iCar,iab,0)
     &                      + Fac1 *B00(i     )*xyz2D(i,iCar,iab-1,0)
 525        Continue
            If (mcdMax-1.eq.1) Then
               Do 520 i = 1, nArg*lRys
                  xyz2D(i,iCar,iab,2)=QCWQ(i,iCar)*xyz2D(i,iCar,iab,1)
     &                              + B01(i     ) *xyz2D(i,iCar,iab,0)
     &                         + Fac1 *B00(i     )*xyz2D(i,iCar,iab-1,1)
 520           Continue
            Else If (mcdMax-1.gt.1) Then
               Fac2 = One
               Do 550 icd = 1, mcdMax-1
                  Do 560 i = 1, nArg*lRys
                     temp1 = QCWQ(i,iCar) *xyz2D(i,iCar,iab,icd)
                     temp2 = Fac2 *B01(i     ) *xyz2D(i,iCar,iab,icd-1)
                     temp3 = Fac1 *B00(i     ) *xyz2D(i,iCar,iab-1,icd)
                     xyz2D(i,iCar,iab,icd+1) = temp1 + temp2 + temp3
 560              Continue
                  Fac2 = Fac2 + One
 550           Continue
            End If
            Fac1 = Fac1 + One
 500     Continue
      End If
 200  Continue
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Do 600 iab = 0, nabMax
            Do 610 icd = 0, ncdMax
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
               Call RecPrt(Label,' ',xyz2D(1,1,iab,icd),nArg,lRys)
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
               Call RecPrt(Label,' ',xyz2D(1,2,iab,icd),nArg,lRys)
               Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
               Call RecPrt(Label,' ',xyz2D(1,3,iab,icd),nArg,lRys)
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
*     Call QExit('Rys2Dm')
      Return
      End
