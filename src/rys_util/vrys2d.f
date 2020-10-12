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
      SubRoutine vRys2D(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,
     &                 B10,laa,B00,lac,B01,lcc)
************************************************************************
*                                                                      *
* Object: to compute the 2-dimensional integrals of the Rys            *
*         quadrature. The z components are assumed to be pre-          *
*         conditioned with the weights of the roots of the             *
*         Rys polynomial.                                              *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
* Modified loop structure for RISC 1991 R. Lindh Dept. of Theoretical  *
* Chemistry, University of Lund, Sweden.                               *
* Further modifications in Jan-Feb. 1994.                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 xyz2D(nArg*lRys*3,0:nabMax,0:ncdMax),
     &       PAWP(nArg*lRys*3), QCWQ(nArg*lRys*3),
     &       B10(nArg*lRys), B00(nArg*lRys),
     &       B01(nArg*lRys)
      Logical lPAWP, lQCWQ
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*30 Label
*
      If (nabMax.gt.0) Call RecPrt('PAWP',' ',PAWP,lRys,nArg*3)
      If (ncdMax.gt.0) Call RecPrt('QCWQ',' ',QCWQ,lRys,nArg*3)
      If (laa.ne.0) Call RecPrt(' B10',' ',B10,lRys,nArg)
      If (lac.ne.0) Call RecPrt(' B00',' ',B00,lRys,nArg)
      If (lcc.ne.0) Call RecPrt(' B01',' ',B01,lRys,nArg)
#endif
*
      iOffy =   nArg*lRys
      iOffz = 2*nArg*lRys
      tTwo = Two
      If (nabMax.eq.0 .and. ncdMax.eq.0) Then
      Else
*
*--------General code
*
*-----Store away PAWPz and QCPQz
*
      lPAWP = nabMax.gt.2
      lQCWQ = ncdMax.gt.2
      If (nabMax.ge.1 .and. ncdMax.ge.1) Then
         If (ncdMax.le.nabMax) Then
            lPAWP = .True.
         Else
            lQCWQ = .True.
         End If
      End If
      If (lPAWP)
     &      call dcopy_(nArg*lRys,PAWP(1+iOffz),1,xyz2D(1      ,0,0),1)
      If (lQCWQ)
     &      call dcopy_(nArg*lRys,QCWQ(1+iOffz),1,xyz2D(1+iOffy,0,0),1)
*
*     Compute 2D integrals with index (i,0)
*
      If (nabMax.eq.0) Then
*
      Else If (nabMax.eq.1) Then
         Do 2000 i = 1, nArg*lRys
            PAWPz = xyz2D(i+iOffz,1,0)
            xyz2D(i+iOffz,1,0) = PAWPz*xyz2D(i+iOffz,0,0)
 2000    Continue
      Else If (nabMax.ge.2) Then
         Do 2001 i = 1, nArg*lRys
            PAWPx = xyz2D(i,1,0)
            xyz2D(i      ,2,0) = PAWPx*xyz2D(i      ,1,0)
     &                         + B10(i)
*
            PAWPy = xyz2D(i+iOffy,1,0)
            xyz2D(i+iOffy,2,0) = PAWPy*xyz2D(i+iOffy,1,0)
     &                         + B10(i)
*
            PAWPz = xyz2D(i+iOffz,1,0)
            xyz2D(i+iOffz,1,0) = PAWPz*xyz2D(i+iOffz,0,0)
            xyz2D(i+iOffz,2,0) = PAWPz*xyz2D(i+iOffz,1,0)
     &                         + B10(i)*xyz2D(i+iOffz,0,0)
 2001    Continue
         If (nabMax.gt.2) Then
            Fact = Two
            Do 250 iab = 2, nabMax-1
               Do 260 i = 1, nArg*lRys
                  PAWPx = xyz2D(i      ,1,0)
                  PAWPy = xyz2D(i+iOffy,1,0)
                  PAWPz = xyz2D(i      ,0,0)
                  temp1x= PAWPx * xyz2D(i      ,iab,0)
                  temp1y= PAWPy * xyz2D(i+iOffy,iab,0)
                  temp1z= PAWPz * xyz2D(i+iOffz,iab,0)
                  temp2x= Fact * B10(i) * xyz2D(i      ,iab-1,0)
                  temp2y= Fact * B10(i) * xyz2D(i+iOffy,iab-1,0)
                  temp2z= Fact * B10(i) * xyz2D(i+iOffz,iab-1,0)
                  xyz2D(i      ,iab+1,0) = temp1x + temp2x
                  xyz2D(i+iOffy,iab+1,0) = temp1y + temp2y
                  xyz2D(i+iOffz,iab+1,0) = temp1z + temp2z
 260           Continue
               Fact = Fact + One
 250        Continue
         End If
      End If
*
*     Compute 2D integrals with index (0,i)
*
      If (ncdMax.eq.1) Then
         Do 2002 i = 1, nArg*lRys
*
            QCWQz = xyz2D(i+iOffz,0,1)
            xyz2D(i+iOffz,0,1) = QCWQz*xyz2D(i+iOffz,0,0)
 2002    Continue
      Else If (ncdMax.ge.2) Then
         Do 2003 i = 1, nArg*lRys
            QCWQx = xyz2D(i,0,1)
            xyz2D(i      ,0,2) = QCWQx*xyz2D(i      ,0,1)
     &                         + B01(i)
*
            QCWQy = xyz2D(i+iOffy,0,1)
            xyz2D(i+iOffy,0,2) = QCWQy*xyz2D(i+iOffy,0,1)
     &                         + B01(i)
*
            QCWQz = xyz2D(i+iOffz,0,1)
            xyz2D(i+iOffz,0,1) = QCWQz*xyz2D(i+iOffz,0,0)
            xyz2D(i+iOffz,0,2) = QCWQz*xyz2D(i+iOffz,0,1)
     &                         + B01(i)*xyz2D(i+iOffz,0,0)
 2003    Continue
         If (ncdMax.gt.2) Then
            Fact = Two
            Do 350 icd = 2, ncdMax-1
               Do 360 i = 1, nArg*lRys
                  QCWQx = xyz2D(i      ,0,1)
                  QCWQy = xyz2D(i+iOffy,0,1)
                  QCWQz = xyz2D(i+iOffy,0,0)
                  temp1x= QCWQx * xyz2D(i      ,0,icd)
                  temp1y= QCWQy * xyz2D(i+iOffy,0,icd)
                  temp1z= QCWQz * xyz2D(i+iOffz,0,icd)
                  temp2x= Fact * B01(i) * xyz2D(i      ,0,icd-1)
                  temp2y= Fact * B01(i) * xyz2D(i+iOffy,0,icd-1)
                  temp2z= Fact * B01(i) * xyz2D(i+iOffz,0,icd-1)
                  xyz2D(i      ,0,icd+1) = temp1x+ temp2x
                  xyz2D(i+iOffy,0,icd+1) = temp1y+ temp2y
                  xyz2D(i+iOffz,0,icd+1) = temp1z+ temp2z
 360           Continue
               Fact = Fact + One
 350        Continue
         End If
      End If
*
*
*     Compute 2D integrals with index (i,j)
*
      If (ncdMax.le.nabMax) Then
         Fac1 = One
         Do 400 icd = 1, ncdMax
            If (icd.eq.1) Then
               Do 425 i = 1, nArg*lRys
                  PAWPx = xyz2D(i      ,1,0)
                  PAWPy = xyz2D(i+iOffy,1,0)
                  PAWPz = xyz2D(i      ,0,0)
                  xyz2D(i      ,1,1) = PAWPx*xyz2D(i      ,0,1)
     &                               +        B00(i)
                  xyz2D(i+iOffy,1,1) = PAWPy*xyz2D(i+iOffy,0,1)
     &                               +        B00(i)
                  xyz2D(i+iOffz,1,1) = PAWPz*xyz2D(i+iOffz,0,1)
     &                               +        B00(i)*xyz2D(i+iOffz,0,0)
 425           Continue
            Else
               Do 4251 i = 1, nArg*lRys
                  PAWPx = xyz2D(i      ,1,0)
                  PAWPy = xyz2D(i+iOffy,1,0)
                  PAWPz = xyz2D(i,      0,0)
                  xyz2D(i      ,1,icd) = PAWPx * xyz2D(i,0,icd)
     &                           + Fac1 * B00(i) * xyz2D(i,0,icd-1)
                  xyz2D(i+iOffy,1,icd) = PAWPy * xyz2D(i+iOffy,0,icd)
     &                           + Fac1* B00(i) * xyz2D(i+iOffy,0,icd-1)
                  xyz2D(i+iOffz,1,icd)  = PAWPz * xyz2D(i+iOffz,0,icd)
     &                           + Fac1* B00(i) * xyz2D(i+iOffz,0,icd-1)
 4251          Continue
            End If
            If (nabMax.eq.2) Then
               Do 420 i = 1, nArg*lRys
                  PAWPx = xyz2D(i      ,1,0)
                  PAWPy = xyz2D(i+iOffy,1,0)
                  PAWPz = xyz2D(i      ,0,0)
                  xyz2D(i,2,icd) =
     &                               PAWPx * xyz2D(i,1,icd)
     &                           +  B10(i) * xyz2D(i,0,icd)
     &                      + Fac1 *B00(i) * xyz2D(i,1,icd-1)
                  xyz2D(i+iOffy,2,icd) =
     &                               PAWPy * xyz2D(i+iOffy,1,icd)
     &                           +  B10(i) * xyz2D(i+iOffy,0,icd)
     &                      + Fac1 *B00(i) * xyz2D(i+iOffy,1,icd-1)
                  xyz2D(i+iOffz,2,icd) =
     &                               PAWPz * xyz2D(i+iOffz,1,icd)
     &                           +  B10(i) * xyz2D(i+iOffz,0,icd)
     &                      + Fac1 *B00(i) * xyz2D(i+iOffz,1,icd-1)
 420              Continue
            Else If (nabMax.gt.2) Then
               Fac2 = One
               Do 450 iab = 1, nabMax-1
                  Do 460 i = 1, nArg*lRys
                     PAWPx = xyz2D(i      ,1,0)
                     PAWPy = xyz2D(i+iOffy,1,0)
                     PAWPz = xyz2D(i      ,0,0)
                     temp1x= PAWPx * xyz2D(i     ,iab,icd)
                     temp1y= PAWPy * xyz2D(i+iOffy,iab,icd)
                     temp1z= PAWPz * xyz2D(i+iOffz,iab,icd)
                     temp2x= Fac2 *B10(i) *xyz2D(i      ,iab-1,icd)
                     temp2y= Fac2 *B10(i) *xyz2D(i+iOffy,iab-1,icd)
                     temp2z= Fac2 *B10(i) *xyz2D(i+iOffz,iab-1,icd)
                     temp3x= Fac1 *B00(i) *xyz2D(i      ,iab,icd-1)
                     temp3y= Fac1 *B00(i) *xyz2D(i+iOffy,iab,icd-1)
                     temp3z= Fac1 *B00(i) *xyz2D(i+iOffz,iab,icd-1)
                     xyz2D(i      ,iab+1,icd) = temp1x+ temp2x+ temp3x
                     xyz2D(i+iOffy,iab+1,icd) = temp1y+ temp2y+ temp3y
                     xyz2D(i+iOffz,iab+1,icd) = temp1z+ temp2z+ temp3z
 460              Continue
                  Fac2 = Fac2 + One
 450           Continue
            End If
            Fac1 = Fac1 + One
 400     Continue
      Else
         Fac1 = One
         Do 500 iab = 1, nabMax
            If (iab.eq.1) Then
               Do 525 i = 1, nArg*lRys
                  QCWQx = xyz2D(i      ,0,1)
                  QCWQy = xyz2D(i+iOffy,0,1)
                  QCWQz = xyz2D(i+iOffy,0,0)
                  xyz2D(i      ,1,1) = QCWQx*xyz2D(i      ,1,0)
     &                               +        B00(i)
                  xyz2D(i+iOffy,1,1) = QCWQy*xyz2D(i+iOffy,1,0)
     &                               +        B00(i)
                  xyz2D(i+iOffz,1,1) = QCWQz*xyz2D(i+iOffz,1,0)
     &                               +        B00(i)*xyz2D(i+iOffz,0,0)
 525           Continue
            Else
               Do 5251 i = 1, nArg*lRys
                  QCWQx = xyz2D(i      ,0,1)
                  QCWQy = xyz2D(i+iOffy,0,1)
                  QCWQz = xyz2D(i+iOffy,0,0)
                  xyz2D(i,iab,1) = QCWQx *xyz2D(i,iab,0)
     &                           + Fac1 *B00(i) *xyz2D(i,iab-1,0)
                  xyz2D(i+iOffy,iab,1) =
     &                            QCWQy *xyz2D(i+iOffy,iab,0)
     &                           + Fac1 *B00(i) *xyz2D(i+iOffy,iab-1,0)
                  xyz2D(i+iOffz,iab,1) =
     &                            QCWQz *xyz2D(i+iOffz,iab,0)
     &                           + Fac1 *B00(i) *xyz2D(i+iOffz,iab-1,0)
 5251          Continue
            End If
            If (ncdMax.eq.2) Then
               Do 520 i = 1, nArg*lRys
                  QCWQx = xyz2D(i      ,0,1)
                  QCWQy = xyz2D(i+iOffy,0,1)
                  QCWQz = xyz2D(i+iOffy,0,0)
                  xyz2D(i,iab,2) =QCWQx   *xyz2D(i,iab,1)
     &                           + B01(i) *xyz2D(i,iab,0)
     &                     + Fac1 *B00(i) *xyz2D(i,iab-1,1)
                  xyz2D(i+iOffy,iab,2) =
     &                            QCWQy   *xyz2D(i+iOffy,iab,1)
     &                           + B01(i) *xyz2D(i+iOffy,iab,0)
     &                     + Fac1 *B00(i) *xyz2D(i+iOffy,iab-1,1)
                  xyz2D(i+iOffz,iab,2) =
     &                            QCWQz   *xyz2D(i+iOffz,iab,1)
     &                           + B01(i) *xyz2D(i+iOffz,iab,0)
     &                     + Fac1 *B00(i) *xyz2D(i+iOffz,iab-1,1)
 520           Continue
            Else If (ncdMax.gt.2) Then
               Fac2 = One
               Do 550 icd = 1, ncdmax-1
                  Do 560 i = 1, nArg*lRys
                     QCWQx = xyz2D(i      ,0,1)
                     QCWQy = xyz2D(i+iOffy,0,1)
                     QCWQz = xyz2D(i+iOffy,0,0)
                     temp1x= QCWQx *xyz2D(i      ,iab,icd)
                     temp1y= QCWQy *xyz2D(i+iOffy,iab,icd)
                     temp1z= QCWQz *xyz2D(i+iOffz,iab,icd)
                     temp2x= Fac2 *B01(i) *xyz2D(i      ,iab,icd-1)
                     temp2y= Fac2 *B01(i) *xyz2D(i+iOffy,iab,icd-1)
                     temp2z= Fac2 *B01(i) *xyz2D(i+iOffz,iab,icd-1)
                     temp3x= Fac1 *B00(i) *xyz2D(i      ,iab-1,icd)
                     temp3y= Fac1 *B00(i) *xyz2D(i+iOffy,iab-1,icd)
                     temp3z= Fac1 *B00(i) *xyz2D(i+iOffz,iab-1,icd)
                     xyz2D(i      ,iab,icd+1) = temp1x+ temp2x+ temp3x
                     xyz2D(i+iOffy,iab,icd+1) = temp1y+ temp2y+ temp3y
                     xyz2D(i+iOffz,iab,icd+1) = temp1z+ temp2z+ temp3z
 560              Continue
                  Fac2 = Fac2 + One
 550           Continue
            End If
            Fac1 = Fac1 + One
 500     Continue
      End If
      End If
*
#ifdef _DEBUGPRINT_
      Do iab = 0, nabMax
         Do icd = 0, ncdMax
            Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
            Call RecPrt(Label,' ',
     &                  xyz2D(1            ,iab,icd),lRys,nArg)
            Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
            Call RecPrt(Label,' ',
     &                  xyz2D(1+  nArg*lRys,iab,icd),lRys,nArg)
            Write (Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
            Call RecPrt(Label,' ',
     &                  xyz2D(1+2*nArg*lRys,iab,icd),lRys,nArg)
         End Do
      End Do
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
