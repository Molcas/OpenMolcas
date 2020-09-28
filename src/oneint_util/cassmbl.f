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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine CAssmbl(Rnxyz,Axyz,la,Bxyz,lb,nZeta,HerW,nHer)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the multipole moment *
*         matrix within the framework of the Gauss-Hermite quadrature. *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 HerW(nHer)
      Complex*16 Rnxyz(nZeta*3,0:la,0:lb),
     &           Axyz(nZeta*3,nHer,0:la),
     &           Bxyz(nZeta*3,nHer,0:lb)
#ifdef _DARWIN_
      Complex*16 vv1, vv2
#endif
      Character*80 Label
*
      iRout = 123
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In CAssmbl:HerW',' ',HerW,1,nHer)
         Call CRecPrt(' In CAssmbl:Axyz',' ',
     &                Axyz,nZeta*3,nHer*(la+1),'R')
         Call CRecPrt(' In CAssmbl:Axyz',' ',
     &                Axyz,nZeta*3,nHer*(la+1),'I')
         Call CRecPrt(' In CAssmbl:Bxyz',' ',
     &                Bxyz,nZeta*3,nHer*(lb+1),'R')
         Call CRecPrt(' In CAssmbl:Bxyz',' ',
     &                Bxyz,nZeta*3,nHer*(lb+1),'I')
      End If
*
*     Initialize to zero
*
      Do ib = 0, lb
         Do ia = 0, la
            Do iZeta = 1, nZeta
               Rnxyz(iZeta,ia,ib)=DCMPLX(Zero,Zero)
            End Do
         End Do
      End Do
*
      Do 100 ia = 0, la
         Do 110 ib = 0, lb
*
*     Generate the cartesian components of the multipole moment
*     matrix as a sum of the value of the integrand, evaluated
*     at a root, times a weight.
*
               Do 30 iHer = 1, nHer
                  Do 10 iZCar = 1, 3*nZeta
c vv. splitted in order to make GNU compiler on Mac more happy.
#ifdef _DARWIN_
                  vv1=Axyz(iZCar,iHer,ia)
                  vv2=Bxyz(iZCar,iHer,ib)*HerW(iHer)
                  Rnxyz(iZCar,ia,ib)= Rnxyz(iZCar,ia,ib)+vv1*vv2
#else
                     Rnxyz(iZCar,ia,ib) = Rnxyz(iZCar,ia,ib) +
     &                             Axyz(iZCar,iHer,ia)*
     &                             Bxyz(iZCar,iHer,ib)*
     &                             HerW(iHer)
#endif
 10               Continue
 30            Continue
*
               If (iPrint.ge.99) Then
                  Write (Label,'(A,I2,A,I2,A)')
     &            ' In CAssmbl: Rnxyz(',ia,',',ib,')'
                  Call CRecPrt(Label,' ',Rnxyz(1,ia,ib),nZeta,3,'R')
                  Call CRecPrt(Label,' ',Rnxyz(1,ia,ib),nZeta,3,'I')
               End If
 110     Continue
 100  Continue
*
      Return
      End
