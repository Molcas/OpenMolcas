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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine Distg2(g2,Hess,nHess,IndGrd,
     &                  IfHss,IndHss,iuvwx,kOp,nop,Tr,IfGr)
************************************************************************
*                                                                      *
* @parameter kop   operators for center generator                      *
*                                                                      *
* Object: trace the gradient of the ERI's with the second order        *
*         density matrix                                               *
*                                                                      *
*     Author: Anders Bernhardsson Dept. of Theoretical Chemistry,      *
*             University of Lund, SWEDEN                               *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl, iChBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 g2(78),  Prmt(0:7),Hess(nHess)
      Logical IfHss(4,3,4,3),Tr(4),IfGr(4)
      Integer IndGrd(3,4,0:(nIrrep-1)),kOp(4),iuvwx(4),
     &          IndHss(4,3,4,3,0:(nIrrep-1)),nop(4)
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
      ix(icn,icar,jcn,jcar)=(((icn-1)*3+icar)*((icn-1)*3+icar-1))/2+
     &                       (jcn-1)*3+jcar
*                                                                      *
************************************************************************
*                                                                      *
c     iRout = 239
c     iPrint = nPrint(iRout)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call recprt('Distg2: g2(raw) ',' ',g2,1,78)
      Call recprt('Distg2: Hess(raw) ',' ',Hess,1,nHess)
#endif
*

*-----Compute some of the contributions via the translational invariance
*
      Do 200 iCn = 1, 4
         Do 210 iCar = 1, 3
           Do 220 jCn = 1,iCn
             If (iCn.eq.jCn) Then
              iStop=iCar
             Else
              iStop=3
             End If
             Do 230 jCar = 1, iStop
               if (Tr(iCn).or.Tr(jCn)) Then
                ij = ix(iCn,iCar,jCn,jCar)
                g2(ij)=zero
*------------------------------------------------------------*
*
*               Both derivatives by translation!
*
*------------------------------------------------------------*
                if (tr(iCn).and.tr(jCn))
     &          Then
                 Do 240 kCn = 1, 4
                  Do 250 lCn = 1,kCn
                     If (lCn.eq.kCn) then
c                      iMax=Max(iCar,jCar)
                      iCa2=Min(iCar,jCar)
                      iCa1=Max(iCar,jCar)
                      If (IfHss(kCn,iCa1,lCn,iCa2)) Then
                        k1=Ix(kCn,iCa1,lCn,iCa2)
                        g2(ij)=g2(ij)+g2(k1)
                      End If
                     Else
                      If (IfHss(kCn,iCar,lCn,jCar)) Then
                        k1=Ix(kCn,iCar,lCn,jCar)
                        k2=Ix(kCn,jCar,lCn,iCar)
                        g2(ij)=g2(ij)+g2(k1)+
     &                                    g2(k2)
                      End If
                     End If
 250              Continue
 240             Continue
*------------------------------------------------------------*
*
*               Centre jCn by translation
*
*------------------------------------------------------------*
                Else If (ifgr(iCn).and.tr(jCn))
     &          Then
                  Do 260 kCn = 1, 4
                    If (kCn.gt.iCn) Then
                      iCn1=kCn
                      iCn2=iCn
                      iCa1=jCar
                      iCa2=iCar
                     Else If (kCn.lt.iCn) Then
                      iCn1=iCn
                      iCn2=kCn
                      iCa1=iCar
                      iCa2=jCar
                     Else
                      iCn1=iCn
                      iCn2=kCn
                      iCa1=Max(iCar,jCar)
                      iCa2=Min(iCar,jCar)
                     End If
                    If (IfHss(iCn1,iCa1,iCn2,iCa2))
     &              Then
                      kl=Ix(iCn1,iCa1,iCn2,iCa2)
                      g2(ij)=g2(ij)-g2(kl)
                    End If
 260              Continue
*------------------------------------------------------------*
*
*               Centre iCn by translation
*
*------------------------------------------------------------*
                Else If (IfGr(jCn).and.tr(iCn))
     &          Then
                  Do 270 kCn = 1, 4
                    If (kCn.gt.jCn) Then
                      iCn1=kCn
                      iCn2=jCn
                      iCa1=iCar
                      iCa2=jCar
                     Else If (kCn.lt.jCn) Then
                      iCn1=jCn
                      iCn2=kCn
                      iCa1=jCar
                      iCa2=iCar
                     Else
                      iCn1=jCn
                      iCn2=kCn
                      iCa1=Max(iCar,jCar)
                      iCa2=Min(iCar,jCar)
                     End If
                    If (IfHss(iCn1,iCa1,iCn2,iCa2))
     &              Then
                      kl=Ix(iCn1,iCa1,iCn2,iCa2)
                      g2(ij)=g2(ij)-g2(kl)
                    End If
 270              Continue
                End If
               End If
 230        Continue
 220       Continue
 210     Continue
 200  Continue
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      Call recprt('Distg2: g2 ',' ',g2,1,78)
#endif
*-----Distribute contribution to the hessian.
*
*----------------------------------------------------------------------*
      Do 90 iIrrep=0,nIrrep-1
       Do 100 iCn = 1, 4
         Do 110 iCar = 1, 3
          Do 120 jCn =1,iCn
           if (iCn.eq.jCn) Then
               iStop=iCar
           Else
               iStop=3
           End If
           Do 130 jCar=1,istop
            If (IndHss(iCn,iCar,jCn,jCar,iIrrep).ne.0) Then
*----------------------------------------------------------------------*
*
*              Get indexes
*
*----------------------------------------------------------------------*
               ij = Ix(iCn,iCar,jCn,jCar)
               iHess = Abs(IndHss(iCn,iCar,jCn,jCar,iIrrep))
*----------------------------------------------------------------------*
*
*            Sign due to integral direction
*
*----------------------------------------------------------------------*
               ps=DBLE(iChTbl(iIrrep,nOp(iCn))*
     &                   iChTbl(iIrrep,nOp(jCn)))
*----------------------------------------------------------------------*
*
*              If over & under triangular integrals are needed
*              multiply by two instead!
*
*----------------------------------------------------------------------*
               if ((iCn.ne.jCn).and.(iCar.eq.jCar).and.
     &               (Abs(indgrd(iCar,iCn,iIrrep)).eq.
     &                Abs(indgrd(jCar,jCn,iIrrep)))) Then
                      ps=ps*Two
               End If
*----------------------------------------------------------------------*
*
*              Sign due to which symmetry group the translation is in
*
*----------------------------------------------------------------------*
               iCh1=iChBas(iCar+1)
               iCh2=iChBas(jCar+1)
               ps = ps*xPrmt(kOp(iCn),iCh1)*xPrmt(kOp(jCn),iCh2)
*----------------------------------------------------------------------*
*
*              Multiply by number of stabilisators.
*
*----------------------------------------------------------------------*
               Fact=ps*DBLE(iuvwx(iCn))/DBLE(nIrrep*nirrep)
     &                * DBLE(iuvwx(jCn))
*----------------------------------------------------------------------*
*
*               Add to hessian
*
*----------------------------------------------------------------------*
               Hess(iHess) = Hess(iHess) + Fact * g2(ij)
            End If
 130       Continue
 120      Continue
 110     Continue
 100   Continue
  90  Continue
#ifdef _DEBUGPRINT_
      Call recprt('Distg2: Hess ',' ',Hess,1,nHess)
#endif
*
      Return
      End
