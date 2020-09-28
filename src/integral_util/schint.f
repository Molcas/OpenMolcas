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
* Copyright (C) 1990-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SchInt(CoorM,
     &                  iAnga,iCmp,mZeta,Zeta,ZInv,
     &                  rKapab,P,rKapcd,Q,nZeta,Wrk,nWork2,
     &                  HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
************************************************************************
*                                                                      *
* Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
*         prescreening. This is done for all unique pairs of centers   *
*         generated from the symmetry unique centers A and B.          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN.                              *
*             June '91, modified to compute zeta, P, kappa and inte-   *
*             grals for Schwartz inequality in a k2 loop.              *
*             April '92 modified from k2Loop to a separate subroutine  *
*              for estimates of the gradient.                          *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
      External TERISq, ModU2, Cff2Dq, xRys2D
#include "real.fh"
#include "print.fh"
      Real*8  CoorM(3,4), CoorAC(3,2), HMtrx(nHrrMtrx,2),
     &       Zeta(mZeta), ZInv(mZeta), rKapab(mZeta), P(nZeta,3),
     &       Q(nZeta,3), rKapcd(mZeta), Wrk(nWork2)
      Integer iAnga(4), iCmp(4)
      Logical EQ, NoSpecial
*
*     Statement function to compute canonical index
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      nElem(i)    = (i+1)*(i+2)/2
*
      iRout = 242
      iPrint = nPrint(iRout)
*     iQ = 1
      la = iAnga(1)
      lb = iAnga(2)
      If (iPrint.ge.19) Then
         Call RecPrt(' In SchInt: CoorM',' ',CoorM,3,4)
         Call RecPrt(' In SchInt: P',' ',P,nZeta,3)
         Call RecPrt(' In SchInt: Q',' ',Q,nZeta,3)
         Write (6,*) 'iAnga=',iAnga
      End If
*
*
*-----Compute primitive integrals to be used in the prescreening
*     by the Schwartz inequality.
*
*
*-----Compute actual size of [a0|c0] block
*
      mabMin=nabSz(Max(la,lb)-1)+1
      If (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
      mabMax=nabSz(la+lb)
      mcdMin=nabSz(Max(la,lb)-1)+1
      If (EQ(CoorM(1,3),CoorM(1,4))) mcdMin = nabSz(la+lb-1)+1
      mcdMax=mabMax
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
*
*-----Find the proper centers to start of with the angular
*     momentum on. If la.eq.lb there will exist an
*     ambiguity to which center that angular momentum should
*     be accumulated on. In that case we will use A and C of
*     the order as defined by the basis functions types.
*
      If (iAnga(1).ge.iAnga(2)) Then
         call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
         call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
      Else
         call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
         call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
      End If
*
*-----Compute [a0|c0], ijkl,a,c
*
      If (iPrint.ge.19) nPrint(13)=99
      nT = mZeta*1
      NoSpecial=.True.
      Call Rys(iAnga,nT,
     &         Zeta,ZInv,mZeta,Zeta,ZInv,mZeta,P,nZeta,Q,nZeta,
     &         rKapab,rKapcd,CoorM,CoorM,CoorAC,
     &         mabMin,mabMax,mcdMin,mcdMax,
     &         Wrk,nWork2,TERISq,ModU2,Cff2Dq,
     &         xRys2D,NoSpecial)
      If (iPrint.ge.19) nPrint(13)=5
*
      If (iPrint.ge.19) Call TrcPrt(' In SchInt: ijkl,[a0|c0]',' ',
     &                              Wrk,mZeta,mabcd)
      If (iPrint.ge.59) Call RecPrt(' In SchInt: ijkl,[a0|c0]',' ',
     &                              Wrk,mZeta,mabcd)
*                                                                      *
************************************************************************
*                                                                      *
*        Generate transformation matrix from intermediate integrals
*        to final angular composition, observe that since the 2nd order
*        particle matrix is transformed from real spherical harmonics
*        to cartesians that the prescreening integrals also are in
*        cartesians!
*
         ne=(mabMax-mabMin+1)
         Call HrrMtrx(HMtrx(1,1),ne,la,lb,CoorM(1,1),CoorM(1,2),
     &                .False.,RSph(ipSph(la)),nElem(la),
     &                .False.,RSph(ipSph(lb)),nElem(lb))
         Call HrrMtrx(HMtrx(1,2),ne,la,lb,CoorM(1,3),CoorM(1,4),
     &                .False.,RSph(ipSph(la)),nElem(la),
     &                .False.,RSph(ipSph(lb)),nElem(lb))
*                                                                      *
************************************************************************
*                                                                      *
*------- Apply a transpose prior to Tnsctl to fake the action of Cntrct.
*
         iW3=1+mZeta*mabcd
         Call DGeTMO(Wrk,mZeta,mZeta,mabcd,Wrk(iW3),mabcd)
         call dcopy_(mabcd*mZeta,Wrk(iW3),1,Wrk,1)
         Call TnsCtl(Wrk,nWork2,CoorM,
     &               mabcd,mZeta,mabMax,mabMin,mabMax,mabMin,
     &               HMtrx(1,1),HMtrx(1,2),
     &               la,lb,la,lb,
     &               nElem(la),nElem(lb),nElem(la),nElem(lb),
     &               iShlla,jShllb,iShlla,jShllb,i_Int)
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.19) Call TrcPrt(' In SchInt',' ',Wrk(i_Int),
     &      mZeta,(nElem(la)*nElem(lb))**2)
      If (iPrint.ge.99) Call RecPrt(' In SchInt',' ',Wrk(i_Int),
     &      mZeta,(nElem(la)*nElem(lb))**2)
*     Call GetMem(' Exit SchInt','CHECK','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iCmp)
      End
