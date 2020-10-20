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
      SubRoutine SchInt_mck(CoorM,iAnga,iCmp,
     &                  nAlpha,nBeta,nMemab,
     &                  Zeta,ZInv,rKapab,P,nZeta,
     &                  Work2,nWork2,Work3,nWork3)
************************************************************************
*                                                                      *
* Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
*         prescreening. This is done for all unique pairs of centers   *
*         generated from the symmetry unique centers A and B.          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             June '91, modified to compute zeta, P, kappa and inte-   *
*             grals for Schwartz inequality in a k2 loop.              *
*             April '92 modified from k2Loop to a separate subroutine  *
*              for estimates of the gradient.                          *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
*     External TERISq, ModU2, Cff2Dq
      External TERIS, ModU2, Cff2DS,rys2d
#include "real.fh"
      Real*8  CoorM(3,4), CoorAC(3,2),
     &       Zeta(nZeta), ZInv(nZeta), rKapab(nZeta), P(nZeta,3),
     &       Q(3),  Work2(nWork2), Work3(nWork3)
      Integer iAnga(4), iCmp(4)
      Logical EQ, NoSpecial
*
*     Statement function to compute canonical index
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      nElem(i)    = (i+1)*(i+2)/2
*
      call dcopy_(3,[One],0,Q,1)
      la = iAnga(1)
      lb = iAnga(2)
      iCmpa = iCmp(1)
      jCmpb = iCmp(2)
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
      mcdMin=mabmin
      mcdMax=mabMax
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
*
*-----Find the proper centers to start of with the angular
*     momentum on. If la.eq.lb there will excist an
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
      mZeta = nAlpha*nBeta
*
*-----Compute [a0|c0], ijkl,a,c
*
      nT = mZeta*1
      NoSpecial=.True.
      Call Rys(iAnga,nT,
     &               Zeta,ZInv,              mZeta,[One],[One],1,
     &               P,         nZeta,Q,1,rKapab,      [One],
     &               CoorM,CoorM,CoorAC,
     &               mabMin,mabMax,mcdMin,mcdMax,
     &               Work2,nWork2,TERIS,ModU2,Cff2DS,
     &               Rys2D,NoSpecial)
*
*-----Apply transfer equation to generate [a0|cd], IJKLa,c,d
*
      nijkla = (mabMax-mabMin+1)*mZeta
      Call HRR(la,lb,CoorM(1,1),CoorM(1,2),Work2,nijkla,nMemab,ipIn)
*
*-----Transform to spherical gaussians [a0|CD], CDIJKL,a. This
*     will also put the integrals in the right position for the
*     transfer equation.
*
      Call CrSph1(Work2(ipIn),nijkla,
     &            Work3,nWork3,
     &            RSph(ipSph(la)),nElem(la),nElem(la),
     &            .False.,.False.,
     &            RSph(ipSph(lb)),nElem(lb),nElem(lb),
     &            .False.,.False.,
     &            Work2,nElem(la)*nElem(lb))
*
*-----Apply transfer equation to generate [ab|CD], CDIJKL,a,b
*
      ijklcd = nElem(la)*nElem(lb)*mZeta
      Call HRR(la,lb,CoorM(1,1),CoorM(1,2),Work2,ijklcd,nMemab,ipIn)
*
*-----Transform to spherical gaussians [AB|CD], IJKL,ABCD
*
      Call CrSph2(Work2(ipIn),mZeta,
     &            nElem(la)*nElem(lb),Work3,nWork3,
     &            RSph(ipSph(la)),nElem(la),nElem(la),
     &            .False.,.False.,
     &            RSph(ipSph(lb)),nElem(lb),nElem(lb),
     &            .False.,.False.,
     &            Work2,nElem(la)*nElem(lb))
*
      Return
      End
