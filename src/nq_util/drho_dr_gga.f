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
* Copyright (C) 2002, Roland Lindh                                     *
************************************************************************
      Subroutine dRho_dR_GGA(Dens,nDens,nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                       nSym,nGrad_Eff,list_g,Maps2p,nShell,
     &                       Grid_Type,Fixed_Grid,Fact,ndc,TabAOMax,T_X,
     &                       list_bas,Index,nIndex)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: Do_Batch                                                *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2002 (at Todai)                        *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer On, Off
      Parameter (On=1, Off=0)
      Integer list_s(2,nlist_s), list_g(3,nlist_s),
     &        ipTabAO(nlist_s), Maps2p(nShell,0:nSym-1),
     &        list_bas(2,nlist_s), Index(nIndex)
      Integer Grid_Type, Fixed_Grid
      Real*8 Dens(nDens,nD), dRho_dR(ndRho_dR,mGrid,nGrad_Eff),
     &       TabAO(nTabAO), Fact(ndc**2), TabAOMax(nlist_s)
      Integer IndGrd_Eff(3,2)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call QEnter('dRho_dR_GGA')
      If (Debug) Then
         Call RecPrt('dRho_dR_GGA:Dens',' ',Dens,nDens,nD)
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nIrrep=',nIrrep
         Write (6,*) 'nlist_s=',nlist_s
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iS = list_s(1,ilist_s)
            iCmp  = iSD( 2,iS)
            iBas  = iSD( 3,iS)
            iBas_Eff = list_bas(1,ilist_s)
            mTabAO=iBas_Eff*iCmp
            mdci  = iSD(10,iS)
            Call RecPrt('dRho_dR_GGA: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      Call FZero(dRho_dR,ndRho_dR*mGrid*nGrad_Eff)
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iS    = list_s(1,ilist_s)
         kDCRE = list_s(2,ilist_s)
         iCmp  = iSD( 2,iS)
         iBas_Eff=list_bas(1,ilist_s)
         ix = iDAMax_(mAO*mGrid*iBas_Eff*iCmp,TabAO(ipTabAO(iList_s)),1)
         TabAOMax(ilist_s)=Abs(TabAO(ipTabAO(ilist_s)-1+ix))
         TMax_i=TabAOMax(ilist_s)
         If (TMax_i.le.T_X) Go To 999
         iBas  = iSD( 3,iS)
         iPrim = iSD( 5,iS)
         mdci  = iSD(10,iS)
         iShell= iSD(11,iS)
         index_i=list_bas(2,ilist_s)
*
         lDCRE=NrOpr(kDCRE)
*
         jNQ=Maps2p(iS,lDCRE)
         Call ICopy(3,list_g(1,ilist_s),1,IndGrd_Eff(1,1),1)
         n1 = IndGrd_Eff(1,1) + IndGrd_Eff(2,1) + IndGrd_Eff(3,1)
*
         Do jlist_s=1,ilist_s
            TMax_j=TabAOMax(jlist_s)
            If (TMax_j.le.T_X) Go To 98
            If (TMax_i*TMax_j.lt.T_X) Go To 98
            jS    = list_s(1,jlist_s)
            kDCRR = list_s(2,jlist_s)
            jCmp  = iSD( 2,jS)
            jBas  = iSD( 3,jS)
            jPrim = iSD( 5,jS)
            mdcj  = iSD(10,jS)
            jShell= iSD(11,jS)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
*
            lDCRR=NrOpr(kDCRR)
*
            kNQ=Maps2p(jS,lDCRR)
            Call ICopy(3,list_g(1,jlist_s),1,IndGrd_Eff(1,2),1)
            n2 = IndGrd_Eff(1,2) + IndGrd_Eff(2,2) + IndGrd_Eff(3,2)
            If (n1+n2.eq.0) Go To 98
*
            mDij=iBas*jBas*iCmp*jCmp
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ipTmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,nD)
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iS,jS
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*ndc + mdci
*
            Deg=Two
            If (ilist_s.eq.jlist_s) Deg=One
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)
*
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
            If (nD.ne.1) Then
               ix=iDAMax_(mDij,DeDe(ip_D_a),1)
               iy=iDAMax_(mDij,DeDe(ip_D_b),1)
               DMax_ij=Half*( Abs(DeDe(ip_D_a-1+ix))
     &                       +Abs(DeDe(ip_D_b-1+iy)) )
            Else
               ix=iDAMax_(mDij,DeDe(ip_D_a),1)
               DMax_ij=Abs(DeDe(ip_D_a-1+ix))
            End If
            DMax_ij=DMax_ij*Fact(ij)
            If (TMax_i*TMax_j*DMax_ij.lt.T_X) Go To 98
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*) 'dRho_dR_GGA'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
            End If
#endif
*
            If (iShell.lt.jShell) Then
               itmp=IndGrd_Eff(1,1)
               IndGrd_Eff(1,1)=IndGrd_Eff(1,2)
               IndGrd_Eff(1,2)=itmp
               itmp=IndGrd_Eff(2,1)
               IndGrd_Eff(2,1)=IndGrd_Eff(2,2)
               IndGrd_Eff(2,2)=itmp
               itmp=IndGrd_Eff(3,1)
               IndGrd_Eff(3,1)=IndGrd_Eff(3,2)
               IndGrd_Eff(3,2)=itmp
            End If
            If (nD.eq.1) Then
               If (iShell.ge.jShell) Then
               Call Do_Rho8da(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho8da(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho8d_(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho8d_(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            End If
*
 98         Continue
         End Do                      ! jlist_s
 999     Continue
      End Do                         ! ilist_s
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('dRho_dR_GGA: dRho_dR',' ',dRho_dR,
     &                        ndRho_dR*mGrid,nGrad_Eff)
*
      Call QExit('dRho_dR_GGA')
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Dens)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(Grid_Type)
         Call Unused_integer(Fixed_Grid)
      End If
      End
      Subroutine Do_Rho8da(dRho_dR,     mGrid,nGrad_Eff,
     &                     DAij,          mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(  4,mGrid,nGrad_Eff), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
*     Rather than to use integers directly we set up a symbolic way of
*     doing the indexation.
*
      Integer F,
     &        dFdx, dFdy, dFdz,
     &        d2Fdx2, d2Fdxdy, d2Fdxdz, d2Fdy2, d2Fdydz, d2Fdz2,
     &                d2Fdydx, d2Fdzdx,         d2Fdzdy
      Parameter(F=1,
     &          dFdx=2, dFdy=3, dFdz=4,
     &          d2Fdx2=5, d2Fdxdy=6, d2Fdxdz=7, d2Fdy2=8, d2Fdydz=9,
     &                    d2Fdydx=6, d2Fdzdx=7,           d2Fdzdy=9,
     &          d2Fdz2=10)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            If (TMax_ij*Abs(DAij_).lt.T_X) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*---------- x component
*
            Ind_x1=IndGrd_Eff(1,1)
            Ind_x2=IndGrd_Eff(1,2)
            If (Ind_x1.ne.0.and.Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x1 =TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_xx1=TabAO1(d2Fdx2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy1=TabAO1(d2Fdxdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_xz1=TabAO1(d2Fdxdz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_x2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdx2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                                   +Rho_x1 *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                                   +Rho_xx1*DAij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                                   +Rho_xy1*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                                   +Rho_xz1*DAij_
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                                   +Rho_x2 *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                                   +Rho_xx2*DAij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                                   +Rho_xy2*DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                                   +Rho_xz2*DAij_
               End Do ! iGrid
            Else If (Ind_x1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x1 =TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_xx1=TabAO1(d2Fdx2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy1=TabAO1(d2Fdxdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_xz1=TabAO1(d2Fdxdz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                                   +Rho_x1 *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                                   +Rho_xx1*DAij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                                   +Rho_xy1*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                                   +Rho_xz1*DAij_
               End Do ! iGrid
            Else If (Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdx2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                                   +Rho_x2 *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                                   +Rho_xx2*DAij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                                   +Rho_xy2*DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                                   +Rho_xz2*DAij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*---------- y component
*
            Ind_y1=IndGrd_Eff(2,1)
            Ind_y2=IndGrd_Eff(2,2)
            If (Ind_y1.ne.0.and.Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y1 =TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_yx1=TabAO1(d2Fdydx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_yy1=TabAO1(d2Fdy2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz1=TabAO1(d2Fdydz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_y2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdy2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                                   +Rho_y1 *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                                   +Rho_yx1*DAij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                                   +Rho_yy1*DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                                   +Rho_yz1*DAij_
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                                   +Rho_y2 *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                                   +Rho_yx2*DAij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                                   +Rho_yy2*DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                                   +Rho_yz2*DAij_
               End Do ! iGrid
            Else If (Ind_y1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y1 =TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_yx1=TabAO1(d2Fdydx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_yy1=TabAO1(d2Fdy2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz1=TabAO1(d2Fdydz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                                   +Rho_y1 *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                                   +Rho_yx1*DAij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                                   +Rho_yy1*DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                                   +Rho_yz1*DAij_
               End Do ! iGrid
            Else If (Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdy2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                                   +Rho_y2 *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                                   +Rho_yx2*DAij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                                   +Rho_yy2*DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                                   +Rho_yz2*DAij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*---------- z component
*
            Ind_z1=IndGrd_Eff(3,1)
            Ind_z2=IndGrd_Eff(3,2)
            If (Ind_z1.ne.0.and.Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z1 =TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_zx1=TabAO1(d2Fdzdx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_zy1=TabAO1(d2Fdzdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_zz1=TabAO1(d2Fdz2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_z2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdz2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                                   +Rho_z1 *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                                   +Rho_zx1*DAij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                                   +Rho_zy1*DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                                   +Rho_zz1*DAij_
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                                   +Rho_z2 *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                                   +Rho_zx2*DAij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                                   +Rho_zy2*DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                                   +Rho_zz2*DAij_
               End Do ! iGrid
            Else If (Ind_z1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z1 =TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_zx1=TabAO1(d2Fdzdx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_zy1=TabAO1(d2Fdzdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_zz1=TabAO1(d2Fdz2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                                   +Rho_z1 *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                                   +Rho_zx1*DAij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                                   +Rho_zy1*DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                                   +Rho_zz1*DAij_
               End Do ! iGrid
            Else If (Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdz2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                                   +Rho_z2 *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                                   +Rho_zx2*DAij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                                   +Rho_zy2*DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                                   +Rho_zz2*DAij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho8d_(dRho_dR,     mGrid,nGrad_Eff,
     &                     DAij,DBij,     mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   8,mGrid,nGrad_Eff),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
*     Rather than to use integers directly we set up a symbolic way of
*     doing the indexation.
*
      Integer F,
     &        dFdx, dFdy, dFdz,
     &        d2Fdx2, d2Fdxdy, d2Fdxdz, d2Fdy2, d2Fdydz, d2Fdz2,
     &                d2Fdydx, d2Fdzdx,         d2Fdzdy
      Parameter(F=1,
     &          dFdx=2, dFdy=3, dFdz=4,
     &          d2Fdx2=5, d2Fdxdy=6, d2Fdxdz=7, d2Fdy2=8, d2Fdydz=9,
     &                    d2Fdydx=6, d2Fdzdx=7,           d2Fdzdy=9,
     &          d2Fdz2=10)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            DBij_=DBij(iCB,jCB)*Fact
            Dij_ =Half*(Abs(DAij_)+Abs(DBij_))
            If (TMax_ij*Abs(Dij_).lt.T_X) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*---------- x component
*
            Ind_x1=IndGrd_Eff(1,1)
            Ind_x2=IndGrd_Eff(1,2)
            If (Ind_x1.ne.0.and.Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x1 =TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_xx1=TabAO1(d2Fdx2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy1=TabAO1(d2Fdxdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_xz1=TabAO1(d2Fdxdz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_x2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdx2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                              +Rho_x1 *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                              +Rho_x1 *DBij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                              +Rho_xx1*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                              +Rho_xy1*DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                              +Rho_xz1*DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                              +Rho_xx1*DBij_
                  dRho_dR(7,iGrid,Ind_x1)=dRho_dR(7,iGrid,Ind_x1)
     &                              +Rho_xy1*DBij_
                  dRho_dR(8,iGrid,Ind_x1)=dRho_dR(8,iGrid,Ind_x1)
     &                              +Rho_xz1*DBij_
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                              +Rho_x2 *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                              +Rho_x2 *DBij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                              +Rho_xx2*DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                              +Rho_xy2*DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                              +Rho_xz2*DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                              +Rho_xx2*DBij_
                  dRho_dR(7,iGrid,Ind_x2)=dRho_dR(7,iGrid,Ind_x2)
     &                              +Rho_xy2*DBij_
                  dRho_dR(8,iGrid,Ind_x2)=dRho_dR(8,iGrid,Ind_x2)
     &                              +Rho_xz2*DBij_
               End Do ! iGrid
            Else If (Ind_x1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x1 =TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_xx1=TabAO1(d2Fdx2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy1=TabAO1(d2Fdxdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_xz1=TabAO1(d2Fdxdz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                              +Rho_x1 *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                              +Rho_x1 *DBij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                              +Rho_xx1*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                              +Rho_xy1*DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                              +Rho_xz1*DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                              +Rho_xx1*DBij_
                  dRho_dR(7,iGrid,Ind_x1)=dRho_dR(7,iGrid,Ind_x1)
     &                              +Rho_xy1*DBij_
                  dRho_dR(8,iGrid,Ind_x1)=dRho_dR(8,iGrid,Ind_x1)
     &                              +Rho_xz1*DBij_
               End Do ! iGrid
            Else If (Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_x2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdx2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_xz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdxdz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                              +Rho_x2 *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                              +Rho_x2 *DBij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                              +Rho_xx2*DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                              +Rho_xy2*DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                              +Rho_xz2*DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                              +Rho_xx2*DBij_
                  dRho_dR(7,iGrid,Ind_x2)=dRho_dR(7,iGrid,Ind_x2)
     &                              +Rho_xy2*DBij_
                  dRho_dR(8,iGrid,Ind_x2)=dRho_dR(8,iGrid,Ind_x2)
     &                              +Rho_xz2*DBij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*---------- y component
*
            Ind_y1=IndGrd_Eff(2,1)
            Ind_y2=IndGrd_Eff(2,2)
            If (Ind_y1.ne.0.and.Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y1 =TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_yx1=TabAO1(d2Fdydx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_yy1=TabAO1(d2Fdy2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz1=TabAO1(d2Fdydz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_y2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdy2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                              +Rho_y1 *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                              +Rho_y1 *DBij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                              +Rho_yx1*DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                              +Rho_yy1*DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                              +Rho_yz1*DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                              +Rho_yx1*DBij_
                  dRho_dR(7,iGrid,Ind_y1)=dRho_dR(7,iGrid,Ind_y1)
     &                              +Rho_yy1*DBij_
                  dRho_dR(8,iGrid,Ind_y1)=dRho_dR(8,iGrid,Ind_y1)
     &                              +Rho_yz1*DBij_
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                              +Rho_y2 *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                              +Rho_y2 *DBij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                              +Rho_yx2*DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                              +Rho_yy2*DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                              +Rho_yz2*DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                              +Rho_yx2*DBij_
                  dRho_dR(7,iGrid,Ind_y2)=dRho_dR(7,iGrid,Ind_y2)
     &                              +Rho_yy2*DBij_
                  dRho_dR(8,iGrid,Ind_y2)=dRho_dR(8,iGrid,Ind_y2)
     &                              +Rho_yz2*DBij_
               End Do ! iGrid
            Else If (Ind_y1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y1 =TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_yx1=TabAO1(d2Fdydx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_yy1=TabAO1(d2Fdy2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz1=TabAO1(d2Fdydz,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                              +Rho_y1 *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                              +Rho_y1 *DBij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                              +Rho_yx1*DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                              +Rho_yy1*DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                              +Rho_yz1*DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                              +Rho_yx1*DBij_
                  dRho_dR(7,iGrid,Ind_y1)=dRho_dR(7,iGrid,Ind_y1)
     &                              +Rho_yy1*DBij_
                  dRho_dR(8,iGrid,Ind_y1)=dRho_dR(8,iGrid,Ind_y1)
     &                              +Rho_yz1*DBij_
               End Do ! iGrid
            Else If (Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_y2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdy2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_yz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdydz,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                              +Rho_y2 *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                              +Rho_y2 *DBij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                              +Rho_yx2*DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                              +Rho_yy2*DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                              +Rho_yz2*DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                              +Rho_yx2*DBij_
                  dRho_dR(7,iGrid,Ind_y2)=dRho_dR(7,iGrid,Ind_y2)
     &                              +Rho_yy2*DBij_
                  dRho_dR(8,iGrid,Ind_y2)=dRho_dR(8,iGrid,Ind_y2)
     &                              +Rho_yz2*DBij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*---------- z component
*
            Ind_z1=IndGrd_Eff(3,1)
            Ind_z2=IndGrd_Eff(3,2)
            If (Ind_z1.ne.0.and.Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z1 =TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_zx1=TabAO1(d2Fdzdx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_zy1=TabAO1(d2Fdzdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_zz1=TabAO1(d2Fdz2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_z2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdz2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                              +Rho_z1 *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                              +Rho_z1 *DBij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                              +Rho_zx1*DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                              +Rho_zy1*DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                              +Rho_zz1*DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                              +Rho_zx1*DBij_
                  dRho_dR(7,iGrid,Ind_z1)=dRho_dR(7,iGrid,Ind_z1)
     &                              +Rho_zy1*DBij_
                  dRho_dR(8,iGrid,Ind_z1)=dRho_dR(8,iGrid,Ind_z1)
     &                              +Rho_zz1*DBij_
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                              +Rho_z2 *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                              +Rho_z2 *DBij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                              +Rho_zx2*DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                              +Rho_zy2*DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                              +Rho_zz2*DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                              +Rho_zx2*DBij_
                  dRho_dR(7,iGrid,Ind_z2)=dRho_dR(7,iGrid,Ind_z2)
     &                              +Rho_zy2*DBij_
                  dRho_dR(8,iGrid,Ind_z2)=dRho_dR(8,iGrid,Ind_z2)
     &                              +Rho_zz2*DBij_
               End Do ! iGrid
            Else If (Ind_z1.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z1 =TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
                  Rho_zx1=TabAO1(d2Fdzdx,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdx   ,iGrid,jCB_Eff)
                  Rho_zy1=TabAO1(d2Fdzdy,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdy   ,iGrid,jCB_Eff)
                  Rho_zz1=TabAO1(d2Fdz2 ,iGrid,iCB_Eff)*
     &                    TabAO2(F      ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                              +Rho_z1 *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                              +Rho_z1 *DBij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                              +Rho_zx1*DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                              +Rho_zy1*DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                              +Rho_zz1*DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                              +Rho_zx1*DBij_
                  dRho_dR(7,iGrid,Ind_z1)=dRho_dR(7,iGrid,Ind_z1)
     &                              +Rho_zy1*DBij_
                  dRho_dR(8,iGrid,Ind_z1)=dRho_dR(8,iGrid,Ind_z1)
     &                              +Rho_zz1*DBij_
               End Do ! iGrid
            Else If (Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  Rho_z2 =TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zx2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdx,iGrid,jCB_Eff)
     &                   +TabAO1(dFdx   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zy2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdzdy,iGrid,jCB_Eff)
     &                   +TabAO1(dFdy   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  Rho_zz2=TabAO1(F      ,iGrid,iCB_Eff)*
     &                    TabAO2(d2Fdz2 ,iGrid,jCB_Eff)
     &                   +TabAO1(dFdz   ,iGrid,iCB_Eff)*
     &                    TabAO2(dFdz   ,iGrid,jCB_Eff)
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                              +Rho_z2 *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                              +Rho_z2 *DBij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                              +Rho_zx2*DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                              +Rho_zy2*DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                              +Rho_zz2*DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                              +Rho_zx2*DBij_
                  dRho_dR(7,iGrid,Ind_z2)=dRho_dR(7,iGrid,Ind_z2)
     &                              +Rho_zy2*DBij_
                  dRho_dR(8,iGrid,Ind_z2)=dRho_dR(8,iGrid,Ind_z2)
     &                              +Rho_zz2*DBij_
               End Do ! iGrid
            End If
*                                                                      *
************************************************************************
*                                                                      *
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
