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
* Copyright (C) 2007, Roland Lindh                                     *
************************************************************************
      Subroutine dRho_dR_meta_GGA2
     &                      (nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                       nSym,nGrad_Eff,list_g,Maps2p,nShell,
     &                       Fact,ndc,TabAOMax,T_X,
     &                       list_bas,Index,nIndex)
************************************************************************
*                                                                      *
* Object: to compute the gradient of rho, grad rho, and nabla rho      *
*                                                                      *
* Called from: Do_Batch                                                *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Theoretical Chemistry,       *
*             Lund university, SWEDEN.  September 2007                 *
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
      Real*8 dRho_dR(ndRho_dR,mGrid,nGrad_Eff),
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
#ifdef _DEBUG_
      Call QEnter('dRho_dR_meta_GGA')
      If (Debug) Then
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
            Call RecPrt('dRho_dR_meta_GGA: TabAO',' ',
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
         lDCRE=NrOpr(kDCRE,iOper,nIrrep)
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
            lDCRR=NrOpr(kDCRR,iOper,nIrrep)
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
#ifdef _DEBUG_
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
            lDCRER=NrOpr(iER,iOper,nIrrep)
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
#ifdef _DEBUG_
            If (Debug) Then
               Write (6,*) 'dRho_dR_meta_GGA2'
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
               Call Do_Rho9da(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho9da(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho9d_(dRho_dR,    mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho9d_(dRho_dR,    mGrid,nGrad_Eff,
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
#ifdef _DEBUG_
      If (Debug) Call RecPrt('dRho_dR_meta_GGA: dRho_dR',' ',dRho_dR,
     &                        ndRho_dR*mGrid,nGrad_Eff)
*
      Call QExit('dRho_dR_meta_GGA')
#endif
      Return
      End
      Subroutine Do_Rho9da(dRho_dR,     mGrid,nGrad_Eff,
     &                     DAij,          mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(  6,mGrid,nGrad_Eff), DAij(iBas*iCmp,jBas*jCmp),
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
     &        dx, dy, dz,
     &        dx2, dxdy, dxdz, dy2, dydz, dz2,
     &             dydx, dzdx,      dzdy,
     &        dx3, dx2dy, dx2dz, dxdy2, dxdydz,
     &             dydx2, dzdx2, dy2dx, dxdzdy,
     &                                  dydxdz,
     &                                  dydzdx,
     &                                  dzdxdy,
     &                                  dzdydz,
     &        dxdz2, dy3, dy2dz, dydz2, dz3,
     &        dz2dx,      dzdy2, dz2dy
      Parameter(F=1,
     &          dx=2, dy=3, dz=4,
     &          dx2=5, dxdy=6, dxdz=7, dy2=8, dydz=9,
     &                 dydx=6, dzdx=7,        dzdy=9,
     &          dz2=10,
     &  dx3=11, dx2dy=12, dx2dz=13, dxdy2=14, dxdydz=15,
     &          dydx2=12, dzdx2=13, dy2dx=14, dxdzdy=15,
     &                                        dydxdz=15,
     &                                        dydzdx=15,
     &                                        dzdxdy=15,
     &                                        dzdydz=15,
     &     dxdz2=16, dy3=17, dy2dz=18, dydz2=19, dz3=20,
     &     dz2dx=16,         dzdy2=18, dz2dy=19
     &         )
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
                  x1   =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  x1_x =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x1_y =TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  x1_z =TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  x1_t =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  x1_r=(TabAO1(dx3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dxdy2,iGrid,iCB_Eff)
     &                 +TabAO1(dxdz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +x1_t
                  x2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
                  x2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dxdy2,iGrid,jCB_Eff)
     &                 +TabAO2(dxdz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
     &                 +x2_t
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                                   +x1  *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                                   +x1_x*DAij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                                   +x1_y*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                                   +x1_z*DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                                   +x1_t*DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                                   +x1_r*DAij_
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                                   +x2   *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                                   +x2_x *DAij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                                   +x2_y *DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                                   +x2_z *DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                                   +x2_t*DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                                   +x2_r*DAij_
               End Do ! iGrid
            Else If (Ind_x1.ne.0) Then
               Do iGrid = 1, mGrid
                  x1   =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  x1_x =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x1_y =TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  x1_z =TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  x1_t =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  x1_r=(TabAO1(dx3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dxdy2,iGrid,iCB_Eff)
     &                 +TabAO1(dxdz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +x1_t
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                                   +x1  *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                                   +x1_x*DAij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                                   +x1_y*DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                                   +x1_z*DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                                   +x1_t*DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                                   +x1_r*DAij_
               End Do ! iGrid
            Else If (Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  x2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
                  x2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dxdy2,iGrid,jCB_Eff)
     &                 +TabAO2(dxdz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
     &                 +x2_t
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                                   +x2   *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                                   +x2_x *DAij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                                   +x2_y *DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                                   +x2_z *DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                                   +x2_t*DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                                   +x2_r*DAij_
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
                  y1   =TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  y1_x =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  y1_y =TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y1_z =TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  y1_t =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  y1_r=(TabAO1(dx2dy,iGrid,iCB_Eff)
     &                 +TabAO1(dy3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dydz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +y1_t
                  y2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydx ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
                  y2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dy,iGrid,jCB_Eff)
     &                 +TabAO2(dy3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dydz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
     &                 +y2_t
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                                   +y1   *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                                   +y1_x *DAij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                                   +y1_y *DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                                   +y1_z *DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                                   +y1_t *DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                                   +y1_r *DAij_
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                                   +y2   *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                                   +y2_x *DAij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                                   +y2_y *DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                                   +y2_z *DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                                   +y2_t *DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                                   +y2_r *DAij_
               End Do ! iGrid
            Else If (Ind_y1.ne.0) Then
               Do iGrid = 1, mGrid
                  y1   =TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  y1_x =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  y1_y =TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y1_z =TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  y1_t =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  y1_r=(TabAO1(dx2dy,iGrid,iCB_Eff)
     &                 +TabAO1(dy3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dydz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +y1_t
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                                   +y1   *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                                   +y1_x *DAij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                                   +y1_y *DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                                   +y1_z *DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                                   +y1_t *DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                                   +y1_r *DAij_
               End Do ! iGrid
            Else If (Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  y2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydx ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
                  y2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dy,iGrid,jCB_Eff)
     &                 +TabAO2(dy3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dydz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
     &                 +y2_t
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                                   +y2   *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                                   +y2_x *DAij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                                   +y2_y *DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                                   +y2_z *DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                                   +y2_t *DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                                   +y2_r *DAij_
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
                  z1   =TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
                  z1_x =TabAO1(dzdx,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  z1_y =TabAO1(dzdy,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  z1_z =TabAO1(dz2 ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z1_t =TabAO1(dzdx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dzdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  z1_r=(TabAO1(dx2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dy2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dz3  ,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +z1_t
                  z2   =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_x =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdx,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_y =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdy,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2  ,iGrid,jCB_Eff)*Two
                  z2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dy2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dz3  ,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
     &                 +z2_t
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                                   +z1  *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                                   +z1_x *DAij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                                   +z1_y *DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                                   +z1_z *DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                                   +z1_t *DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                                   +z1_r *DAij_
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                                   +z2   *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                                   +z2_x *DAij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                                   +z2_y *DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                                   +z2_z *DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                                   +z2_t *DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                                   +z2_r *DAij_
               End Do ! iGrid
            Else If (Ind_z1.ne.0) Then
               Do iGrid = 1, mGrid
                  z1   =TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
                  z1_x =TabAO1(dzdx,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  z1_y =TabAO1(dzdy,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  z1_z =TabAO1(dz2 ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z1_t =TabAO1(dzdx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dzdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  z1_r=(TabAO1(dx2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dy2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dz3  ,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +z1_t
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                                   +z1  *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                                   +z1_x *DAij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                                   +z1_y *DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                                   +z1_z *DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                                   +z1_t *DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                                   +z1_r *DAij_
               End Do ! iGrid
            Else If (Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  z2   =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_x =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdx,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_y =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdy,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2  ,iGrid,jCB_Eff)*Two
                  z2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dy2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dz3  ,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
     &                 +z2_t
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                                   +z2   *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                                   +z2_x *DAij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                                   +z2_y *DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                                   +z2_z *DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                                   +z2_t *DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                                   +z2_r *DAij_
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
      Subroutine Do_Rho9d_(dRho_dR,     mGrid,nGrad_Eff,
     &                     DAij,DBij,     mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(  12,mGrid,nGrad_Eff),
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
     &        dx, dy, dz,
     &        dx2, dxdy, dxdz, dy2, dydz, dz2,
     &             dydx, dzdx,      dzdy,
     &        dx3, dx2dy, dx2dz, dxdy2, dxdydz,
     &             dydx2, dzdx2, dy2dx, dxdzdy,
     &                                  dydxdz,
     &                                  dydzdx,
     &                                  dzdxdy,
     &                                  dzdydz,
     &        dxdz2, dy3, dy2dz, dydz2, dz3,
     &        dz2dx,      dzdy2, dz2dy
      Parameter(F=1,
     &          dx=2, dy=3, dz=4,
     &          dx2=5, dxdy=6, dxdz=7, dy2=8, dydz=9,
     &                 dydx=6, dzdx=7,        dzdy=9,
     &          dz2=10,
     &  dx3=11, dx2dy=12, dx2dz=13, dxdy2=14, dxdydz=15,
     &          dydx2=12, dzdx2=13, dy2dx=14, dxdzdy=15,
     &                                        dydxdz=15,
     &                                        dydzdx=15,
     &                                        dzdxdy=15,
     &                                        dzdydz=15,
     &     dxdz2=16, dy3=17, dy2dz=18, dydz2=19, dz3=20,
     &     dz2dx=16,         dzdy2=18, dz2dy=19
     &         )
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
                  x1   =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  x1_x =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x1_y =TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  x1_z =TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  x1_t =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  x1_r=(TabAO1(dx3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dxdy2,iGrid,iCB_Eff)
     &                 +TabAO1(dxdz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +x1_t
                  x2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
                  x2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dxdy2,iGrid,jCB_Eff)
     &                 +TabAO2(dxdz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
     &                 +x2_t
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                              +x1   *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                              +x1   *DBij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                              +x1_x *DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                              +x1_y *DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                              +x1_z *DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                              +x1_x *DBij_
                  dRho_dR(7,iGrid,Ind_x1)=dRho_dR(7,iGrid,Ind_x1)
     &                              +x1_y *DBij_
                  dRho_dR(8,iGrid,Ind_x1)=dRho_dR(8,iGrid,Ind_x1)
     &                              +x1_z *DBij_
                  dRho_dR(9,iGrid,Ind_x1)=dRho_dR(9,iGrid,Ind_x1)
     &                              +x1_t *DAij_
                  dRho_dR(10,iGrid,Ind_x1)=dRho_dR(10,iGrid,Ind_x1)
     &                              +x1_t *DBij_
                  dRho_dR(11,iGrid,Ind_x1)=dRho_dR(11,iGrid,Ind_x1)
     &                              +x1_r *DAij_
                  dRho_dR(12,iGrid,Ind_x1)=dRho_dR(12,iGrid,Ind_x1)
     &                              +x1_r *DBij_
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                              +x2   *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                              +x2   *DBij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                              +x2_x *DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                              +x2_y *DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                              +x2_z *DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                              +x2_x *DBij_
                  dRho_dR(7,iGrid,Ind_x2)=dRho_dR(7,iGrid,Ind_x2)
     &                              +x2_y *DBij_
                  dRho_dR(8,iGrid,Ind_x2)=dRho_dR(8,iGrid,Ind_x2)
     &                              +x2_z *DBij_
                  dRho_dR(9,iGrid,Ind_x2)=dRho_dR(9,iGrid,Ind_x2)
     &                              +x2_t *DAij_
                  dRho_dR(10,iGrid,Ind_x2)=dRho_dR(10,iGrid,Ind_x2)
     &                              +x2_t *DBij_
                  dRho_dR(11,iGrid,Ind_x2)=dRho_dR(11,iGrid,Ind_x2)
     &                              +x2_r *DAij_
                  dRho_dR(12,iGrid,Ind_x2)=dRho_dR(12,iGrid,Ind_x2)
     &                              +x2_r *DBij_
               End Do ! iGrid
            Else If (Ind_x1.ne.0) Then
               Do iGrid = 1, mGrid
                  x1   =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  x1_x =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x1_y =TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  x1_z =TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  x1_t =TabAO1(dx2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dxdz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  x1_r=(TabAO1(dx3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dxdy2,iGrid,iCB_Eff)
     &                 +TabAO1(dxdz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +x1_t
                  dRho_dR(1,iGrid,Ind_x1)=dRho_dR(1,iGrid,Ind_x1)
     &                              +x1   *DAij_
                  dRho_dR(2,iGrid,Ind_x1)=dRho_dR(2,iGrid,Ind_x1)
     &                              +x1   *DBij_
                  dRho_dR(3,iGrid,Ind_x1)=dRho_dR(3,iGrid,Ind_x1)
     &                              +x1_x *DAij_
                  dRho_dR(4,iGrid,Ind_x1)=dRho_dR(4,iGrid,Ind_x1)
     &                              +x1_y *DAij_
                  dRho_dR(5,iGrid,Ind_x1)=dRho_dR(5,iGrid,Ind_x1)
     &                              +x1_z *DAij_
                  dRho_dR(6,iGrid,Ind_x1)=dRho_dR(6,iGrid,Ind_x1)
     &                              +x1_x *DBij_
                  dRho_dR(7,iGrid,Ind_x1)=dRho_dR(7,iGrid,Ind_x1)
     &                              +x1_y *DBij_
                  dRho_dR(8,iGrid,Ind_x1)=dRho_dR(8,iGrid,Ind_x1)
     &                              +x1_z *DBij_
                  dRho_dR(9,iGrid,Ind_x1)=dRho_dR(9,iGrid,Ind_x1)
     &                              +x1_t *DAij_
                  dRho_dR(10,iGrid,Ind_x1)=dRho_dR(10,iGrid,Ind_x1)
     &                              +x1_t *DBij_
                  dRho_dR(11,iGrid,Ind_x1)=dRho_dR(11,iGrid,Ind_x1)
     &                              +x1_r *DAij_
                  dRho_dR(12,iGrid,Ind_x1)=dRho_dR(12,iGrid,Ind_x1)
     &                              +x1_r *DBij_
               End Do ! iGrid
            Else If (Ind_x2.ne.0) Then
               Do iGrid = 1, mGrid
                  x2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  x2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
                  x2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dxdy2,iGrid,jCB_Eff)
     &                 +TabAO2(dxdz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
     &                 +x2_t
                  dRho_dR(1,iGrid,Ind_x2)=dRho_dR(1,iGrid,Ind_x2)
     &                              +x2   *DAij_
                  dRho_dR(2,iGrid,Ind_x2)=dRho_dR(2,iGrid,Ind_x2)
     &                              +x2   *DBij_
                  dRho_dR(3,iGrid,Ind_x2)=dRho_dR(3,iGrid,Ind_x2)
     &                              +x2_x *DAij_
                  dRho_dR(4,iGrid,Ind_x2)=dRho_dR(4,iGrid,Ind_x2)
     &                              +x2_y *DAij_
                  dRho_dR(5,iGrid,Ind_x2)=dRho_dR(5,iGrid,Ind_x2)
     &                              +x2_z *DAij_
                  dRho_dR(6,iGrid,Ind_x2)=dRho_dR(6,iGrid,Ind_x2)
     &                              +x2_x *DBij_
                  dRho_dR(7,iGrid,Ind_x2)=dRho_dR(7,iGrid,Ind_x2)
     &                              +x2_y *DBij_
                  dRho_dR(8,iGrid,Ind_x2)=dRho_dR(8,iGrid,Ind_x2)
     &                              +x2_z *DBij_
                  dRho_dR(9,iGrid,Ind_x2)=dRho_dR(9,iGrid,Ind_x2)
     &                              +x2_t *DAij_
                  dRho_dR(10,iGrid,Ind_x2)=dRho_dR(10,iGrid,Ind_x2)
     &                              +x2_t *DBij_
                  dRho_dR(11,iGrid,Ind_x2)=dRho_dR(11,iGrid,Ind_x2)
     &                              +x2_r *DAij_
                  dRho_dR(12,iGrid,Ind_x2)=dRho_dR(12,iGrid,Ind_x2)
     &                              +x2_r *DBij_
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
                  y1   =TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  y1_x =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  y1_y =TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y1_z =TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  y1_t =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  y1_r=(TabAO1(dx2dy,iGrid,iCB_Eff)
     &                 +TabAO1(dy3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dydz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +y1_t
                  y2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydx ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
                  y2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dy,iGrid,jCB_Eff)
     &                 +TabAO2(dy3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dydz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
     &                 +y2_t
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                              +y1   *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                              +y1   *DBij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                              +y1_x *DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                              +y1_y *DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                              +y1_z *DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                              +y1_x *DBij_
                  dRho_dR(7,iGrid,Ind_y1)=dRho_dR(7,iGrid,Ind_y1)
     &                              +y1_y *DBij_
                  dRho_dR(8,iGrid,Ind_y1)=dRho_dR(8,iGrid,Ind_y1)
     &                              +y1_z *DBij_
                  dRho_dR(9,iGrid,Ind_y1)=dRho_dR(9,iGrid,Ind_y1)
     &                              +y1_t *DAij_
                  dRho_dR(10,iGrid,Ind_y1)=dRho_dR(10,iGrid,Ind_y1)
     &                              +y1_t *DBij_
                  dRho_dR(11,iGrid,Ind_y1)=dRho_dR(11,iGrid,Ind_y1)
     &                              +y1_r *DAij_
                  dRho_dR(12,iGrid,Ind_y1)=dRho_dR(12,iGrid,Ind_y1)
     &                              +y1_r *DBij_
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                              +y2   *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                              +y2   *DBij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                              +y2_x *DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                              +y2_y *DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                              +y2_z *DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                              +y2_x *DBij_
                  dRho_dR(7,iGrid,Ind_y2)=dRho_dR(7,iGrid,Ind_y2)
     &                              +y2_y *DBij_
                  dRho_dR(8,iGrid,Ind_y2)=dRho_dR(8,iGrid,Ind_y2)
     &                              +y2_z *DBij_
                  dRho_dR(9,iGrid,Ind_y2)=dRho_dR(9,iGrid,Ind_y2)
     &                              +y2_t *DAij_
                  dRho_dR(10,iGrid,Ind_y2)=dRho_dR(10,iGrid,Ind_y2)
     &                              +y2_t *DBij_
                  dRho_dR(11,iGrid,Ind_y2)=dRho_dR(11,iGrid,Ind_y2)
     &                              +y2_r *DAij_
                  dRho_dR(12,iGrid,Ind_y2)=dRho_dR(12,iGrid,Ind_y2)
     &                              +y2_r *DBij_
               End Do ! iGrid
            Else If (Ind_y1.ne.0) Then
               Do iGrid = 1, mGrid
                  y1   =TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
                  y1_x =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  y1_y =TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y1_z =TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  y1_t =TabAO1(dydx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dydz ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  y1_r=(TabAO1(dx2dy,iGrid,iCB_Eff)
     &                 +TabAO1(dy3  ,iGrid,iCB_Eff)
     &                 +TabAO1(dydz2,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +y1_t
                  dRho_dR(1,iGrid,Ind_y1)=dRho_dR(1,iGrid,Ind_y1)
     &                              +y1   *DAij_
                  dRho_dR(2,iGrid,Ind_y1)=dRho_dR(2,iGrid,Ind_y1)
     &                              +y1   *DBij_
                  dRho_dR(3,iGrid,Ind_y1)=dRho_dR(3,iGrid,Ind_y1)
     &                              +y1_x *DAij_
                  dRho_dR(4,iGrid,Ind_y1)=dRho_dR(4,iGrid,Ind_y1)
     &                              +y1_y *DAij_
                  dRho_dR(5,iGrid,Ind_y1)=dRho_dR(5,iGrid,Ind_y1)
     &                              +y1_z *DAij_
                  dRho_dR(6,iGrid,Ind_y1)=dRho_dR(6,iGrid,Ind_y1)
     &                              +y1_x *DBij_
                  dRho_dR(7,iGrid,Ind_y1)=dRho_dR(7,iGrid,Ind_y1)
     &                              +y1_y *DBij_
                  dRho_dR(8,iGrid,Ind_y1)=dRho_dR(8,iGrid,Ind_y1)
     &                              +y1_z *DBij_
                  dRho_dR(9,iGrid,Ind_y1)=dRho_dR(9,iGrid,Ind_y1)
     &                              +y1_t *DAij_
                  dRho_dR(10,iGrid,Ind_y1)=dRho_dR(10,iGrid,Ind_y1)
     &                              +y1_t *DBij_
                  dRho_dR(11,iGrid,Ind_y1)=dRho_dR(11,iGrid,Ind_y1)
     &                              +y1_r *DAij_
                  dRho_dR(12,iGrid,Ind_y1)=dRho_dR(12,iGrid,Ind_y1)
     &                              +y1_r *DBij_
               End Do ! iGrid
            Else If (Ind_y2.ne.0) Then
               Do iGrid = 1, mGrid
                  y2   =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_x =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydx ,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_y =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_z =TabAO1(F    ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  y2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdy ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy2  ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
                  y2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dy,iGrid,jCB_Eff)
     &                 +TabAO2(dy3  ,iGrid,jCB_Eff)
     &                 +TabAO2(dydz2,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
     &                 +y2_t
                  dRho_dR(1,iGrid,Ind_y2)=dRho_dR(1,iGrid,Ind_y2)
     &                              +y2   *DAij_
                  dRho_dR(2,iGrid,Ind_y2)=dRho_dR(2,iGrid,Ind_y2)
     &                              +y2   *DBij_
                  dRho_dR(3,iGrid,Ind_y2)=dRho_dR(3,iGrid,Ind_y2)
     &                              +y2_x *DAij_
                  dRho_dR(4,iGrid,Ind_y2)=dRho_dR(4,iGrid,Ind_y2)
     &                              +y2_y *DAij_
                  dRho_dR(5,iGrid,Ind_y2)=dRho_dR(5,iGrid,Ind_y2)
     &                              +y2_z *DAij_
                  dRho_dR(6,iGrid,Ind_y2)=dRho_dR(6,iGrid,Ind_y2)
     &                              +y2_x *DBij_
                  dRho_dR(7,iGrid,Ind_y2)=dRho_dR(7,iGrid,Ind_y2)
     &                              +y2_y *DBij_
                  dRho_dR(8,iGrid,Ind_y2)=dRho_dR(8,iGrid,Ind_y2)
     &                              +y2_z *DBij_
                  dRho_dR(9,iGrid,Ind_y2)=dRho_dR(9,iGrid,Ind_y2)
     &                              +y2_t *DAij_
                  dRho_dR(10,iGrid,Ind_y2)=dRho_dR(10,iGrid,Ind_y2)
     &                              +y2_t *DBij_
                  dRho_dR(11,iGrid,Ind_y2)=dRho_dR(11,iGrid,Ind_y2)
     &                              +y2_r *DAij_
                  dRho_dR(12,iGrid,Ind_y2)=dRho_dR(12,iGrid,Ind_y2)
     &                              +y2_r *DBij_
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
                  z1   =TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
                  z1_x =TabAO1(dzdx,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  z1_y =TabAO1(dzdy,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  z1_z =TabAO1(dz2 ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z1_t =TabAO1(dzdx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dzdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  z1_r=(TabAO1(dx2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dy2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dz3  ,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +z1_t
                  z2   =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_x =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdx,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_y =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdy,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2  ,iGrid,jCB_Eff)*Two
                  z2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dy2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dz3  ,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
     &                 +z2_t
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                              +z1   *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                              +z1   *DBij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                              +z1_x *DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                              +z1_y *DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                              +z1_z *DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                              +z1_x *DBij_
                  dRho_dR(7,iGrid,Ind_z1)=dRho_dR(7,iGrid,Ind_z1)
     &                              +z1_y *DBij_
                  dRho_dR(8,iGrid,Ind_z1)=dRho_dR(8,iGrid,Ind_z1)
     &                              +z1_z *DBij_
                  dRho_dR(9,iGrid,Ind_z1)=dRho_dR(9,iGrid,Ind_z1)
     &                              +z1_t *DBij_
                  dRho_dR(10,iGrid,Ind_z1)=dRho_dR(10,iGrid,Ind_z1)
     &                              +z1_t *DBij_
                  dRho_dR(11,iGrid,Ind_z1)=dRho_dR(11,iGrid,Ind_z1)
     &                              +z1_r *DAij_
                  dRho_dR(12,iGrid,Ind_z1)=dRho_dR(12,iGrid,Ind_z1)
     &                              +z1_r *DBij_
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                              +z2   *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                              +z2   *DBij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                              +z2_x *DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                              +z2_y *DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                              +z2_z *DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                              +z2_x *DBij_
                  dRho_dR(7,iGrid,Ind_z2)=dRho_dR(7,iGrid,Ind_z2)
     &                              +z2_y *DBij_
                  dRho_dR(8,iGrid,Ind_z2)=dRho_dR(8,iGrid,Ind_z2)
     &                              +z2_z *DBij_
                  dRho_dR(9,iGrid,Ind_z2)=dRho_dR(9,iGrid,Ind_z2)
     &                              +z2_t *DAij_
                  dRho_dR(10,iGrid,Ind_z2)=dRho_dR(10,iGrid,Ind_z2)
     &                              +z2_t *DBij_
                  dRho_dR(11,iGrid,Ind_z2)=dRho_dR(11,iGrid,Ind_z2)
     &                              +z2_r *DAij_
                  dRho_dR(12,iGrid,Ind_z2)=dRho_dR(12,iGrid,Ind_z2)
     &                              +z2_r *DBij_
               End Do ! iGrid
            Else If (Ind_z1.ne.0) Then
               Do iGrid = 1, mGrid
                  z1   =TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
                  z1_x =TabAO1(dzdx,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)
                  z1_y =TabAO1(dzdy,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)
                  z1_z =TabAO1(dz2 ,iGrid,iCB_Eff)*
     &                  TabAO2(F      ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z1_t =TabAO1(dzdx ,iGrid,iCB_Eff)*
     &                  TabAO2(dx   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dzdy ,iGrid,iCB_Eff)*
     &                  TabAO2(dy   ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)*Two
                  z1_r=(TabAO1(dx2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dy2dz,iGrid,iCB_Eff)
     &                 +TabAO1(dz3  ,iGrid,iCB_Eff))*
     &                  TabAO2(F    ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dy2  ,iGrid,jCB_Eff)
     &                 +TabAO2(dz2  ,iGrid,jCB_Eff))
     &                 +z1_t
                  dRho_dR(1,iGrid,Ind_z1)=dRho_dR(1,iGrid,Ind_z1)
     &                              +z1   *DAij_
                  dRho_dR(2,iGrid,Ind_z1)=dRho_dR(2,iGrid,Ind_z1)
     &                              +z1   *DBij_
                  dRho_dR(3,iGrid,Ind_z1)=dRho_dR(3,iGrid,Ind_z1)
     &                              +z1_x *DAij_
                  dRho_dR(4,iGrid,Ind_z1)=dRho_dR(4,iGrid,Ind_z1)
     &                              +z1_y *DAij_
                  dRho_dR(5,iGrid,Ind_z1)=dRho_dR(5,iGrid,Ind_z1)
     &                              +z1_z *DAij_
                  dRho_dR(6,iGrid,Ind_z1)=dRho_dR(6,iGrid,Ind_z1)
     &                              +z1_x *DBij_
                  dRho_dR(7,iGrid,Ind_z1)=dRho_dR(7,iGrid,Ind_z1)
     &                              +z1_y *DBij_
                  dRho_dR(8,iGrid,Ind_z1)=dRho_dR(8,iGrid,Ind_z1)
     &                              +z1_z *DBij_
                  dRho_dR(9,iGrid,Ind_z1)=dRho_dR(9,iGrid,Ind_z1)
     &                              +z1_t *DBij_
                  dRho_dR(10,iGrid,Ind_z1)=dRho_dR(10,iGrid,Ind_z1)
     &                              +z1_t *DBij_
                  dRho_dR(11,iGrid,Ind_z1)=dRho_dR(11,iGrid,Ind_z1)
     &                              +z1_r *DAij_
                  dRho_dR(12,iGrid,Ind_z1)=dRho_dR(12,iGrid,Ind_z1)
     &                              +z1_r *DBij_
               End Do ! iGrid
            Else If (Ind_z2.ne.0) Then
               Do iGrid = 1, mGrid
                  z2   =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_x =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdx,iGrid,jCB_Eff)
     &                 +TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_y =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dzdy,iGrid,jCB_Eff)
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_z =TabAO1(F      ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2 ,iGrid,jCB_Eff)
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
                  z2_t =TabAO1(dx   ,iGrid,iCB_Eff)*
     &                  TabAO2(dxdz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dy   ,iGrid,iCB_Eff)*
     &                  TabAO2(dydz ,iGrid,jCB_Eff)*Two
     &                 +TabAO1(dz   ,iGrid,iCB_Eff)*
     &                  TabAO2(dz2  ,iGrid,jCB_Eff)*Two
                  z2_r =TabAO1(F    ,iGrid,iCB_Eff)*
     &                 (TabAO2(dx2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dy2dz,iGrid,jCB_Eff)
     &                 +TabAO2(dz3  ,iGrid,jCB_Eff))
     &                +(TabAO1(dx2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dy2  ,iGrid,iCB_Eff)
     &                 +TabAO1(dz2  ,iGrid,iCB_Eff))*
     &                  TabAO2(dz   ,iGrid,jCB_Eff)
     &                 +z2_t
                  dRho_dR(1,iGrid,Ind_z2)=dRho_dR(1,iGrid,Ind_z2)
     &                              +z2   *DAij_
                  dRho_dR(2,iGrid,Ind_z2)=dRho_dR(2,iGrid,Ind_z2)
     &                              +z2   *DBij_
                  dRho_dR(3,iGrid,Ind_z2)=dRho_dR(3,iGrid,Ind_z2)
     &                              +z2_x *DAij_
                  dRho_dR(4,iGrid,Ind_z2)=dRho_dR(4,iGrid,Ind_z2)
     &                              +z2_y *DAij_
                  dRho_dR(5,iGrid,Ind_z2)=dRho_dR(5,iGrid,Ind_z2)
     &                              +z2_z *DAij_
                  dRho_dR(6,iGrid,Ind_z2)=dRho_dR(6,iGrid,Ind_z2)
     &                              +z2_x *DBij_
                  dRho_dR(7,iGrid,Ind_z2)=dRho_dR(7,iGrid,Ind_z2)
     &                              +z2_y *DBij_
                  dRho_dR(8,iGrid,Ind_z2)=dRho_dR(8,iGrid,Ind_z2)
     &                              +z2_z *DBij_
                  dRho_dR(9,iGrid,Ind_z2)=dRho_dR(9,iGrid,Ind_z2)
     &                              +z2_t *DAij_
                  dRho_dR(10,iGrid,Ind_z2)=dRho_dR(10,iGrid,Ind_z2)
     &                              +z2_t *DBij_
                  dRho_dR(11,iGrid,Ind_z2)=dRho_dR(11,iGrid,Ind_z2)
     &                              +z2_r *DAij_
                  dRho_dR(12,iGrid,Ind_z2)=dRho_dR(12,iGrid,Ind_z2)
     &                              +z2_r *DBij_
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
