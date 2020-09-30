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
* Copyright (C) 2000,2002, Roland Lindh                                *
************************************************************************
      Subroutine dRho_dR_LDA(Dens,nDens,nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                       nGrad_Eff,list_g,Maps2p,nShell,
     &                       Grid_type,Fixed_Grid,Fact,ndc,TabAOMax,T_X,
     &                       list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
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
      Real*8 Dens(nDens,nD), TabAO(nTabAO), Fact(ndc**2),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff),
     &       TabAOMax(nlist_s)
      Integer IndGrd_Eff(3,2)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      If (Debug) Then
         Call RecPrt('dRho_dR_LDA:Dens',' ',Dens,nDens,nD)
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nlist_s=',nlist_s
         Write (6,*) 'nD=',nD
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iS = list_s(1,ilist_s)
            iCmp  = iSD( 2,iS)
            iBas_Eff  = list_bas(1,ilist_s)
            mTabAO=iBas_Eff*iCmp
            mdci  = iSD(10,iS)
            Call RecPrt('dRho_dR_LDA: TabAO',' ',
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
         iS      =list_s(1,ilist_s)
         iCmp       =iSD( 2,iS)
         iBas_Eff  = list_bas(1,ilist_s)
         ix = iDAMax_(mAO*mGrid*iBas_Eff*iCmp,TabAO(ipTabAO(iList_s)),1)
         TabAOMax(ilist_s)=Abs(TabAO(ipTabAO(ilist_s)-1+ix))
         TMax_i=TabAOMax(ilist_s)
         If (TMax_i.le.T_X) Go To 999
         iBas       =iSD( 3,iS)
         kDCRE   =list_s(2,ilist_s)
         iShell     =iSD(11,iS)
         mdci       =iSD(10,iS)
         iPrim      =iSD( 5,iS)
         index_i    =list_bas(2,ilist_s)
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
            jS      =list_s(1,jlist_s)
            kDCRR   =list_s(2,jlist_s)
            jCmp       =iSD( 2,jS)
            jBas       =iSD( 3,jS)
            jBas_Eff   =list_bas(1,jlist_s)
            jPrim      =iSD( 5,jS)
            mdcj       =iSD(10,jS)
            jShell     =iSD(11,jS)
            index_j    =list_bas(2,jlist_s)
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
#ifdef _DEBUG_
            If (Debug) Then
               Write (6,*) 'dRho_dR_LDA'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jShell=',iShell,jShell
               Write (6,*) 'kDCRE,kDCRR=',kDCRE,kDCRR
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
               Call Do_Rho2da(dRho_dR,   mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho2da(dRho_dR,   mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho2d_(dRho_dR,   mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho2d_(dRho_dR,   mGrid,nGrad_Eff,
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
*
         End Do                      ! jlist_s
 999     Continue
      End Do                         ! ilist_s
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      If (Debug) Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       ndRho_dR*mGrid,nGrad_Eff)
*
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Dens)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(Grid_type)
         Call Unused_integer(Fixed_Grid)
      End If
      End
      Subroutine Do_Rho2da(dRho_dR,   mGrid,nGrad_Eff,
     &                     DAij,          mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   mGrid,nGrad_Eff), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind(3)
      Data Ind/2,3,4/
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB=Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            If (TMax_ij*Abs(DAij_).lt.T_X) Go To 99
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3
               Ind_1=IndGrd_Eff(iCar,1)
               Ind_2=IndGrd_Eff(iCar,2)
               j = Ind(iCar)
               If (Ind_1.ne.0.and.Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_1)=dRho_dR(iGrid,Ind_1)
     &                                   +Prod_1*DAij_
                     dRho_dR(iGrid,Ind_2)=dRho_dR(iGrid,Ind_2)
     &                                   +Prod_2*DAij_
                  End Do ! iGrid
               Else If (Ind_1.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_1)=dRho_dR(iGrid,Ind_1)
     &                                   +Prod_1*DAij_
                  End Do ! iGrid
               Else If (Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(iGrid,Ind_2)=dRho_dR(iGrid,Ind_2)
     &                                   +Prod_2*DAij_
                  End Do ! iGrid
               End If
*
            End Do  ! iCar
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho2d_(dRho_dR,   mGrid,nGrad_Eff,
     &                     DAij,DBij,     mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   2,mGrid,nGrad_Eff),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind(3)
      Data Ind/2,3,4/
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB=Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            DBij_=DBij(iCB,jCB)*Fact
            Dij_ =Half*(Abs(DAij_)+Abs(DBij_))
            If (TMax_ij*Abs(Dij_).lt.T_X) Go To 99
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3
               Ind_1=IndGrd_Eff(iCar,1)
               Ind_2=IndGrd_Eff(iCar,2)
               j = Ind(iCar)
               If (Ind_1.ne.0.and.Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_1)=dRho_dR(1,iGrid,Ind_1)
     &                                     +Prod_1*DAij_
                     dRho_dR(2,iGrid,Ind_1)=dRho_dR(2,iGrid,Ind_1)
     &                                     +Prod_1*DBij_
                     dRho_dR(1,iGrid,Ind_2)=dRho_dR(1,iGrid,Ind_2)
     &                                     +Prod_2*DAij_
                     dRho_dR(2,iGrid,Ind_2)=dRho_dR(2,iGrid,Ind_2)
     &                                     +Prod_2*DBij_
                  End Do ! iGrid
               Else If (Ind_1.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_1=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_1)=dRho_dR(1,iGrid,Ind_1)
     &                                     +Prod_1*DAij_
                     dRho_dR(2,iGrid,Ind_1)=dRho_dR(2,iGrid,Ind_1)
     &                                     +Prod_1*DBij_
                  End Do ! iGrid
               Else If (Ind_2.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,Ind_2)=dRho_dR(1,iGrid,Ind_2)
     &                                     +Prod_2*DAij_
                     dRho_dR(2,iGrid,Ind_2)=dRho_dR(2,iGrid,Ind_2)
     &                                     +Prod_2*DBij_
                  End Do ! iGrid
               End If
*
            End Do ! iCar
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
