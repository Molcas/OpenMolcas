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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Rho_meta_GGA1(nD,Rho,nRho,mGrid,
     &                         list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                         Fact,mdc,TabAOMax,list_bas,Index,nIndex)
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
*             of Lund, SWEDEN.  2000                                   *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "k2.fh"
      Integer list_s(2,nlist_s), ipTabAO(nlist_s), list_bas(2,nlist_s),
     &        Index(nIndex)
      Real*8 Rho(nRho,mGrid), Fact(mdc**2),
     &       TabAO(nTabAO), TabAOMax(nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*define _TIME_
#ifdef _TIME_
      Call QEnter('Rho_meta_GGA')
#endif
*define _DEBUG_
#ifdef _DEBUG_
      Debug=.True.
      If (Debug) Then
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nIrrep=',nIrrep
         Write (6,*) 'nlist_s=',nlist_s
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iSkal = list_s(1,ilist_s)
            iCmp  = iSD( 2,iSkal)
            iBas  = iSD( 3,iSkal)
            mTabAO=iBas*iCmp
            Call RecPrt('Rho_meta_GGA: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      Call FZero(Rho,nRho*mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         ix = iDAMax_(mAO*mGrid*iBas_Eff*iCmp,TabAO(ipTabAO(iList_s)),1)
         TabAOMax(ilist_s)=Abs(TabAO(ipTabAO(ilist_s)-1+ix))
         TMax_i=TabAOMax(ilist_s)
         If (TMax_i.le.T_X) Go To 999
         kDCRE=list_s(2,ilist_s)
         index_i=list_bas(2,ilist_s)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
         nFunc_i=iBas*iCmp
*
         mDij=nFunc_i*nFunc_i
*
*------- Get the Density
*
         ijS=iTri(iShell,iShell)
         ip_Tmp=ipDijs
         Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
*
         ij = (mdci-1)*mdc + mdci
*
         iER=iEOr(kDCRE,kDCRE)
         lDCRER=NrOpr(iER,iOper,nIrrep)
*
         ip_D_a=ipDij+lDCRER*mDij
         ip_D_b=ip_D_a
         If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
         If (nD.ne.1) Then
            ix=iDAMax_(mDij,Work(ip_D_a),1)
            iy=iDAMax_(mDij,Work(ip_D_b),1)
            DMax_ii=Half*( Abs(Work(ip_D_a-1+ix))
     &                    +Abs(Work(ip_D_b-1+iy)) )
         Else
            ix=iDAMax_(mDij,Work(ip_D_a),1)
            DMax_ii=Abs(Work(ip_D_a-1+ix))
         End If
         If (TMax_i*TMax_i*DMax_ii.ge.T_X) Then
            If (nD.eq.1) Then
               Call Do_Rho5a_d(Rho,nRho,mGrid,
     &                         Work(ip_D_a),mAO,
     &                         TabAO(ipTabAO(iList_s)),
     &                         iBas,iBas_Eff,iCmp,
     &                         Fact(ij),T_X,TMax_i*TMax_i,
     &                         Index(index_i))
            Else
               Call Do_Rho5_d(Rho,nRho,mGrid,
     &                        Work(ip_D_a),Work(ip_D_b),mAO,
     &                        TabAO(ipTabAO(iList_s)),
     &                        iBas,iBas_Eff,iCmp,
     &                        Fact(ij),T_X,TMax_i*TMax_i,
     &                        Index(index_i))
            End If
         End If

*
         Do jlist_s=1,ilist_s-1
            TMax_j=TabAOMax(jlist_s)
            If (TMax_i*TMax_j.lt.T_X) Go To 998
            jSkal = list_s(1,jlist_s)
            kDCRR=list_s(2,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
            nFunc_j=jBas*jCmp
*
            mDij=nFunc_i*nFunc_j
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ip_Tmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
#ifdef _DEBUG_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iSkal,jSkal
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*mdc + mdci
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER,iOper,nIrrep)
*
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
            If (nD.ne.1) Then
               ix=iDAMax_(mDij,Work(ip_D_a),1)
               iy=iDAMax_(mDij,Work(ip_D_b),1)
               DMax_ij=Half*( Abs(Work(ip_D_a-1+ix))
     &                       +Abs(Work(ip_D_b-1+iy)) )
            Else
               ix=iDAMax_(mDij,Work(ip_D_b),1)
               DMax_ij=Abs(Work(ip_D_a-1+ix))
            End If
            If (TMax_i*TMax_j*DMax_ij.lt.T_X) Go To 998
#ifdef _DEBUG_
            If (Debug) Then
               Write (6,*) 'Rho_meta_GGA1'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jshell=', iShell,jshell
               Write (6,*) 'kDCRE,kDCRR=', kDCRE,kDCRR
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Write (6,*) 'MaxDe,MaxDCR=',MaxDe,MaxDCR
               Write (6,*) 'ipDeDe,nDeDe_DFT=',ipDeDe, nDeDe_DFT
               Write (6,*) 'ip_D_a=',ip_D_a
               Call RecPrt('DAij',' ',Work(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',Work(ip_D_b),nBB,nCC)
            End If
#endif
*
            If (nD.eq.1) Then
               If (iShell.ge.jShell) Then
               Call Do_Rho5a(Rho,nRho,mGrid,
     &                       Work(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho5a(Rho,nRho,mGrid,
     &                       Work(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho5_(Rho,nRho,mGrid,
     &                       Work(ip_D_a),Work(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho5_(Rho,nRho,mGrid,
     &                       Work(ip_D_a),Work(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            End If
*
 998        Continue
         End Do                      ! jlist_s
 999     Continue
      End Do                         ! ilist_s
*
#ifdef _DEBUG_
      If (Debug) Then
c        Do iGrid=1,mGrid_Eff
c           Write (*,*) (Rho(iRho,iGrid),iRho=1,nRho)
c        End Do
         Call RecPrt('Rho_meta_GGA: Rho',' ',Rho,nRho,mGrid)
      End If
*
#endif
#ifdef _TIME_
      Call QExit('Rho_GGA')
#endif
      Return
      End
      Subroutine Do_Rho5a(Rho,nRho,mGrid,
     &                    DAij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,T_X,TMax_ij,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QEnter('Do_Rho5a')
#endif
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*

         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            If (TMax_ij*Abs(DAij_).lt.T_X) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid)
     &                     + (Prod_22+Prod_33+Prod_44)*DAij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
      Call QExit('Do_Rho5a')
#endif
      Return
      End
      Subroutine Do_Rho5_(Rho,nRho,mGrid,
     &                    DAij,DBij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,T_X,TMax_ij,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
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
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Rho( 1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho( 2,iGrid)=Rho(2,iGrid) +     Prod_11      *DBij_
               Rho( 3,iGrid)=Rho(3,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho( 4,iGrid)=Rho(4,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho( 5,iGrid)=Rho(5,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho( 6,iGrid)=Rho(6,iGrid) + (Prod_21+Prod_12)*DBij_
               Rho( 7,iGrid)=Rho(7,iGrid) + (Prod_31+Prod_13)*DBij_
               Rho( 8,iGrid)=Rho(8,iGrid) + (Prod_41+Prod_14)*DBij_
               Rho( 9,iGrid)=Rho( 9,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DAij_
               Rho(10,iGrid)=Rho(10,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DBij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho5a_d(Rho,nRho,mGrid,
     &                    DAii,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                    Fact,T_X,TMax_ii,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid), DAii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QEnter('Do_Rho5a_d')
#endif
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB=Index_i(jCB_Eff)
*
         DAii_=DAii(jCB,jCB)*Fact
         If (TMax_ii*Abs(DAii_).ge.T_X) Then
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,jCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,jCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,jCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11*DAii_
               Rho(2,iGrid)=Rho(2,iGrid) + Two*Prod_21*DAii_
               Rho(3,iGrid)=Rho(3,iGrid) + Two*Prod_31*DAii_
               Rho(4,iGrid)=Rho(4,iGrid) + Two*Prod_41*DAii_
               Rho(5,iGrid)=Rho(5,iGrid)
     &                     +(Prod_22+Prod_33+Prod_44)*DAii_
            End Do    ! iGrid
         End If
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
            If (TMax_ii*Abs(DAij_).lt.T_X) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid)
     &                     +(Prod_22+Prod_33+Prod_44)*DAij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
      Call QExit('Do_Rho5a_d')
#endif
      Return
      End
      Subroutine Do_Rho5_d(Rho,nRho,mGrid,
     &                     DAii,DBii,
     &                     mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                     Fact,T_X,TMax_ii,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(nRho,mGrid),
     &       DAii(iBas*iCmp,iBas*iCmp), DBii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB=Index_i(jCB_Eff)
*
         DAii_=DAii(jCB,jCB)*Fact
         DBii_=DBii(jCB,jCB)*Fact
         Dii_ =Half*(Abs(DAii_)+Abs(DBii_))
         If (TMax_ii*Abs(Dii_).ge.T_X) Then
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,jCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,jCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,jCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11*DAii_
               Rho(2,iGrid)=Rho(2,iGrid) +     Prod_11*DBii_
               Rho(3,iGrid)=Rho(3,iGrid) + Two*Prod_21*DAii_
               Rho(4,iGrid)=Rho(4,iGrid) + Two*Prod_31*DAii_
               Rho(5,iGrid)=Rho(5,iGrid) + Two*Prod_41*DAii_
               Rho(6,iGrid)=Rho(6,iGrid) + Two*Prod_21*DBii_
               Rho(7,iGrid)=Rho(7,iGrid) + Two*Prod_31*DBii_
               Rho(8,iGrid)=Rho(8,iGrid) + Two*Prod_41*DBii_
               Rho( 9,iGrid)=Rho( 9,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DAii_
               Rho(10,iGrid)=Rho(10,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DBii_
            End Do    ! iGrid
         End If
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB=Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
            DBij_=DBii(iCB,jCB)*Fact*Two
            Dij_ =Half*(Abs(DAij_)+Abs(DBij_))
            If (TMax_ii*Abs(Dij_).lt.T_X ) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_11= TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_21= TabAO1(2,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_12= TabAO1(1,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_31= TabAO1(3,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_13= TabAO1(1,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_41= TabAO1(4,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Prod_14= TabAO1(1,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO1(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO1(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO1(4,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) +     Prod_11      *DBij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(6,iGrid)=Rho(6,iGrid) + (Prod_21+Prod_12)*DBij_
               Rho(7,iGrid)=Rho(7,iGrid) + (Prod_31+Prod_13)*DBij_
               Rho(8,iGrid)=Rho(8,iGrid) + (Prod_41+Prod_14)*DBij_
               Rho( 9,iGrid)=Rho( 9,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DAij_
               Rho(10,iGrid)=Rho(10,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DBij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
