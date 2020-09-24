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
      Subroutine Rho_LDA(Dens,nDens,nD,Rho,nRho,mGrid,
     &                   list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                   Fact,mdc,TabAOMax,list_bas,Index,nIndex)
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
      use k2_arrays, only: DeDe, ipDijS
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer list_s(2,nlist_s), ipTabAO(nlist_s), list_bas(2,nlist_s),
     &        Index(nIndex)
      Real*8 Dens(nDens,nD), TabAO(nTabAO), Rho(nRho,mGrid),
     &       Fact(mdc**2), TabAOMax(nlist_s)
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
      Call QEnter('Rho_LDA')
#endif
*define _DEBUGPRINT_DEBUG_
#ifdef _DEBUGPRINT_DEBUG_
      Call RecPrt('Rho_LDA:Dens',' ',Dens,nDens,nD)
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
         mdci  = iSD(10,iSkal)
         Call RecPrt('Rho_LDA: TabAO',' ',
     &               TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
      End Do
#endif
*
      Call FZero(Rho,nRho*mGrid)
*                                                                      *
************************************************************************
*                                                                      *
*
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
         iShell= iSD(11,iSkal)
         mdci  = iSD(10,iSkal)
         index_i=list_bas(2,ilist_s)
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
         lDCRER=NrOpr(iER)
*
         ip_D_a=ipDij+lDCRER*mDij
         ip_D_b=ip_D_a
         If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
         If (nD.ne.1) Then
            ix=iDAMax_(mDij,DeDe(ip_D_a),1)
            iy=iDAMax_(mDij,DeDe(ip_D_b),1)
            DMax_ii=Half*( Abs(DeDe(ip_D_a-1+ix))
     &                    +Abs(DeDe(ip_D_b-1+iy)) )
         Else
            ix=iDAMax_(mDij,DeDe(ip_D_a),1)
            DMax_ii=Abs(DeDe(ip_D_a-1+ix))
         End If
         If (TMax_i*TMax_i*DMax_ii.ge.T_X) Then
            If (nD.eq.1) Then
               Call Do_Rho2a_d(Rho,     mGrid,
     &                         DeDe(ip_D_a),     mAO,
     &                         TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,
     &                         iCmp,Fact(ij),T_X,TMax_i*TMax_i,
     &                         Index(index_i))
            Else
               Call Do_Rho2_d(Rho,     mGrid,
     &                        DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                        TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,
     &                        iCmp,Fact(ij),T_X,TMax_i*TMax_i,
     &                        Index(index_i))
            End If
         End If
*
         Do jlist_s=1,ilist_s-1
            TMax_j=TabAOMax(jlist_s)
            If (TMax_i*TMax_j.lt.T_X) Go To 998
            jSkal = list_s(1,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            kDCRR=list_s(2,jlist_s)
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
#ifdef _DEBUGPRINT_DEBUG_
            Write (6,*)
            Write (6,*) 'iS,jS=',iSkal,jSkal
            Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
            Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
#endif
*
            ij = (mdcj-1)*mdc + mdci
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)
*
            Deg=Two
            If (ilist_s.eq.jlist_s) Deg=One
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
            If (TMax_i*TMax_j*DMax_ij.lt.T_X) Go To 998
#ifdef _DEBUGPRINT_DEBUG_
            Write (6,*) 'Rho_LDA'
            nBB = iBas*jBas
            nCC = iCmp*jCmp
            Write (6,*) 'iShell,jshell=', iShell,jshell
            Write (6,*) 'kDCRE,kDCRR=', kDCRE,kDCRR
            Write (6,*) 'iER,lDCRER=',iER,lDCRER
            Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
            If (nD.ne.1) Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
#endif
*                                                                      *
************************************************************************
*                                                                      *
            If (nD.eq.1) Then
*
               If (iShell.ge.jShell) Then
*
               Call Do_Rho2a(Rho,     mGrid,
     &                       DeDe(ip_D_a),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
*
               Else
*
               Call Do_Rho2a(Rho,     mGrid,
     &                       DeDe(ip_D_a),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
*
               End If
*
            Else
*
               If (iShell.ge.jShell) Then
*
               Call Do_Rho2_(Rho,     mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
*
               Else
*
               Call Do_Rho2_(Rho,     mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Two,T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
*
               End If
            End If
*                                                                      *
************************************************************************
*                                                                      *
 998        Continue
         End Do                      ! jlist_s
 999     Continue
      End Do                         ! ilist_s
*
#ifdef _DEBUGPRINT_DEBUG_
      Call RecPrt('Rho_LDA: Rho','(10F15.6)',Rho,nRho,mGrid)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Dens)
#endif
#ifdef _TIME_
      Call QExit('Rho_LDA')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nSym)
      End




      Subroutine Do_Rho2a(Rho,     mGrid,DAij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,T_X,TMax_ij,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "debug.fh"
      Real*8 Rho(     mGrid), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QEnter('Do_Rho2a')
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
               Prod_ij=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Rho(  iGrid)=Rho(  iGrid)+Prod_ij*DAij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
      Call QExit('Do_Rho2a')
#endif
      Return
      End



      Subroutine Do_Rho2_(Rho,     mGrid,DAij,DBij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,T_X,TMax_ij,Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(   2,mGrid),
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
            Dij_ = Half*(Abs(DAij_)+Abs(DBij_))
            If (TMax_ij*Abs(Dij_).lt.T_X ) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_ij=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Rho(1,iGrid)=Rho(1,iGrid)+Prod_ij*DAij_
               Rho(2,iGrid)=Rho(2,iGrid)+Prod_ij*DBij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End




      Subroutine Do_Rho2a_d(Rho,     mGrid,DAii,
     &                      mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                      Fact,T_X,TMax_ii,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "debug.fh"
      Real*8 Rho(     mGrid), DAii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QEnter('Do_Rho2a_d')
#endif
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB = Index_i(jCB_Eff)
*
         DAjj_=DAii(jCB,jCB)*Fact
         If (TMax_ii*Abs(DAjj_).ge.T_X) Then
            Do iGrid = 1, mGrid
               Prod_jj=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Rho(  iGrid)=Rho(  iGrid)+Prod_jj*DAjj_
            End Do    ! iGrid
         End If
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
            If (TMax_ii*Abs(DAij_).lt.T_X) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_ij=TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Rho(  iGrid)=Rho(  iGrid)+Prod_ij*DAij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
#ifdef _TIME_
      Call QExit('Do_Rho2a_d')
#endif
      Return
      End




      Subroutine Do_Rho2_d(Rho,     mGrid,DAii,DBii,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                    Fact,T_X,TMax_ii,Index_i)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Rho(   2,mGrid),
     &       DAii(iBas*iCmp,iBas*iCmp), DBii(iBas*iCmp,iBas*iCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp)
      Integer Index_i(iBas_Eff*iCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, iBas_Eff*iCmp
         jCB = Index_i(jCB_Eff)
*
         DAjj_=DAii(jCB,jCB)*Fact
         DBjj_=DBii(jCB,jCB)*Fact
         Djj_ =Half*(Abs(DAjj_)+Abs(DBjj_))
         If (TMax_ii*Abs(Djj_).ge.T_X) Then
            Do iGrid = 1, mGrid
               Prod_jj=TabAO1(1,iGrid,jCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Rho(1,iGrid)=Rho(1,iGrid)+Prod_jj*DAjj_
               Rho(2,iGrid)=Rho(2,iGrid)+Prod_jj*DBjj_
            End Do    ! iGrid
         End If
*
         Do iCB_Eff = 1, jCB_Eff-1
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAii(iCB,jCB)*Fact*Two
            DBij_=DBii(iCB,jCB)*Fact*Two
            Dij_ = Half*(Abs(DAij_)+Abs(DBij_))
            If (TMax_ii*Abs(Dij_).lt.T_X) Go To 99
*
            Do iGrid = 1, mGrid
               Prod_ij=TabAO1(1,iGrid,iCB_Eff)*TabAO1(1,iGrid,jCB_Eff)
               Rho(1,iGrid)=Rho(1,iGrid)+Prod_ij*DAij_
               Rho(2,iGrid)=Rho(2,iGrid)+Prod_ij*DBij_
            End Do    ! iGrid
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
