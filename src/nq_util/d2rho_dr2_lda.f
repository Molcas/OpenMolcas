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
      Subroutine d2Rho_dR2_LDA(Dens,nDens,nD,dRho_dR,d2Rho_dR2,ndRho_dR,
     &                         mGrid,list_s,nlist_s,
     &                         TabAO,ipTabAO,mAO,nTabAO,
     &                         nGrad_Eff,list_g,
     &                         Grid_type,Fixed_Grid,Fact,ndc,TabAOMax,
     &                         T_X,list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      use Center_Info
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
     &        ipTabAO(nlist_s), list_bas(2,nlist_s), Index(nIndex)
      Integer Grid_Type, Fixed_Grid
      Real*8 Dens(nDens,nD), TabAO(nTabAO), Fact(ndc**2),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff), Phase(3,2),
     &       TabAOMax(nlist_s),
     &       d2Rho_dR2(ndRho_dR,mGrid,nGrad_Eff,nGrad_Eff)
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
      If (Debug) Then
         Call RecPrt('d2Rho_dR2_LDA:Dens',' ',Dens,nDens,nD)
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
            Call RecPrt('d2Rho_dR2_LDA: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      Call FZero(dRho_dR,ndRho_dR*mGrid*nGrad_Eff)
      Call FZero(d2Rho_dR2,ndRho_dR*mGrid*nGrad_Eff**2)
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
         index_i    =list_bas(2,ilist_s)
*
         Phase(1,1)=DBLE(dc(mdci)%nStab)
         Phase(2,1)=DBLE(dc(mdci)%nStab)
         Phase(3,1)=DBLE(dc(mdci)%nStab)
*
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
            mdcj       =iSD(10,jS)
            jShell     =iSD(11,jS)
            index_j    =list_bas(2,jlist_s)
*
            Phase(1,2)=DBLE(dc(mdcj)%nStab)
            Phase(2,2)=DBLE(dc(mdcj)%nStab)
            Phase(3,2)=DBLE(dc(mdcj)%nStab)
*
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
            ix=iDAMax_(mDij,DeDe(ip_D_a),1)
            DMax_ij=Abs(DeDe(ip_D_a-1+ix))
            If (nD.ne.1) Then
               ix=iDAMax_(mDij,DeDe(ip_D_b),1)
               DMax_ij=DMax_ij+Abs(DeDe(ip_D_a-1+ix))
            End If
            If (TMax_i*TMax_j*DMax_ij.lt.T_X) Go To 98
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*) 'd2Rho_dR2_LDA'
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
               tmp = Phase(1,1)
               Phase(1,1)=Phase(1,2)
               Phase(1,2)=tmp
               tmp = Phase(2,1)
               Phase(2,1)=Phase(2,2)
               Phase(2,2)=tmp
               tmp = Phase(3,1)
               Phase(3,1)=Phase(3,2)
               Phase(3,2)=tmp
            End If
            If (nD.eq.1) Then
               If (iShell.ge.jShell) Then
               Call Do_Rho2ha(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,Phase,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho2ha(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,Phase,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho2h_(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,Phase,
     &                       T_X,TMax_i*TMax_j,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho2h_(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Deg,IndGrd_Eff,Phase,
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
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('d2Rho_dR2_LDA: dRho_dR',' ',dRho_dR,
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
      Subroutine Do_Rho2ha(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                     DAij,          mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,Phase,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(   mGrid,nGrad_Eff), DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp),
     &       Phase(3,2), d2Rho_dR2(  mGrid,nGrad_Eff,nGrad_Eff)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind1(3), Ind2(3,3)
      Data Ind1/2,3,4/, Ind2/5, 6, 8,
     &                       6, 7, 9,
     &                       8, 9,10/
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
               I_Ai=IndGrd_Eff(iCar,1)
               I_Bi=IndGrd_Eff(iCar,2)
               i = Ind1(iCar)
*
*              Do the gradient term
*
               If (I_Ai.ne.0.and.I_Bi.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_1=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_ij_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                     dRho_dR(iGrid,I_Ai)=dRho_dR(iGrid,I_Ai)
     &                                   +Prod_ij_1*DAij_*Phase(iCar,1)
                     dRho_dR(iGrid,I_Bi)=dRho_dR(iGrid,I_Bi)
     &                                   +Prod_ij_2*DAij_*Phase(iCar,2)
                  End Do ! iGrid
               Else If (I_Ai.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_1=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(iGrid,I_Ai)=dRho_dR(iGrid,I_Ai)
     &                                   +Prod_ij_1*DAij_*Phase(iCar,1)
                  End Do ! iGrid
               Else If (I_Bi.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                     dRho_dR(iGrid,I_Bi)=dRho_dR(iGrid,I_Bi)
     &                                   +Prod_ij_2*DAij_*Phase(iCar,2)
                  End Do ! iGrid
               End If
*
*              Do the hessian term
*
               Do jCar = 1, 3
                  I_Aj=IndGrd_Eff(jCar,1)
                  I_Bj=IndGrd_Eff(jCar,2)
                  j = Ind1(jCar)
                  ij = Ind2(iCar,jCar)
*
                  If (I_Ai*I_Aj.ne.0) Then
                     Do iGrid = 1, mGrid
                        AAij=
     &                  TabAO1(ij,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                        d2Rho_dR2(iGrid,I_Ai,I_Aj) =
     &                    d2Rho_dR2(iGrid,I_Ai,I_Aj)
     &                    + AAij*DAij_*Phase(iCar,1)*Phase(jCar,1)
                     End Do
                  End if
                  If (I_Ai*I_Bj.ne.0) Then
                     Do iGrid = 1, mGrid
                        ABij=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                        d2Rho_dR2(iGrid,I_Ai,I_Bj) =
     &                    d2Rho_dR2(iGrid,I_Ai,I_Bj)
     &                    + ABij*DAij_*Phase(iCar,1)*Phase(jCar,2)
                     End Do
                  End if
                  If (I_Aj*I_Bi.ne.0) Then
                     Do iGrid = 1, mGrid
                        ABji=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                        d2Rho_dR2(iGrid,I_Aj,I_Bi) =
     &                    d2Rho_dR2(iGrid,I_Aj,I_Bi)
     &                    + ABji*DAij_*Phase(jCar,1)*Phase(iCar,2)
                     End Do
                  End if
                  If (I_Bi*I_Bj.ne.0) Then
                     Do iGrid = 1, mGrid
                        BBij=
     &                  TabAO1(1,iGrid,iCB_Eff)*TabAO2(ij,iGrid,jCB_Eff)
                        d2Rho_dR2(iGrid,I_Bj,I_Bi) =
     &                    d2Rho_dR2(iGrid,I_Bj,I_Bi)
     &                    + BBij*DAij_*Phase(iCar,2)*Phase(jCar,2)
                     End Do
                  End if
*
               End Do  ! jCar
*
            End Do     ! iCar
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho2h_(dRho_dR,d2Rho_dR2,mGrid,nGrad_Eff,
     &                     DAij,DBij,     mAO,
     &                     TabAO1,iBas,iBas_Eff,iCmp,
     &                     TabAO2,jBas,jBas_Eff,jCmp,
     &                     Fact,IndGrd_Eff,Phase,T_X,TMax_ij,
     &                     Index_i,Index_j)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dRho_dR(2,mGrid,nGrad_Eff),
     &       DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp),
     &       Phase(3,2), d2Rho_dR2(2,mGrid,nGrad_Eff,nGrad_Eff)
      Integer IndGrd_Eff(3,2), Index_i(iBas_Eff*iCmp),
     &                         Index_j(jBas_Eff*jCmp)
      Integer Ind1(3), Ind2(3,3)
      Data Ind1/2,3,4/, Ind2/5, 6, 8,
     &                       6, 7, 9,
     &                       8, 9,10/
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
            If (TMax_ij*Abs(DAij_).lt.T_X .and.
     &          TMax_ij*Abs(DBij_).lt.T_X) Go To 99
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3
               I_Ai=IndGrd_Eff(iCar,1)
               I_Bi=IndGrd_Eff(iCar,2)
               i = Ind1(iCar)
*
*              Do the gradient term
*
               If (I_Ai.ne.0.and.I_Bi.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_1=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     Prod_ij_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,I_Ai)=dRho_dR(1,iGrid,I_Ai)
     &                                   +Prod_ij_1*DAij_*Phase(iCar,1)
                     dRho_dR(2,iGrid,I_Ai)=dRho_dR(2,iGrid,I_Ai)
     &                                   +Prod_ij_1*DBij_*Phase(iCar,1)
                     dRho_dR(1,iGrid,I_Bi)=dRho_dR(1,iGrid,I_Bi)
     &                                   +Prod_ij_2*DAij_*Phase(iCar,2)
                     dRho_dR(2,iGrid,I_Bi)=dRho_dR(2,iGrid,I_Bi)
     &                                   +Prod_ij_2*DBij_*Phase(iCar,2)
                  End Do ! iGrid
               Else If (I_Ai.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_1=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,I_Ai)=dRho_dR(1,iGrid,I_Ai)
     &                                   +Prod_ij_1*DAij_*Phase(iCar,1)
                     dRho_dR(2,iGrid,I_Ai)=dRho_dR(2,iGrid,I_Ai)
     &                                   +Prod_ij_1*DBij_*Phase(iCar,1)
                  End Do ! iGrid
               Else If (I_Bi.ne.0) Then
                  Do iGrid = 1, mGrid
                     Prod_ij_2=
     &                   TabAO1(1,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                     dRho_dR(1,iGrid,I_Bi)=dRho_dR(1,iGrid,I_Bi)
     &                                   +Prod_ij_2*DAij_*Phase(iCar,2)
                     dRho_dR(2,iGrid,I_Bi)=dRho_dR(2,iGrid,I_Bi)
     &                                   +Prod_ij_2*DBij_*Phase(iCar,2)
                  End Do ! iGrid
               End If
*
*              Do the hessian term
*
               Do jCar = 1, 3
                  I_Aj=IndGrd_Eff(jCar,1)
                  I_Bj=IndGrd_Eff(jCar,2)
                  j = Ind1(jCar)
                  ij = Ind2(iCar,jCar)
*
                  If (I_Ai*I_Aj.ne.0) Then
                     Do iGrid = 1, mGrid
                        AAij=
     &                  TabAO1(ij,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
                        d2Rho_dR2(1,iGrid,I_Ai,I_Aj) =
     &                    d2Rho_dR2(1,iGrid,I_Ai,I_Aj)
     &                    + AAij*DAij_*Phase(iCar,1)*Phase(jCar,1)
                        d2Rho_dR2(2,iGrid,I_Ai,I_Aj) =
     &                    d2Rho_dR2(2,iGrid,I_Ai,I_Aj)
     &                    + AAij*DBij_*Phase(iCar,1)*Phase(jCar,1)
                     End Do
                  End if
                  If (I_Ai*I_Bj.ne.0) Then
                     Do iGrid = 1, mGrid
                        ABij=
     &                   TabAO1(i,iGrid,iCB_Eff)*TabAO2(j,iGrid,jCB_Eff)
                        d2Rho_dR2(1,iGrid,I_Ai,I_Bj) =
     &                    d2Rho_dR2(1,iGrid,I_Ai,I_Bj)
     &                    + ABij*DAij_*Phase(iCar,1)*Phase(jCar,2)
                        d2Rho_dR2(2,iGrid,I_Ai,I_Bj) =
     &                    d2Rho_dR2(2,iGrid,I_Ai,I_Bj)
     &                    + ABij*DBij_*Phase(iCar,1)*Phase(jCar,2)
                     End Do
                  End if
                  If (I_Aj*I_Bi.ne.0) Then
                     Do iGrid = 1, mGrid
                        ABji=
     &                   TabAO1(j,iGrid,iCB_Eff)*TabAO2(i,iGrid,jCB_Eff)
                        d2Rho_dR2(1,iGrid,I_Aj,I_Bi) =
     &                    d2Rho_dR2(1,iGrid,I_Aj,I_Bi)
     &                    + ABji*DAij_*Phase(jCar,1)*Phase(iCar,2)
                        d2Rho_dR2(2,iGrid,I_Aj,I_Bi) =
     &                    d2Rho_dR2(2,iGrid,I_Aj,I_Bi)
     &                    + ABji*DBij_*Phase(jCar,1)*Phase(iCar,2)
                     End Do
                  End if
                  If (I_Bi*I_Bj.ne.0) Then
                     Do iGrid = 1, mGrid
                        BBij=
     &                  TabAO1(1,iGrid,iCB_Eff)*TabAO2(ij,iGrid,jCB_Eff)
                        d2Rho_dR2(1,iGrid,I_Bj,I_Bi) =
     &                    d2Rho_dR2(1,iGrid,I_Bj,I_Bi)
     &                    + BBij*DAij_*Phase(iCar,2)*Phase(jCar,2)
                        d2Rho_dR2(2,iGrid,I_Bj,I_Bi) =
     &                    d2Rho_dR2(2,iGrid,I_Bj,I_Bi)
     &                    + BBij*DBij_*Phase(iCar,2)*Phase(jCar,2)
                     End Do
                  End if
*
               End Do  ! jCar
*
            End Do     ! iCar
*
 99         Continue
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
