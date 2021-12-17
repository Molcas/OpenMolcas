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
* Copyright (C) 2000,2021, Roland Lindh                                *
************************************************************************
      Subroutine Mk_Rho(nD,mGrid,list_s,nlist_s,Fact,mdc,list_bas,Index,
     &                  nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      use nq_grid, only: Rho, TabAO, Dens_AO, Grid_AO
      use nq_grid, only: GradRho, Sigma, Tau, Lapl
#ifdef _DEBUGPRINT_
      use nq_grid, only: nRho
#endif
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer list_s(2,nlist_s), list_bas(2,nlist_s), Index(nIndex)
      Real*8 Fact(mdc**2)
      Integer ipD(2)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
      nAO = SIZE(Dens_AO,1)
      Dens_AO(:,:,:)=Zero
      mAO = SIZE(TabAO,1)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'mAO=',mAO
      Write (6,*) 'mGrid=',mGrid
      Write (6,*) 'nlist_s=',nlist_s
      Call RecPrt('Rho: TabAO',' ',TabAO,mAO*mGrid,nAO)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
      iOff = 0
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         kDCRE=list_s(2,ilist_s)
         iShell= iSD(11,iSkal)
         mdci  = iSD(10,iSkal)
         index_i=list_bas(2,ilist_s)
         nFunc_i=iBas*iCmp
         n_iBas=iBas_Eff*iCmp
*
         jOff = 0
         Do jlist_s=1,ilist_s
            jSkal = list_s(1,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            kDCRR=list_s(2,jlist_s)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
            nFunc_j=jBas*jCmp
            n_jBas=jBas_Eff*jCmp
*
            mDij=nFunc_i*nFunc_j
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ip_Tmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
#ifdef _DEBUGPRINT_
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
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
            ipD(1)=ip_D_a
            ipD(2)=ip_D_b
*
#ifdef _DEBUGPRINT_
            Write (6,*) 'Rho_LDA'
            nBB = iBas*jBas
            nCC = iCmp*jCmp
            Write (6,*) 'iShell,jshell=', iShell,jshell
            Write (6,*) 'kDCRE,kDCRR=', kDCRE,kDCRR
            Write (6,*) 'iER,lDCRER=',iER,lDCRER
            Write (6,*) 'Fact(ij)=',Fact(ij)
            Call RecPrt('DAij',' ',DeDe(ip_D_a),nFunc_i,nFunc_j)
            If (nD.ne.1) Call RecPrt('DBij',' ',DeDe(ip_D_b),nFunc_i,
     &                                                       nFunc_j)
#endif
*                                                                      *
************************************************************************
*                                                                      *
            Do iD = 1, nD
               Do j_R = 1, n_jBas
                  jCB = Index(index_j-1+j_R)    ! Real index
                  j_A = j_R + jOff            ! Absolute index
*
                  Do i_R = 1, n_iBas
                     iCB = Index(index_i-1+i_R)
                     i_A = i_R + iOff
*
                     ij_D = (jCB-1)*nFunc_i + iCB - 1
                     DAij =DeDe(ipD(iD)+ij_D)*Fact(ij)
                     Dens_AO(i_A,j_A,iD) = DAij
                     Dens_AO(j_A,i_A,iD) = DAij
*
                  End Do          ! iCB
               End Do             ! jCB
            End Do
*                                                                      *
************************************************************************
*                                                                      *
            jOff = jOff + jBas_Eff*jCmp
         End Do                      ! jlist_s
         iOff = iOff + iBas_Eff*iCmp
      End Do                         ! ilist_s
*
      Call DGEMM_('N','N',mAO*mGrid,nAO*nD,nAO,
     &            One,TabAO,mAO*mGrid,
     &                Dens_AO,nAO,
     &            Zero,Grid_AO,mAO*mGrid)

      If (Functional_Type.eq.LDA_Type) Then
         Rho(:,1:mGrid)=Zero
         Do iD = 1, nD
            Do iAO = 1, nAO
               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
               End Do
            End Do
         End Do
      Else If (Functional_Type.eq.GGA_Type) Then
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO
               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
               End Do
            End Do
         End Do
      Else If (Functional_Type.eq.meta_GGA_Type1) Then
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Tau(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO
               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  Tau(iD,iGrid) = Tau(iD,iGrid)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
               End Do
            End Do
         End Do
      Else If (Functional_Type.eq.meta_GGA_Type2) Then
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Tau(:,1:mGrid)=Zero
         Lapl(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO
               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  Tau(iD,iGrid) = Tau(iD,iGrid)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
                  Lapl(iD,iGrid) = Lapl(iD,iGrid)
     &                          + TabAO( 1,iGrid,iAO) *
     &                          ( Grid_AO( 5,iGrid,iAO,iD)
     &                          + Grid_AO( 8,iGrid,iAO,iD)
     &                          + Grid_AO(10,iGrid,iAO,iD) )
     &                          + Two *
     &                          ( Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD) )
     &                          +(TabAO( 5,iGrid,iAO)
     &                          + TabAO( 8,iGrid,iAO)
     &                          + TabAO(10,iGrid,iAO) )
     &                          * Grid_AO(1,iGrid,iAO,iD)
               End Do
            End Do
         End Do
      Else If (Functional_Type.eq.CASDFT_Type) Then
         Rho(:,1:mGrid)=Zero
         GradRho(:,1:mGrid)=Zero
         Tau(:,1:mGrid)=Zero
         Do iD = 1, nD
            ix = (iD-1)*3 + 1
            iy = (iD-1)*3 + 2
            iz = (iD-1)*3 + 3
            Do iAO = 1, nAO
               Do iGrid = 1, mGrid
                  Rho(iD,iGrid) = Rho(iD,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(ix,iGrid)=GradRho(ix,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(2,iGrid,iAO)
     &                          + Grid_AO(2,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iy,iGrid)=GradRho(iy,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(3,iGrid,iAO)
     &                          + Grid_AO(3,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  GradRho(iz,iGrid)=GradRho(iz,iGrid)
     &                          + Grid_AO(1,iGrid,iAO,iD)
     &                          * TabAO(4,iGrid,iAO)
     &                          + Grid_AO(4,iGrid,iAO,iD)
     &                          * TabAO(1,iGrid,iAO)
                  Tau(iD,iGrid) = Tau(iD,iGrid)
     &                          + Grid_AO(5,iGrid,iAO,iD)
     &                          * TabAO(5,iGrid,iAO)
     &                          + Grid_AO(8,iGrid,iAO,iD)
     &                          * TabAO(8,iGrid,iAO)
     &                          + Grid_AO(10,iGrid,iAO,iD)
     &                          * TabAO(10,iGrid,iAO)
               End Do
            End Do
         End Do
      Else
         Call abend()
      End If
*
*     Generate the sigma vectors
*
      If (Allocated(Sigma)) Then
         If (nD.eq.1) Then
            Do iGrid=1, mGrid
               Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                       +GradRho(2,iGrid)**2
     &                       +GradRho(3,iGrid)**2
            End Do
         Else
            Do iGrid=1, mGrid
               Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                       +GradRho(2,iGrid)**2
     &                       +GradRho(3,iGrid)**2
               Sigma(2,iGrid)=GradRho(1,iGrid)*GradRho(4,iGrid)
     &                       +GradRho(2,iGrid)*GradRho(5,iGrid)
     &                       +GradRho(3,iGrid)*GradRho(6,iGrid)
               Sigma(3,iGrid)=GradRho(4,iGrid)**2
     &                       +GradRho(5,iGrid)**2
     &                       +GradRho(6,iGrid)**2
            End Do
         End If
      End If
#ifdef _DEBUGPRINT_
      Do iD = 1, nD
         Call RecPrt('Dens_AO',' ',Dens_AO(:,:,iD),nAO,nAO)
         Call RecPrt('Grid_AO',' ',Grid_AO(:,:,:,iD),mAO*mGrid,
     &                                               nAO)
      End Do
#endif

      Return
      End
