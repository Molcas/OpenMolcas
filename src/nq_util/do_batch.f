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
*               2021, Jie Bao                                          *
************************************************************************
      Subroutine Do_Batch(Kernel,Func,mGrid,
     &                    list_s,nlist_s,List_Exp,List_Bas,
     &                    Index,nIndex,
     &                    FckInt,nFckDim,nFckInt,
     &                    ipTabAO,mAO,nSym,
     &                    Dens,nDens,nD,
     &                    nP2_ontop,nShell,
     &                    Do_Mo,Do_TwoEl,l_Xhol,
     &                    TmpPUVX,nTmpPUVX,TabMO,TabSO,
     &                    nMOs,CMOs,nCMO,DoIt,
     &                    P2mo,P2unzip,np2act,D1mo,D1Unzip,nd1mo,
     &                    P2_ontop,
     &                    Do_Grad,Grad,nGrad,ndRho_dR,nGrad_Eff,
     &                    list_g,IndGrd,iTab,Temp,F_xc,dW_dR,iNQ,Maps2p,
     &                    DFTFOCK,LTEG_DB,
     &                    PDFTPot1,PDFTFocI,PDFTFocA)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use iSD_data
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use KSDFT_Info
      use nq_Grid, only: Grid, Weights, Rho, GradRho, Sigma, nRho
      use nq_Grid, only: vRho, vSigma, vTau, vLapl
      use nq_Grid, only: l_CASDFT, TabAO, TabAO_Pack, dRho_dR
      use nq_pdft
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "SysDef.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "debug.fh"
#include "ksdft.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "pamint.fh"
#include "grid_on_disk.fh"
      Integer list_s(2,nlist_s),List_Exp(nlist_s),DoIt(nMOs),
     &        ipTabAO(nlist_s+1,2), IndGrd(nGrad_Eff),
     &        list_g(3,nlist_s), iTab(4,nGrad_Eff), Index(nIndex),
     &        Maps2p(nShell,0:nSym-1), List_Bas(2,nlist_s)
      Real*8 A(3), RA(3), Grad(nGrad),
     &       FckInt(nFckInt,nFckDim), Dens(nDens,nD),
     &       TabMO(mAO,mGrid,nMOs),TabSO(mAO,mGrid,nMOs),
     &       CMOs(nCMO),P2mo(np2act),D1mo(nd1mo),
     &       P2_ontop(nP2_ontop,mGrid) , Temp(nGrad),
     &       F_xc(mGrid),
     &       dW_dR(nGrad_Eff,mGrid),
     &       PDFTPot1(nPot1),PDFTFocI(nPot1),PDFTFocA(nPot1)
      Real*8 TmpPUVX(nTmpPUVX)
      Logical Do_Grad,Do_Mo,Do_TwoEl,Unpack
      Logical l_Xhol, l_tanhr
      Character*4 DFTFOCK
      Integer nAOs
      Real*8 P2_ontop_d(nP2_ontop,nGrad_Eff,mGrid)
      Real*8,DIMENSION(:),ALLOCATABLE::P2MOCube,P2MOCubex,P2MOCubey,
     &                                 P2MOCubez,MOs,MOx,MOy,MOz
*     MOs,MOx,MOy and MOz are for active MOs.
*     MOas is for all MOs.
      Real*8,DIMENSION(NASHT4)::P2Unzip
      Real*8,DIMENSION(NASHT**2)::D1Unzip
      Integer LTEG_DB,nPMO3p
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Logical Debug_Save
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      iRout = 112
      iPrint = nPrint(iRout)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_Twoel.and.nTmpPUVX.eq.1) Then
         TmpPUVX(1)=0.0D0
         Write (6,*) DFTFOCK
         Call abend()
      End If
************************************************************************
*                                                                      *
      nTabAO=Size(TabAO)
#ifdef _DEBUGPRINT_
      Debug_Save=Debug
      Debug=Debug.or.iPrint.ge.99
*
      If (Debug) Then
         Write (6,*) ' In Do_Batch'
         Write (6,*) ' nRho=',nRho
         Write (6,*) 'Grid=',DDot_(mGrid,Grid(1,1),3,Grid(1,1),3),
     &                       DDot_(mGrid,Grid(2,1),3,Grid(2,1),3),
     &                       DDot_(mGrid,Grid(3,1),3,Grid(3,1),3),
     &                       DDot_(mGrid,Weights  ,1,Weights  ,1), mGrid
      End If
#endif
*                                                                      *
      mRho=-1
      l_tanhr=.false.

      CALL PDFTMemAlloc(mGrid,nOrbt)
      If (Functional_Type.eq.CASDFT_Type) Then
         mRho = nP2_ontop
      Else If(l_casdft) then !GLM
         mRho = nP2_ontop
      End If
*
      If (mRho.ne.-1) Then
         Call GetMem('Rho_I','Allo','Real',ipRhoI,mGrid*mRho)
         Call GetMem('Rho_A','Allo','Real',ipRhoA,mGrid*mRho)
         call dcopy_(mGrid*mRho,[Zero],0 ,Work(ipRhoI), 1)
         call dcopy_(mGrid*mRho,[Zero],0 ,Work(ipRhoA), 1)
      Else
         ipRhoI = ip_Dummy
         ipRhoA = ip_Dummy
      End If
************************************************************************
*                                                                      *
      ipSOS=0 ! dummy initialize
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the AOs on the grid points.                             *
*                                                                      *
************************************************************************
*
      TabAO(:,:,:)=Zero
      UnPack=.False.
      If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Use_Old .and.
     &      .Not.Do_Grad         .and.
     &    Functional_Type.eq.Old_Functional_Type)) Then
*
*------- Retrieve the AOs from disc
*
         Call iDaFile(Lu_Grid,2,ipTabAO,2*(nlist_s+1),iDisk_Grid)
         mTabAO=ipTabAO(nlist_s+1,2)-1
         Call dDaFile(Lu_Grid,2,TabAO,mTabAO,iDisk_Grid)
         Unpack=Packing.eq.On
*
      Else
*
*------- Generate the values of the AOs on the grid
*
         Call FZero(Work(ipMem),nMem)
         ipxyz=ipMem
*
         iOff = 1
         Do ilist_s=1,nlist_s
            ish=list_s(1,ilist_s)

            iShll = iSD( 0,iSh)
            iAng  = iSD( 1,iSh)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iPrim_Eff=List_Exp(ilist_s)
            iAO   = iSD( 7,iSh)
            mdci  = iSD(10,iSh)
            iShll = iSD(0,iSh)
            iCnttp= iSD(13,iSh)
            iCnt  = iSD(14,iSh)
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
            nDrv     = mRad - 1
            nForm    = 0
            Do iDrv  = 0, nDrv
               nForm = nForm + nElem(iDrv)
            End Do
            nTerm    = 2**nDrv
            nxyz     = mGrid*3*(iAng+mRad)
            nRadial  = iBas_Eff*mGrid*mRad
            ipRadial = ipxyz + nxyz
            ipAng_   = ipRadial + nRadial
            ipAng    = ip_of_iWork_d(Work(ipAng_))
*
            iR=list_s(2,ilist_s)
*
            ipx=iPhase(1,iR)
            ipy=iPhase(2,iR)
            ipz=iPhase(3,iR)
            px=DBLE(iPhase(1,iR))
            py=DBLE(iPhase(2,iR))
            pz=DBLE(iPhase(3,iR))
            RA(1) = px*A(1)
            RA(2) = py*A(2)
            RA(3) = pz*A(3)
            iSym=NrOpr(iR)
*
*---------- Evaluate AOs at RA
*
            ipTabAO(iList_s,1)=iOff
*                                                                      *
            Call AOEval(iAng,mGrid,Grid,Work(ipxyz),RA,
     &                  Shells(iShll)%Transf,
     &                  RSph(ipSph(iAng)),nElem(iAng),iCmp,
     &                  iWork(ipAng),nTerm,nForm,T_X,mRad,
     &                  iPrim,iPrim_Eff,Shells(iShll)%Exp,
     &                  Work(ipRadial),
     &                  iBas_Eff,
     &                  Shells(iShll)%pCff(1,iBas-iBas_Eff+1),
     &                  TabAO_Pack(iOff:),
     &                  mAO,px,py,pz,ipx,ipy,ipz)
            iOff = iOff + mAO*mGrid*iBas_Eff*iCmp
*
         End Do
         ipTabAO(nList_s+1,1)=iOff
*
*        AOs are packed and written to disk.
*
         If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Regenerate .and.
     &       .Not.Do_Grad)) Then
*
            If (Packing.eq.On) Then
               Unpack=.True.
*
*------------- Pack before they are put on disc
*
               jOff = 1
               Do ilist_s=1,nlist_s
                  ish=list_s(1,ilist_s)
                  iCmp  = iSD( 2,iSh)
                  iBas_Eff = List_Bas(1,ilist_s)
                  nData=mAO*mGrid*iBas_Eff*iCmp
*
*                 Check if we should store any AOs at all!
*
                  iOff = ipTabAO(ilist_s,1)
                  If (nData.gt.nTmp) Then
                     Call WarningMessage(2,'nData.gt.nTmp')
                     Call Abend()
                  End If
                  call dcopy_(nData,TabAO_Pack(iOff:),1,
     &                              Work(ipTmp),1)
                  Call PkR8(0,nData,nByte,Work(ipTmp),
     &                                    TabAO_Pack(jOff:))
                  mData = (nByte+RtoB-1)/RtoB
                  If (mData.gt.nData) Then
                     Call WarningMessage(2,'mData.gt.nData')
                     Write (6,*) 'nData=',nData
                     Write (6,*) 'nData=',nData
                     Call Abend()
                  End If
                  ipTabAO(iList_s,2)=nByte
                  jOff = jOff + mData
               End Do
               ipTabAO(nList_s+1,2)=jOff
            Else
               ipTabAO(nList_s+1,2)=ipTabAO(nList_s+1,1)
            End If
*
            Call iDaFile(Lu_Grid,1,ipTabAO,2*(nlist_s+1),iDisk_Grid)
            mTabAO=ipTabAO(nList_s+1,2)-1
            Call dDaFile(Lu_Grid,1,TabAO,mTabAO,iDisk_Grid)
*
         End If
*
      End If
*
*---- Unpack AOs
*
      If (Unpack) Then
*
         jOff = ipTabAO(nlist_s+1,2)
         Do ilist_s=nlist_s,1,-1
            ish=list_s(1,ilist_s)
            iCmp  = iSD( 2,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            nData=mAO*mGrid*iBas_Eff*iCmp
            nByte=ipTabAO(ilist_s,2)
*
            iOff = ipTabAO(ilist_s,1)
            If (nByte.gt.0) Then
               mData = (nByte+RtoB-1)/RtoB
               jOff = jOff - mData
               If (mData.gt.nTmp) Then
                  Call WarningMessage(2,'mData.gt.nTmp')
                  Call Abend()
               End If
               Call UpkR8(0,nData,nByte,TabAO_Pack(jOff:),Work(ipTmp))
               call dcopy_(nData,Work(ipTmp),1,TabAO_Pack(iOff:),1)
            Else
               mData=0
               TabAO_Pack(1:nData)=Zero
            End If
         End Do
*
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Evaluate some MOs on the grid                                    *
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_MO) Then
         Call FZero(TabMO,mAO*mGrid*nMOs)
         Call FZero(TabSO,mAO*mGrid*nMOs)
*
         Do ilist_s=1,nlist_s
            ish=list_s(1,ilist_s)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iAO   = iSD( 7,iSh)
            mdci  = iSD(10,iSh)
*
*---------- Allocate memory for SO and MO
*
            kAO   = iCmp*iBas*mGrid
            nDeg  = nSym/dc(mdci)%nStab
            nSO   = kAO*nDeg*mAO
            ipSOs = ipMem
            Call FZero(Work(ipSOs),nSO)
*
            iR=list_s(2,ilist_s)
            iSym=NrOpr(iR)
*
*---------- Distribute contributions of AOs if this particular shell
*           on to the SOs of this shell. The SOs are only stored
*           temporarily!
*
            Call SOAdpt_NQ(TabAO_Pack(ipTabAO(iList_s,1):),mAO,mGrid,
     &                     iBas,iBas_Eff,iCmp,iSym,Work(ipSOs),nDeg,
     &                     iAO)
*
            Call GetMem('TmpCM','Allo','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Allo','Inte',ipTDoIt,nMOs)
            Call  SODist2(Work(ipSOs),mAO,mGrid,iBas,
     &                   iCmp,nDeg,TabSO,
     &                   nMOs,iAO,Work(ipTmpCMO),
     &                   nCMO,iWork(ipTDoIt))
            Call GetMem('TmpCM','Free','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Free','Inte',ipTDoIt,nMOs)
*
            Call  SODist(Work(ipSOs),mAO,mGrid,iBas,iCmp,nDeg,TabMO,
     &                  nMOs,iAO,CMOs,nCMO,DoIt)
*
         End Do
      End If
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*---- Compute Rho, Grad Rho, Tau, Laplacian, and the Sigma vectors.
*     In case of gradient calculations compute Cartesian derivatives
*     of Rho, Grad Rho, Tau, and the Laplacian.
*                                                                      *
      Call Mk_Rho(list_s,nlist_s,Work(ip_Fact),ndc,list_bas,
     &            Index,nIndex,list_g,Do_Grad)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      If (l_casdft) then
         T_Rho=T_X*1.0D-4
         Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
         Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
         Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)
      End If

      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
**************************************************************************
* TLSDA,TLSDA5                                                           *
**************************************************************************
       If(KSDFA(1:5).eq.'TLSDA'.or.KSDFA(1:6).eq.'TLSDA5') then !GLM
        if(Debug) write(6,*) 'in do_batch.f for TLSDA option'

       nPMO3p=1
       IF(lft.and.lGGA) THEN
        nPMO3p=mGrid*NASHT
       END IF

       CALL mma_allocate(P2MOCube,mGrid*NASHT)
       CALL mma_allocate(P2MOCubex,nPMO3p)
       CALL mma_allocate(P2MOCubey,nPMO3p)
       CALL mma_allocate(P2MOCubez,nPMO3p)
       CALL mma_allocate(MOs,mGrid*NASHT)
       CALL mma_allocate(MOx,mGrid*NASHT)
       CALL mma_allocate(MOy,mGrid*NASHT)
       CALL mma_allocate(MOz,mGrid*NASHT)

       CALL CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,
     &                   MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,
     &                   Do_Grad)
       Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,
     &               nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &               Work(ipRhoA),mRho,Do_Grad,
     &               P2MOCube,MOs,MOx,MOy,MOz)
       Else !AO-based run for gradients
!        nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
        P2_ontop_d(:,:,:) = 0
        !Determine number of AOs:
        nAOs = nMOs
        Call  Do_Pi2grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &                   P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &                   list_s,nlist_s,list_bas,Index,nIndex,
     &                   D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &                   Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,
     &                   nAOs,nCMO,TabSO,nsym,lft,
     &                   P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                   nPMO3p,MOs,MOx,MOy,MOz)
       End If
       CALL TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,
     &                     l_tanhr,nRho,mGrid,nP2_OnTop,
     &                     ndRho_dR,nGrad_Eff,Do_Grad)
       CALL mma_deallocate(P2MOCube)
       CALL mma_deallocate(P2MOCubex)
       CALL mma_deallocate(P2MOCubey)
       CALL mma_deallocate(P2MOCubez)
       CALL mma_deallocate(MOs)
       CALL mma_deallocate(MOx)
       CALL mma_deallocate(MOy)
       CALL mma_deallocate(MOz)
**************************************************************************
       End if !tlsda

cRKCft
************************************************************************
* FTLSDA                                                               *
************************************************************************

       If(KSDFA(1:6).eq.'FTLSDA') then !GLM

       nPMO3p=1
       IF(lft.and.lGGA) THEN
        nPMO3p=mGrid*NASHT
       END IF
       CALL mma_allocate(P2MOCube,mGrid*NASHT)
       CALL mma_allocate(P2MOCubex,nPMO3p)
       CALL mma_allocate(P2MOCubey,nPMO3p)
       CALL mma_allocate(P2MOCubez,nPMO3p)
       CALL mma_allocate(MOs,mGrid*NASHT)
       CALL mma_allocate(MOx,mGrid*NASHT)
       CALL mma_allocate(MOy,mGrid*NASHT)
       CALL mma_allocate(MOz,mGrid*NASHT)

       CALL CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,
     &                   MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,
     &                   Do_Grad)
       Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,
     &               nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &               Work(ipRhoA),mRho,Do_Grad,
     &               P2MOCube,MOs,MOx,MOy,MOz)
       Else !AO-based run for gradients
        P2_ontop_d(:,:,:) = 0
        !Determine number of AOs:
        nAOs = nMOs
        Call  Do_Pi2grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &                   P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &                   list_s,nlist_s,list_bas,Index,nIndex,
     &                   D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &                   Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,
     &                   nAOs,nCMO,TabSO,nsym,lft,
     &                   P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                   nPMO3p,MOs,MOx,MOy,MOz)
       End If
       CALL TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,
     &                    l_tanhr,nRho,mGrid,nP2_OnTop,
     &                    ndRho_dR,nGrad_Eff,Do_Grad)
       CALL mma_deallocate(P2MOCube)
       CALL mma_deallocate(P2MOCubex)
       CALL mma_deallocate(P2MOCubey)
       CALL mma_deallocate(P2MOCubez)
       CALL mma_deallocate(MOs)
       CALL mma_deallocate(MOx)
       CALL mma_deallocate(MOy)
       CALL mma_deallocate(MOz)

       End if


*      ^ end if for GLM stuff
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*======================================================================*
************************************************************************
* TBLYP,TPBE,TREVPBE,TSSBSW,TSSBD,TS12G                                *
************************************************************************
      if(KSDFA(1:5).eq.'TBLYP'.or. !GLM
     &   KSDFA(1:4).eq.'TPBE'.or.
     &   KSDFA(1:6).eq.'TSSBSW'.or.
     &   KSDFA(1:5).eq.'TSSBD'.or.
     &   KSDFA(1:5).eq.'TS12G'.or.
     &   KSDFA(1:5).eq.'TOPBE'.or.
     &   KSDFA(1:7).eq.'TREVPBE') then

       nPMO3p=1
       IF(lft.and.lGGA) THEN
        nPMO3p=mGrid*NASHT
       END IF
       CALL mma_allocate(P2MOCube,mGrid*NASHT)
       CALL mma_allocate(P2MOCubex,nPMO3p)
       CALL mma_allocate(P2MOCubey,nPMO3p)
       CALL mma_allocate(P2MOCubez,nPMO3p)
       CALL mma_allocate(MOs,mGrid*NASHT)
       CALL mma_allocate(MOx,mGrid*NASHT)
       CALL mma_allocate(MOy,mGrid*NASHT)
       CALL mma_allocate(MOz,mGrid*NASHT)

       CALL CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,
     &                   MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,
     &                   Do_Grad)
       Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,
     &               nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &               Work(ipRhoA),mRho,Do_Grad,
     &               P2MOCube,MOs,MOx,MOy,MOz)
       Else !AO-based run for gradients
!        nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
        P2_ontop_d(:,:,:) = 0
        !Determine number of AOs:
        nAOs = nMOs
        Call  Do_Pi2grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &                   P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &                   list_s,nlist_s,list_bas,Index,nIndex,
     &                   D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &                   Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,
     &                   nAOs,nCMO,TabSO,nsym,lft,
     &                   P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                   nPMO3p,MOs,MOx,MOy,MOz)
       End If

        CALL TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,
     &                     l_tanhr,nRho,mGrid,nP2_OnTop,
     &                     ndRho_dR,nGrad_Eff,Do_Grad)
       CALL mma_deallocate(P2MOCube)
       CALL mma_deallocate(P2MOCubex)
       CALL mma_deallocate(P2MOCubey)
       CALL mma_deallocate(P2MOCubez)
       CALL mma_deallocate(MOs)
       CALL mma_deallocate(MOx)
       CALL mma_deallocate(MOy)
       CALL mma_deallocate(MOz)
      end if
*======================================================================*
*======================================================================*
************************************************************************
* FTBLYP,FTPBE,FTREVPBE                                                *
************************************************************************
      if(KSDFA(1:6).eq.'FTBLYP'.or. !GLM
     &   KSDFA(1:5).eq.'FTPBE'.or.
     &   KSDFA(1:6).eq.'FTOPBE'.or.
     &   KSDFA(1:8).eq.'FTREVPBE') then
*  *
       nPMO3p=1
       IF((lft.and.lGGA).and.Do_Grad) THEN
        nPMO3p=mGrid*NASHT
       END IF
       CALL mma_allocate(P2MOCube,mGrid*NASHT)
       CALL mma_allocate(P2MOCubex,nPMO3p)
       CALL mma_allocate(P2MOCubey,nPMO3p)
       CALL mma_allocate(P2MOCubez,nPMO3p)
       CALL mma_allocate(MOs,mGrid*NASHT)
       CALL mma_allocate(MOx,mGrid*NASHT)
       CALL mma_allocate(MOy,mGrid*NASHT)
       CALL mma_allocate(MOz,mGrid*NASHT)

       CALL CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,
     &                   MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,
     &                   Do_Grad)
       Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,
     &               nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &               Work(ipRhoA),mRho,Do_Grad,
     &               P2MOCube,MOs,MOx,MOy,MOz)
       Else !AO-based run for gradients
!        nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
        P2_ontop_d(:,:,:) = 0
        !Determine number of AOs:
        nAOs = nMOs
        Call  Do_Pi2grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &                   P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &                   list_s,nlist_s,list_bas,Index,nIndex,
     &                   D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &                   Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,
     &                   nAOs,nCMO,TabSO,nsym,lft,
     &                   P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                   nPMO3p,MOs,MOx,MOy,MOz)

       End If
*
          CALL TranslateDens(P2_OnTop,dRho_dr,P2_OnTop_d,
     &                       l_tanhr,nRho,mGrid,nP2_OnTop,
     &                       ndRho_dR,nGrad_Eff,Do_Grad)

       CALL mma_deallocate(P2MOCube)
       CALL mma_deallocate(P2MOCubex)
       CALL mma_deallocate(P2MOCubey)
       CALL mma_deallocate(P2MOCubez)
       CALL mma_deallocate(MOs)
       CALL mma_deallocate(MOx)
       CALL mma_deallocate(MOy)
       CALL mma_deallocate(MOz)
       end if
*======================================================================*
*======================================================================*

      If (l_casdft) Then
      If (nD.eq.1) Then
         Do iGrid=1, mGrid
            Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                    +GradRho(2,iGrid)**2
     &                    +GradRho(3,iGrid)**2
         End Do
      Else
         Do iGrid=1, mGrid
            Sigma(1,iGrid)=GradRho(1,iGrid)**2
     &                    +GradRho(2,iGrid)**2
     &                    +GradRho(3,iGrid)**2
            Sigma(2,iGrid)=GradRho(1,iGrid)*GradRho(4,iGrid)
     &                    +GradRho(2,iGrid)*GradRho(5,iGrid)
     &                    +GradRho(3,iGrid)*GradRho(6,iGrid)
            Sigma(3,iGrid)=GradRho(4,iGrid)**2
     &                    +GradRho(5,iGrid)**2
     &                    +GradRho(6,iGrid)**2
         End Do
      End If
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.CASDFT_type) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*------- Compute P2_OnTop at the grid
         Call Do_P2new(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*     Integrate out the number of electrons
*
      if(l_casdft) then
        T_Rho=T_X*1.0D-4
        Dens_t2=Dens_t2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
        Dens_a2=Dens_a2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
        Dens_b2=Dens_b2+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)
      end if
*
*     Integrate out the number of electrons, |grad|, and tau
*
      T_Rho=T_X*1.0D-4
      If (Functional_type.eq.LDA_Type) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.GGA_type) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type1) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type2) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,nD,T_Rho)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-- (A.Ohrn): Here I add the routine which constructs the kernel for
*   the Xhole application. A bit 'cheating' but hey what da hey!
*
      If(l_Xhol) then
#ifdef _NOT_USED_TESTED_OR_MAINTAINED_
        Call Xhole(nRho,mGrid,Rho,Grid,mAO,nMOs,TabMO,ndF_dRho,nD,
     &             dF_dRho,Weights,ip_OrbDip,Func)
#endif
        Go To 1979
      Endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Evaluate the functional on the grid                              *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      vRho(:,1:mGrid)=Zero
      If (Allocated(vSigma)) vSigma(:,1:mGrid)=Zero
      If (Allocated(vTau)) vTau(:,1:mGrid)=Zero
      If (Allocated(vLapl)) vLapl(:,1:mGrid)=Zero
      F_xc(1:mGrid)=Zero
      If (l_casdft) Then
         F_xca(1:mGrid)=Zero
         F_xcb(1:mGrid)=Zero
         tmpB(1:mGrid)=Zero
      End If
*
*1)   evaluate the energy density, the derivative of the functional with
*     respect to rho and grad rho.
*
      Call Kernel(mGrid,nD)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*     Integrate the energy of the functional
*
      Func=Func+DDot_(mGrid,Weights,1,F_xc,1)
      If (l_casdft) Then
         Funcaa=Funcaa+DDot_(mGrid,Weights,1,F_xca,1)
         Funcbb=Funcbb+DDot_(mGrid,Weights,1,F_xcb,1)
         Funccc=Funccc+DDot_(mGrid,Weights,1,tmpB,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
1979  Continue  !Jump here and skip the call to the kernel.
      IF(do_pdftPot) THEN
       CALL mma_allocate(MOs ,mGrid*NASHT)
      END IF

      If (.Not.Do_Grad) Then
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Compute the DFT contribution to the Fock matrix                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         If (Functional_type.eq.LDA_Type) Then
*
*2)         form matrix elements for the potential, from derivative of
*           the functional with  respect to rho, respectively.
*
!First, calculate some sizes:
!             NFINT=nTmpPUVX
!             NTOT1=nFckInt

             If(KSDFA(1:5).eq.'TLSDA') then
               If(do_pdftPot) then
               lft=.false.
               lGGA=.false.

               CALL TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
               IF(lft.and.lGGA) THEN
                CALL TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
                CALL TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
                CALL TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
               END IF
               CALL TransActMO(MOs, TabMO,mAO,mGrid,nMOs)

               Call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,
     &                        nP2_ontop,MOas)
               Call Calc_Pot2(Work(LTEG_DB),mGrid,P2_ontop,nP2_ontop)
               Call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
              end if
             else If(KSDFA(1:6).eq.'FTLSDA') then
               If(do_pdftPot) then
               lft=.true.
               lGGA=.false.

               CALL TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
               IF(lft.and.lGGA) THEN
                CALL TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
                CALL TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
                CALL TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
               END IF
               CALL TransActMO(MOs, TabMO,mAO,mGrid,nMOs)

               Call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,
     &                        nP2_ontop,MOas)
               Call Calc_Pot2(Work(LTEG_DB),mGrid,P2_ontop,nP2_ontop)
               Call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
              end if
             end if

             If(KSDFA(1:5).ne.'TLSDA'.and.KSDFA(1:6).ne.'FTLSDA') then
                 Call DFT_Int(list_s,nlist_s,FckInt,nFckInt,
     &                        nD,Work(ip_Fact),ndc,list_bas)
             End If
*                                                                      *
************************************************************************
*                                                                      *
        Else If (Functional_type.eq.GGA_type   .or.
     &           Functional_type.eq.CASDFT_type     ) Then     ! CGG
*
*2)        form matrix elements for the potential, from derivative of
*          the functional with  respect to rho, respectively.
*
*          and
*
*3)        form contributions to the matrix elements of the potenial
*          from the derivative of the functional with respect to grad
*          rho.
*
             If(KSDFA(1:4).eq.'TPBE'.or.
     &               KSDFA(1:5).eq.'TOPBE'.or.
     &               KSDFA(1:5).eq.'TBLYP'.or.
     &               KSDFA(1:7).eq.'TREVPBE') then

              If(do_pdftPot) then
               lft=.false.
               lGGA=.true.

               CALL TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
               IF(lft.and.lGGA) THEN
                CALL TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
                CALL TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
                CALL TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
               END IF
               CALL TransActMO(MOs, TabMO,mAO,mGrid,nMOs)

               Call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,
     &                        nP2_ontop,MOas)

               Call Calc_Pot2(Work(LTEG_DB),mGrid,P2_ontop,nP2_ontop)

               Call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
              end if
             Else If(KSDFA(1:5).eq.'FTPBE'.or.
     &               KSDFA(1:6).eq.'FTOPBE'.or.
     &               KSDFA(1:6).eq.'FTBLYP'.or.
     &               KSDFA(1:8).eq.'FTREVPBE') then
               If(do_pdftPot) then
               lft=.true.
               lGGA=.true.

               CALL TransferMO(MOas,TabMO,mAO,mGrid,nMOs,1)
               IF(lft.and.lGGA) THEN
                CALL TransferMO(MOax,TabMO,mAO,mGrid,nMOs,2)
                CALL TransferMO(MOay,TabMO,mAO,mGrid,nMOs,3)
                CALL TransferMO(MOaz,TabMO,mAO,mGrid,nMOs,4)
               END IF
               CALL TransActMO(MOs, TabMO,mAO,mGrid,nMOs)

               Call Calc_Pot1(PDFTPot1,TabMO,mAO,mGrid,nMOs,P2_ontop,
     &                        nP2_ontop,MOas)
               Call Calc_Pot2(Work(LTEG_DB),mGrid,P2_ontop,nP2_ontop)
               Call PDFTFock(PDFTFocI,PDFTFocA,D1Unzip,mGrid,MOs)
              end if
             end if
             If(.not.l_casdft) then
               Call DFT_Int(list_s,nlist_s,FckInt,nFckInt,
     &                      nD,Work(ip_Fact),ndc,list_bas)
             end if
*                                                                      *
************************************************************************
*                                                                      *
        Else If (Functional_type.eq.meta_GGA_type1 .or.
     &           Functional_Type.eq.meta_GGA_type2) Then
*
*2)        form matrix elements for the potential, from derivative of
*          the functional with  respect to rho, respectively.
*
*          and
*
*3)        form contributions to the matrix elements of the potenial
*          from the derivative of the functional with respect to grad
*          rho.
*
*         and
*
*4)       form contributions to the matrix elements of the potential
*         from the derivatives of the functional with respect to nabla
*         rho and/or tau.
*
          Call DFT_Int(list_s,nlist_s,FckInt,nFckInt,
     &                 nD,Work(ip_Fact),ndc,list_bas)
*                                                                      *
************************************************************************
*                                                                      *
         Else
            Call WarningMessage(2,'Wrong type of functional')
            Call Abend()
         End If    !  Functional Type
*                                                                      *
************************************************************************
*                                                                      *
*    Compute the DFT contribution to the gradient                      *
*                                                                      *
************************************************************************
*                                                                      *
      Else
*
         Call DFT_Grad(Grad,nGrad,nD,Grid,mGrid,
     &                 dRho_dR,ndRho_dR,nGrad_Eff,IndGrd,
     &                 Weights,iTab,Temp,F_xc,dW_dR,iNQ)
*
      End If
*                                                                      *
      If(do_pdftPot) then
       CALL mma_deallocate(MOs   )
      End If
************************************************************************
*                                                                      *
      CALL PDFTMemDeAlloc()
      If(mRho.ne.-1) Then
        Call GetMem('Rho_I','Free','Real',ipRhoI,mGrid*mRho)
        Call GetMem('Rho_A','Free','Real',ipRhoA,mGrid*mRho)
      End If
#ifdef _DEBUGPRINT_
      Debug=Debug_Save
#endif
      Return
* Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(Maps2p)
      If (.False.) Call Unused_real_array(Dens)
      End
