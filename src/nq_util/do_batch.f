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
      Subroutine Do_Batch(Kernel,Func,
     &                    Grid,Weights,Rho,mGrid,nRho,
     &                    list_s,nlist_s,List_Exp,List_Bas,
     &                    Index,nIndex,AOInt,nAOInt,
     &                    FckInt,nFckDim,nFckInt,SOTemp,nSOTemp,
     &                    TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                    Dens,nDens,nD,
     &                    ndF_dRho,nP2_ontop,ndF_dP2ontop,nShell,
     &                    Do_Mo,Do_TwoEl,l_Xhol,
     &                    TmpPUVX,nTmpPUVX,TabMO,TabSO,
     &                    nMOs,CMOs,nCMO,DoIt,
     &                    P2mo,np2act,D1mo,nd1mo,P2_ontop,
     &                    Do_Grad,Grad,nGrad,dRho_dR,ndRho_dR,nGrad_Eff,
     &                    list_g,IndGrd,iTab,Temp,F_xc,dW_dR,iNQ,Maps2p,
     &                    dF_dRho,dF_dP2ontop,DFTFOCK,LOE_DB,LTEG_DB)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: SubBlock                                                *
*                                                                      *
* Calling    : QEnter                                                  *
*              AOEval                                                  *
*              Do_Rho_*                                                *
*              Kernel                                                  *
*              DFT_Int                                                 *
*              Do_DFT_Grad                                             *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use iSD_data
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "SysDef.fh"
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
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
      Real*8 A(3), RA(3), Grid(3,mGrid), Weights(mGrid),
     &       Rho(nRho,mGrid), AOInt(nAOInt*nAOInt,nD),
     &       dF_dRho(ndF_dRho,mGrid), TabAO(nTabAO), Grad(nGrad),
     &       FckInt(nFckInt,nFckDim),
     &       Dens(nDens,nD), SOTemp(nSOTemp,nD),
     &       TabMO(mAO,mGrid,nMOs),TabSO(mAO,mGrid,nMOs),
     &       CMOs(nCMO),P2mo(np2act),D1mo(nd1mo),
     &       P2_ontop(nP2_ontop,mGrid) , Temp(nGrad),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff), F_xc(mGrid),
     &       dW_dR(nGrad_Eff,mGrid),dF_dP2ontop(ndF_dP2ontop,mGrid)
      Real*8 TmpPUVX(nTmpPUVX)
      Logical Do_Grad,Do_Mo,Do_TwoEl,Unpack
      Logical l_Xhol, l_tanhr,l_casdft
      Logical ft
      Character*4 DFTFOCK
      Integer dindex
      Real*8 dTot_d,ratio_d,Zeta_d
      Integer nAOs
      Real*8 P2_ontop_d(nP2_ontop,nGrad_Eff,mGrid)
      Integer ntot1
      Integer LOE_DB,LTEG_DB
*define _DEBUG_
#ifdef _DEBUG_
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
      iRout = 112
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*        Initializing MCPDFT global variable                           *
************************************************************************
*                                                                      *
      l_casdft = KSDFA(1:5).eq.'TLSDA'   .or.
     &           KSDFA(1:6).eq.'TLSDA5'  .or.
     &           KSDFA(1:5).eq.'TBLYP'   .or.
     &           KSDFA(1:6).eq.'TSSBSW'  .or.
     &           KSDFA(1:5).eq.'TSSBD'  .or.
     &           KSDFA(1:5).eq.'TS12G'  .or.
     &           KSDFA(1:4).eq.'TPBE'    .or.
     &           KSDFA(1:5).eq.'FTPBE'   .or.
     &           KSDFA(1:7).eq.'TREVPBE' .or.
     &           KSDFA(1:8).eq.'FTREVPBE'.or.
     &           KSDFA(1:6).eq.'FTLSDA'  .or.
     &           KSDFA(1:5).eq.'TOPBE'   .or.
     &           KSDFA(1:6).eq.'FTOPBE'  .or.
     &           KSDFA(1:6).eq.'FTBLYP'
************************************************************************
#ifdef _TIME_
      Call qEnter('Do_Batch ')
#endif
#ifdef _DEBUG_
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
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      thrspi = 1.0d-15
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
      l_tanhr=.false.

      If (Functional_Type.eq.CASDFT_Type) Then
         mRho = nP2_ontop
      Else If(DFTFOCK.eq.'DIFF'.and.nD.eq.2) Then
         mRho = nRho/nD
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
#ifdef _TIME_
      Call QEnter('AO')
#endif
*
      Call FZero(TabAO,nTabAO)
c      write(6,*) 'nTabAO value in do_batch.f=', nTabAO
      UnPack=.False.
      If (NQ_Direct.eq.Off .and. (Grid_Status.eq.Use_Old .and.
     &      .Not.Do_Grad         .and.
     &    Functional_Type.eq.Old_Functional_Type)) Then
C        Write (*,*) 'Read from disk!'
*
*------- Retrieve the AOs from disc
*
C        Write (*,*) 'iDisk_Grid=',iDisk_Grid
         Call iDaFile(Lu_Grid,2,ipTabAO,2*(nlist_s+1),iDisk_Grid)
C        Write (*,*) 'ipTabAO(*,1)=',(ipTabAO(i,1),i=1,nlist_s+1)
C        Write (*,*) 'ipTabAO(*,2)=',(ipTabAO(i,2),i=1,nlist_s+1)
         mTabAO=ipTabAO(nlist_s+1,2)-1
C        Write (*,*) 'iDisk_Grid,mTabAO=',iDisk_Grid,mTabAO
         Call dDaFile(Lu_Grid,2,TabAO,mTabAO,iDisk_Grid)
C        Call RecPrt('TabAO from disk',' ',TabAO,1,mTabAO)
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

            iAng  = iSD( 1,iSh)
            iCmp  = iSD( 2,iSh)
            iBas  = iSD( 3,iSh)
            iBas_Eff = List_Bas(1,ilist_s)
            iPrim = iSD( 5,iSh)
            iPrim_Eff=List_Exp(ilist_s)
            iCff  = iSD( 4,iSh)
            iCff_Eff=iCff+iPrim*(iBas-iBas_Eff)
            iExp  = iSD( 6,iSh)
            iAO   = iSD( 7,iSh)
            ixyz  = iSD( 8,iSh)
            mdci  = iSD(10,iSh)
            iShell= iSD(11,iSh)
            iShll = iSD(0,iSh)
            call dcopy_(3,Work(ixyz),1,A,1)
*
            nDrv     = mRad - 1
            nForm    = 0
            Do iDrv  = 0, nDrv
               nForm = nForm + nElem(iDrv)
            End Do
            nTerm    = 2**nDrv
            nAngular = 5*nForm*nTerm
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
            iSym=NrOpr(iR,iOper,nSym)
#ifdef _DEBUG_
            If (debug) Write (6,*) 'mAO=',mAO
            If (iPrim_Eff.le.0 .or. iPrim_Eff.gt.iPrim) Then
               Call WarningMessage(2,'Do_batch: error in iPrim_Eff!')
               Write (6,*) '          iPrim=',iPrim
               Write (6,*) '          iPrim_Eff=',iPrim_Eff
               Call Abend()
            End If
#endif
*
*---------- Evaluate AOs at RA
*
#ifdef _DEBUG_
            iPrint_132=nPrint(132)
            If (Debug) nPrint(132)=99
#endif
            ipTabAO(iList_s,1)=iOff
*                                                                      *
c            write(6,*) 'iOff =', iOff
            Call AOEval(iAng,mGrid,Grid,Work(ipxyz),RA,
     &                  Transf(iShll),
     &                  RSph(ipSph(iAng)),nElem(iAng),iCmp,
     &                  iWork(ipAng),nTerm,nForm,T_X,mRad,
     &                  iPrim,iPrim_Eff,Work(iExp),Work(ipRadial),
     &                  iBas_Eff,Work(iCff_Eff),TabAO(iOff),
     &                  mAO,px,py,pz,ipx,ipy,ipz)
            iOff = iOff + mAO*mGrid*iBas_Eff*iCmp
#ifdef _DEBUG_
            nPrint(132)=iPrint_132
#endif
*
         End Do
*                write(6,*) 'AOValue in do_batch:'
*         write(6,'(10G20.10)') (TabAO(i),i=1,nTabAO)
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
                  iBas  = iSD( 3,iSh)
                  iBas_Eff = List_Bas(1,ilist_s)
                  nData=mAO*mGrid*iBas_Eff*iCmp
*
*                 Check if we should store any AO's at all!
*
                  iOff = ipTabAO(ilist_s,1)
                  ix=iDAMax_(nData,TabAO(iOff),1)
                  AOMax=Abs(TabAO(iOff-1+ix))
                  If (AOMax.ge.T_X) Then
                     If (nData.gt.nTmp) Then
                        Call WarningMessage(2,'nData.gt.nTmp')
                        Call Abend()
                     End If
                     call dcopy_(nData,TabAO(iOff),1,Work(ipTmp),1)
                     Call PkR8(0,nData,nByte,Work(ipTmp),TabAO(jOff))
                     mData = (nByte+RtoB-1)/RtoB
                     If (mData.gt.nData) Then
                        Call WarningMessage(2,'mData.gt.nData')
                        Write (6,*) 'nData=',nData
                        Write (6,*) 'nData=',nData
                        Call Abend()
                     End If
                  Else
                     nByte=0
                     mData = 0
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
            iBas  = iSD( 3,iSh)
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
               Call UpkR8(0,nData,nByte,TabAO(jOff),Work(ipTmp))
               call dcopy_(nData,Work(ipTmp),1,TabAO(iOff),1)
            Else
               mData=0
               Call FZero(TabAO(iOff),nData)
            End If
         End Do
*
      End If
*                                                                      *
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
            iShell= iSD(11,iSh)
*
*---------- Allocate memory for SO and MO
*
cGLM            kAO   = iCmp*iBas_Eff*mGrid
            kAO   = iCmp*iBas*mGrid
            nDeg  = nSym/nStab(mdci)
            nSO   = kAO*nDeg*mAO
            ipSOs = ipMem
            Call FZero(Work(ipSOs),nSO)
*
            iR=list_s(2,ilist_s)
            iSym=NrOpr(iR,iOper,nSym)
*
*---------- Distribute contributions of AOs if this particular shell
*           on to the SOs of this shell. The SOs are only stored
*           temporarily!
*
            Call SOAdpt_NQ(TabAO(ipTabAO(iList_s,1)),mAO,mGrid,iBas,
     &                  iBas_Eff,iCmp,iSym,Work(ipSOs),nDeg,iShell)
*
            Call GetMem('TmpCM','Allo','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Allo','Inte',ipTDoIt,nMOs)
            Call  SODist2(Work(ipSOs),mAO,mGrid,iBas,
     &                   iCmp,nDeg,TabSO,
     &                   iShell,nMOs,iAO,Work(ipTmpCMO),
     &                   nCMO,iWork(ipTDoIt))
            Call GetMem('TmpCM','Free','Real',ipTmpCMO,nCMO)
            Call GetMem('TDoIt','Free','Inte',ipTDoIt,nMOs)
*
            Call  SODist(Work(ipSOs),mAO,mGrid,iBas,iCmp,nDeg,TabMO,
     &                  iShell,nMOs,iAO,CMOs,nCMO,DoIt)
*
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QExit('AO')
      Call QEnter('Rho')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute density and grad_density
*                                                                      *
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
*
         Call Rho_LDA(Dens,nDens,nD,Rho,nRho,mGrid,
     &                list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                Work(ip_Fact),ndc,Work(ipTabAOMax),
     &                list_bas,Index,nIndex)
*
         If (Do_Grad)
     &      Call dRho_dR_LDA(Dens,nDens,nD,dRho_dR,ndRho_dr,
     &                       mGrid,list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                       nGrad_Eff,list_g,Maps2p,
     &                       nShell,Grid_Type,Fixed_Grid,
     &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
     &                       list_bas,Index,nIndex)

**************************************************************************
* TLSDA,TLSDA5                                                           *
**************************************************************************
       If(KSDFA(1:5).eq.'TLSDA'.or.KSDFA(1:6).eq.'TLSDA5') then !GLM
        if(Debug) write(6,*) 'in do_batch.f for TLSDA option'

**************************************************************************
* Comp_d is a function to integrate densities or gradients or whatever...*
* as long as it is worth integrating it.. to the limit of a horse! hehe  *
* iSwitch (last entry) will dictate what quantity will be integrated:    *
*  iSwitch = 0  total density                                            *
*  iSwitch = 1  alpha density                                            *
*  iSwitch = 2  beta density                                             *
*  ......                                                                *
*  iSwitch = 10 the horse!                                               *
**************************************************************************
        T_Rho=T_X*1.0D-4
        Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
        Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
        Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)
**************************************************************************
         Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run

         Call Do_P2glm(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho,Do_Grad)

       Else !AO-based run for gradients
         nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
         P2_ontop_d(:,:,:) = 0
         !Determine number of AOs:
         nAOs = nMOs
      ft=.false.
        Call  Do_P2GLM_grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &          list_s,nlist_s,list_bas,Index,nIndex,
     &          P2mo,np2act,D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &      Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,nAOs,nCMO,
     &      TabSO,nsym,ft)
       End If

       If(.not.Do_Grad) then
         do iGrid=0,mGrid-1
          dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
          ratio = 0.0d0
            write(LuMT,'(3(F10.6,A),5(F17.10,A))')
     &       Grid(1,iGrid+1),',',
     &       Grid(2,iGrid+1),',',
     &       Grid(3,iGrid+1),',',
     &       Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &       Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &       dTot*Weights(iGrid+1),',',
     &       Weights(iGrid+1),',',
     &       dTot
cGLM          if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
          if(dTot.ge.thrsrho) then
            ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
            if(l_tanhr) ratio = tanh(ratio)
            if((1.0d0-ratio).gt.thrsrho2) then
             Zeta  = sqrt(1.0d0-ratio)
* Compute alpha and beta densities when ratio < 1
             Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
             Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
            else
             Zeta  = 0.0d0
* Compute alpha and beta densities when ratio > 1
             Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
             Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
            end if
*           ^ end conditional on ratio values
* Compute gradients
          end if
*         ^ end if for little rho
          write(LuMC,'(3(F10.6,A),7(F17.10,A))')
     &          Grid(1,iGrid+1),',',
     &          Grid(2,iGrid+1),',',
     &          Grid(3,iGrid+1),',',
     &          Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &          Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &          dTot*Weights(iGrid+1),',',
     &          Weights(iGrid+1),',',
     &          dTot,',',
     &          P2_ontop(1,iGrid+1),',',
     &          ratio
         end do
*        ^ end loop over grid points
       Else !GRADIENT CALCULATION
           do iGrid=0,mGrid-1
           dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
            do dindex=1,nGrad_Eff
              dTot_d=dRho_dr(1,iGrid+1,dindex)+dRho_dr(2,iGrid+1,dindex)
              ratio = 0.0d0
              ratio_d = 0.0d0
cGLM           if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
           if(dTot.ge.thrsrho) then
             ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
             ratio_d = 4.0d0*P2_ontop_d(1,dindex,iGrid+1)/(dTot**2.0d0)
     &                - 8*P2_ontop(1,iGrid+1)*dTot_d/(dTot**3.0d0)

             if((1.0d0-ratio).gt.thrsrho2) then
              Zeta  = sqrt(1.0d0-ratio)
              Zeta_d = -ratio_d/(2*Zeta)

* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0

              dRho_dr(1,iGrid+1,dindex) =
     &        Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
              dRho_dr(2,iGrid+1,dindex) =
     &        -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0

             else
              Zeta  = 0.0d0
              Zeta_d  = 0.0d0

              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
              dRho_dr(1,iGrid+1,dindex)= dTot_d/2.0d0
              dRho_dr(2,iGrid+1,dindex)= dTot_d/2.0d0


             end if
           end if
            end do!ngrad_eff
#ifdef _DEBUG_
          !if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
!         write(*,*) Grid(1,iGrid+1),Grid(2,iGrid+1),Grid(3,iGrid+1)
          if(Grid(1,iGrid+1).eq.0d0.and.
     &    Grid(2,iGrid+1).eq.0d0) then
          write(97,'(3F12.6,12E12.4)')
     &(Grid(i,iGrid+1),i=1,3),dRho_dr(1,iGrid+1,6),dRho_dr(1,iGrid+1,9),
     &dRho_dr(2,iGrid+1,6),dRho_dr(2,iGrid+1,9),p2_ontop_d(1,6,iGrid+1),
     &p2_ontop_d(1,9,iGrid+1),
     &dTot,dTot_d,Rho(1,iGrid+1),Rho(2,iGrid+1),Zeta,p2_ontop(1,iGrid+1)
          write(98,'(3F12.6,12E12.4)')
     &Grid(1,iGrid+1),-Grid(2,iGrid+1),-Grid(3,iGrid+1),
     &dRho_dr(1,iGrid+1,6),dRho_dr(1,iGrid+1,9),
     &dRho_dr(2,iGrid+1,6),dRho_dr(2,iGrid+1,9),
     &-p2_ontop_d(1,6,iGrid+1),-p2_ontop_d(1,9,iGrid+1),
     &dTot,dTot_d,Rho(1,iGrid+1),Rho(2,iGrid+1),Zeta,p2_ontop(1,iGrid+1)
          end if
          !end if
#endif

!          write(LuMC,'(3F12.6,5E12.4)')
!     &(Grid(i,iGrid+1),i=1,3),Rho(1,iGrid+1),
!     &Rho(2,iGrid+1),dTot,P2_ontop_d(1,iGrid+1),ratio
           end do!igrid
       End if !not gradient or gradient
       End if !tlsda
#ifdef _DEBUG_
!       Close(97)
!       Close(98)
#endif


cRKCft
************************************************************************
* FTLSDA                                                               *
************************************************************************

       If(KSDFA(1:6).eq.'FTLSDA') then !GLM

         T_Rho=T_X*1.0D-4
         Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
         Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
         Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)

         Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_P2glm(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho,Do_Grad)

!      Call Do_P2GLM_AO(P2AO,np2AO,TabAO,mAO,mGrid,nAOs,
!     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
!     &          list_s,nlist_s,list_bas,Index,nIndex)


       Else !AO-based run for gradients
!        write(*,*) 'nlist_s',nlist_s
         nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
!         Call GetMem('P2_ontop_d','Allo','Real',ipP2_d,nP2_ontop_d)
!         Call FZero(Work(ipP2_d),nP2_ontop_d)
         P2_ontop_d(:,:,:) = 0
         !Determine number of AOs:
         nAOs = nMOs
      ft=.true.
        Call  Do_P2GLM_grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &          list_s,nlist_s,list_bas,Index,nIndex,
     &          P2mo,np2act,D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &      Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,nAOs,nCMO,
     &      TabSO,nsym,ft)
       End If

         If(.not.Do_Grad) then
          do iGrid=0,mGrid-1
           dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
           ratio = 0.0d0
           write(LuMT,'(3(F10.6,A),5(F17.10,A))')
     &       Grid(1,iGrid+1),',',
     &       Grid(2,iGrid+1),',',
     &       Grid(3,iGrid+1),',',
     &       Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &       Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &       dTot*Weights(iGrid+1),',',
     &       Weights(iGrid+1),',',
     &       dTot
cGLM       if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
           if(dTot.ge.thrsrho) then
             ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
             if(((1.0d0-ratio).gt.thrsrho2).and.(ratio.lt.thrsrho3))then
              Zeta  = sqrt(1.0d0-ratio)
* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
             end if
*           ^ end if
              if((ratio.ge.thrsrho3).and.(ratio.le.thrsrho4)) then
                Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &         (Bb1*(ratio-1.15d0)**4.0d0) + (Cb1*(ratio-1.15d0)**3.0d0)
* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
             end if
*         ^ end if over spline
              if (ratio.gt.thrsrho4) then
              Zeta  = 0.0d0
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
             end if
           end if
*         ^ end if for little rho
          write(LuMC,'(3(F10.6,A),7(F17.10,A))')
     &          Grid(1,iGrid+1),',',
     &          Grid(2,iGrid+1),',',
     &          Grid(3,iGrid+1),',',
     &          Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &          Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &          dTot*Weights(iGrid+1),',',
     &          Weights(iGrid+1),',',
     &          dTot,',',
     &          P2_ontop(1,iGrid+1),',',
     &          ratio
          end do

!******************************************************************
         Else !Gradient calculation
!           do i=1,nP2_ontop_d
!             write(*,*) Work(ipP2_d-1+i)
!           end do
!         Call Fzero(Work(ipP2_d),mGrid*nGrad_Eff)
           do iGrid=0,mGrid-1
!      if(.false.)then
           dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
            do dindex=1,nGrad_Eff
              dTot_d=dRho_dr(1,iGrid+1,dindex)+dRho_dr(2,iGrid+1,dindex)
              ratio = 0.0d0
              ratio_d = 0.0d0
cGLM      if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
           if(dTot.ge.thrsrho) then
             ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
!             ratio_d= 4.0d0*P2_ontop_d(dindex,iGrid+1)/(dTot**2.0d0)
!             ratio_d=4.0d0*Work(ipP2_d - 1 +
!     &               (dindex-1)*mGrid+iGrid)/(dTot**2.0d0)
              ratio_d = 4*P2_ontop_d(1,dindex,iGrid+1)/(dTot**2.0d0)
     &                - 8*P2_ontop(1,iGrid+1)*dTot_d/(dTot**3.0d0)
!             ratio_d= 0


             if(((1.0d0-ratio).gt.thrsrho2).and.(ratio.lt.thrsrho3))then
              Zeta  = sqrt(1.0d0-ratio)
              Zeta_d = -ratio_d/(2*Zeta)

* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0

              dRho_dr(1,iGrid+1,dindex) =
     &        Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
              dRho_dr(2,iGrid+1,dindex) =
     &        -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0

             end if

!      end if
              if((ratio.ge.thrsrho3).and.(ratio.le.thrsrho4)) then
                Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &         (Bb1*(ratio-1.15d0)**4.0d0) + (Cb1*(ratio-1.15d0)**3.0d0)
                Zeta_d = 5*Ab1*(ratio-1.15d0)**4 * ratio_d
     &        + 4*Bb1*(ratio-1.15d0)**3 * ratio_d
     &        + 3*Cb1*(ratio-1.15d0)**2 * ratio_d

* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
              dRho_dr(1,iGrid+1,dindex) =
     &        Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
              dRho_dr(2,iGrid+1,dindex) =
     &        -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0

             end if
*         ^ end if over spline
              if (ratio.gt.thrsrho4) then
              Zeta  = 0.0d0
              Zeta_d  = 0.0d0

              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
              dRho_dr(1,iGrid+1,dindex)= dTot_d/2.0d0
              dRho_dr(2,iGrid+1,dindex)= dTot_d/2.0d0
             end if
           end if
*          ^ end if for little rho
          write(LuMC,'(3(F10.6,A),7(F17.10,A))')
     &          Grid(1,iGrid+1),',',
     &          Grid(2,iGrid+1),',',
     &          Grid(3,iGrid+1),',',
     &          Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &          Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &          dTot*Weights(iGrid+1),',',
     &          Weights(iGrid+1),',',
     &          dTot,',',
     &          P2_ontop(1,iGrid+1),',',
     &          ratio
          end do!loop over effective coordinates
          end do!loop over gridpts
!         Call GetMem('P2_ontop_d','Free','Real',ipP2_d,mGrid*nGrad_Eff)
         End If!Gradient calculation

*        ^ end loop over grid points
c         write(6,*)'spinDens aft on-top density'
c         do iGrid=1,mGrid
c           write(6,'(2E12.6)') Rho(1,iGrid), Rho(2,iGrid)
c         end do
       End if


*      ^ end if for GLM stuff
C        If (Do_Hess)
C    &      Call d2Rho_dR2_LDA(Dens,nDens,nD,dRho_dR,d2Rho_dr2,
C    &                         ndRho_dr,mGrid,list_s,nlist_s,
C    &                         TabAO,ipTabAO,mAO,nTabAO,nSym,
C    &                         nGrad_Eff,list_g,Maps2p,
C    &                         nShell,Grid_Type,Fixed_Grid,
C    &                         Work(ip_Fact),ndc,Work(ipTmp),T_X,
C    &                         list_bas,Index,nIndex)
*
      Else If (Functional_type.eq.GGA_type) Then
*
         Call Rho_GGA(Dens,nDens,nD,Rho,nRho,mGrid,
     &                list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                Work(ip_Fact),ndc,Work(ipTabAOMax),
     &                list_bas,Index,nIndex)
*
         If (Do_Grad)
     &      Call dRho_dR_GGA(Dens,nDens,nD,dRho_dR,ndRho_dr,
     &                       mGrid,list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                       nGrad_Eff,list_g,Maps2p,
     &                       nShell,Grid_Type,Fixed_Grid,
     &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
     &                       list_bas,Index,nIndex)

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
*  iSwitch = 0  total density
*  iSwitch = 1  alpha density
*  iSwitch = 2  beta density
        T_Rho=T_X*1.0D-4
        Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
        Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
        Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)
c        write(6,*) 'in do_batch.f for GLM=GGA options'
C         If (Do_MO) Then
C            Call Do_RhoIA(nRho,mGrid,Rho,Work(ipRhoI),
C     &                    Work(ipRhoA),mRho,TabMO,mAO,
C     &                    nMOs,D1mo,nd1mo)
C         Else
C            Call WarningMessage(2,' We need MO for RhoIA')
C            Call Abend()
C         End If

         Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_P2glm(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho,Do_Grad)
       Else !AO-based run for gradients
         nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
         P2_ontop_d(:,:,:) = 0
         !Determine number of AOs:
         nAOs = nMOs
      ft=.false.
        Call  Do_P2GLM_grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &          list_s,nlist_s,list_bas,Index,nIndex,
     &          P2mo,np2act,D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &      Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,nAOs,nCMO,
     &      TabSO,nsym,ft)
       End If


c         write(6,*) 'X -- Y -- Z -- P2 -- gradX -- gradY -- gradZ'
c         do iGrid=1,mGrid
c             write(6,'(3F8.3,4E12.4)') (Grid(i,iGrid),i=1,3),
c     &       P2_ontop(1,iGrid),
c     &       P2_ontop(2,iGrid),P2_ontop(3,iGrid),P2_ontop(4,iGrid)
c         end do

C         write(6,*)'X Y Z d_a d_b grXa grYa grZa grXb grYb grZb bf P2'
C         do iGrid=1,mGrid
C           write(6,'(3F8.3,8E12.4)')
C     &           (Grid(i,iGrid),i=1,3),(Rho(j,iGrid),j=1,8)
C         end do
c         write(6,*) 'thrsrho', thrsrho
c         write(6,*) 'thrsrho2', thrsrho2

       If(.not.Do_Grad) then
         do iGrid=0,mGrid-1
          dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
          grad_x = Rho(3,iGrid+1)+Rho(6,iGrid+1)
          grad_y = Rho(4,iGrid+1)+Rho(7,iGrid+1)
          grad_z = Rho(5,iGrid+1)+Rho(8,iGrid+1)
          ratio = 0.0d0
          write(LuMT,'(3(F10.6,A),5(F17.10,A))')
     &       Grid(1,iGrid+1),',',
     &       Grid(2,iGrid+1),',',
     &       Grid(3,iGrid+1),',',
     &       Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &       Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &       dTot*Weights(iGrid+1),',',
     &       Weights(iGrid+1),',',
     &       dTot
cGLM       if(dTot.ge.thrsrho.and.abs(P2_ontop(1,iGrid+1)).ge.
cGLM     &                               0.25*thrsrho**3.0d0)then
          if(dTot.ge.thrsrho) then
            ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
            if(l_tanhr) ratio = tanh(ratio)
            if((1.0d0-ratio).gt.thrsrho2) then
             Zeta  = sqrt(1.0d0-ratio)
             pi_p_x = (1.0d0 - zeta**2.0d0)*dTot*grad_x/2.0d0
             pi_p_y = (1.0d0 - zeta**2.0d0)*dTot*grad_y/2.0d0
             pi_p_z = (1.0d0 - zeta**2.0d0)*dTot*grad_z/2.0d0
* Compute alpha and beta densities for ratio < 1
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
* Compute  gradients
*              write(6,*) 'zeta, ratio, dTot',Zeta, ratio, dTot
*              write(6,*) 'Rho_grad = ', grad_x, grad_y, grad_z
*              write(6,*) 'P2_grad  = ',
*     &     P2_ontop(2,iGrid+1),P2_ontop(3,iGrid+1),P2_ontop(4,iGrid+1)
c              write(6,*) 'P2_p-pi_p =', P2_ontop(2,iGrid+1)-pi_p_x,
c     &        P2_ontop(3,iGrid+1)-pi_p_y, P2_ontop(4,iGrid+1)-pi_p_z
*              write(6,*) 'Old Rho = ',
*     &        Rho(3,iGrid+1),Rho(4,iGrid+1),Rho(5,iGrid+1),
*     &        Rho(6,iGrid+1),Rho(7,iGrid+1),Rho(8,iGrid+1)
        if(l_tanhr) then
          Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
     &    + ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(2,iGrid+1)-(grad_x*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
c
          Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
     &    + ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(3,iGrid+1)-(grad_y*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
c
          Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
     &    + ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(4,iGrid+1)-(grad_z*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
c
          Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
     &    - ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(2,iGrid+1)-(grad_x*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
c
          Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
     &    - ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(3,iGrid+1)-(grad_y*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
c
          Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
     &    - ((ratio**2.0d0 - 1.0d0)
     &    *(P2_ontop(4,iGrid+1)-(grad_z*4.0d0*P2_ontop(1,iGrid+1)/dTot))
     &    /(dTot*zeta))
        else
              Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
*     &                         + ratio*grad_x/(2.0d0*Zeta)
*     &                         - P2_ontop(2,iGrid+1)/(dTot*Zeta)
              Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
*     &                         + ratio*grad_y/(2.0d0*Zeta)
*     &                         - P2_ontop(3,iGrid+1)/(dTot*Zeta)
              Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
*     &                         + ratio*grad_z/(2.0d0*Zeta)
*     &                         - P2_ontop(4,iGrid+1)/(dTot*Zeta)
              Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
*     &                         - ratio*grad_x/(2.0d0*Zeta)
*     &                         + P2_ontop(2,iGrid+1)/(dTot*Zeta)
              Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
*     &                         - ratio*grad_y/(2.0d0*Zeta)
*     &                         + P2_ontop(3,iGrid+1)/(dTot*Zeta)
              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
*     &                         - ratio*grad_z/(2.0d0*Zeta)
*     &                         + P2_ontop(4,iGrid+1)/(dTot*Zeta)
        end if ! End tanh option
*              write(6,*) 'New Rho = ',
*     &        Rho(3,iGrid+1),Rho(4,iGrid+1),Rho(5,iGrid+1),
*     &        Rho(6,iGrid+1),Rho(7,iGrid+1),Rho(8,iGrid+1)
           else
             Zeta  = 0.0d0
             pi_p_x = (1.0d0 - zeta**2.0d0)*dTot*grad_x/2.0d0
             pi_p_y = (1.0d0 - zeta**2.0d0)*dTot*grad_y/2.0d0
             pi_p_z = (1.0d0 - zeta**2.0d0)*dTot*grad_z/2.0d0
* Compute alpha and beta densities if ratio > 1
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
* Compute  gradients
              Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
              Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
              Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
              Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
              Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
           end if
*          ^ end conditional regarding the ratio
          end if
*         ^ end if over little density
          write(LuMC,'(3(F10.6,A),7(F17.10,A))')
     &          Grid(1,iGrid+1),',',
     &          Grid(2,iGrid+1),',',
     &          Grid(3,iGrid+1),',',
     &          Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &          Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &          dTot*Weights(iGrid+1),',',
     &          Weights(iGrid+1),',',
     &          dTot,',',
     &          P2_ontop(1,iGrid+1),',',
     &          ratio
        end do!end loop over grid points

       Else !GRADIENT CALCULATION
           do iGrid=0,mGrid-1
           dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
           grad_x = Rho(3,iGrid+1)+Rho(6,iGrid+1)
           grad_y = Rho(4,iGrid+1)+Rho(7,iGrid+1)
           grad_z = Rho(5,iGrid+1)+Rho(8,iGrid+1)
            do dindex=1,nGrad_Eff
              dTot_d=dRho_dr(1,iGrid+1,dindex)+dRho_dr(2,iGrid+1,dindex)
              dTot_dx=dRho_dr(3,iGrid+1,dindex)
     &               +dRho_dr(6,iGrid+1,dindex)
              dTot_dy=dRho_dr(4,iGrid+1,dindex)
     &               +dRho_dr(7,iGrid+1,dindex)
              dTot_dz=dRho_dr(5,iGrid+1,dindex)
     &               +dRho_dr(8,iGrid+1,dindex)
              ratio = 0.0d0
              ratio_d = 0.0d0
cGLM          if(P2_ontop(1,iGrid+1).lt.0.0d0) P2_ontop(1,iGrid+1)
cGLM     &                                           = 0.0d0
cGLM       if(dTot.ge.thrsrho.and.abs(P2_ontop(1,iGrid+1))
cGLM     &                           .ge.0.25*thrsrho**3.0d0) then
          if(dTot.ge.thrsrho) then
              ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
              ratio_d = 4.0d0*P2_ontop_d(1,dindex,iGrid+1)/(dTot**2.0d0)
     &                - 8*P2_ontop(1,iGrid+1)*dTot_d/(dTot**3.0d0)

              if((1.0d0-ratio).gt.thrsrho2) then
               Zeta  = sqrt(1.0d0-ratio)
               Zeta_d = -ratio_d/(2*Zeta)

* Compute alpha and beta densities
               Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
               Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
               Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0

               dRho_dr(1,iGrid+1,dindex) =
     &         Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
               dRho_dr(2,iGrid+1,dindex) =
     &         -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0
               dRho_dr(3,iGrid+1,dindex) =
     &         Zeta_d*grad_x/2.0d0+(1.0d0+Zeta)*dTot_dx/2.0d0
               dRho_dr(4,iGrid+1,dindex) =
     &         Zeta_d*grad_y/2.0d0+(1.0d0+Zeta)*dTot_dy/2.0d0
               dRho_dr(5,iGrid+1,dindex) =
     &         Zeta_d*grad_z/2.0d0+(1.0d0+Zeta)*dTot_dz/2.0d0
               dRho_dr(6,iGrid+1,dindex) =
     &         -Zeta_d*grad_x/2.0d0+(1.0d0-Zeta)*dTot_dx/2.0d0
               dRho_dr(7,iGrid+1,dindex) =
     &         -Zeta_d*grad_y/2.0d0+(1.0d0-Zeta)*dTot_dy/2.0d0
               dRho_dr(8,iGrid+1,dindex) =
     &         -Zeta_d*grad_z/2.0d0+(1.0d0-Zeta)*dTot_dz/2.0d0
              else
               Zeta  = 0.0d0
               Zeta_d  = 0.0d0

               Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
               Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
               Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0

               dRho_dr(1,iGrid+1,dindex) =
     &         Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
               dRho_dr(2,iGrid+1,dindex) =
     &         -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0
               dRho_dr(3,iGrid+1,dindex) =
     &         Zeta_d*grad_x/2.0d0+(1.0d0+Zeta)*dTot_dx/2.0d0
               dRho_dr(4,iGrid+1,dindex) =
     &         Zeta_d*grad_y/2.0d0+(1.0d0+Zeta)*dTot_dy/2.0d0
               dRho_dr(5,iGrid+1,dindex) =
     &         Zeta_d*grad_z/2.0d0+(1.0d0+Zeta)*dTot_dz/2.0d0
               dRho_dr(6,iGrid+1,dindex) =
     &         -Zeta_d*grad_x/2.0d0+(1.0d0-Zeta)*dTot_dx/2.0d0
               dRho_dr(7,iGrid+1,dindex) =
     &         -Zeta_d*grad_y/2.0d0+(1.0d0-Zeta)*dTot_dy/2.0d0
               dRho_dr(8,iGrid+1,dindex) =
     &         -Zeta_d*grad_z/2.0d0+(1.0d0-Zeta)*dTot_dz/2.0d0
              end if!Threshrho_2
            end if!little density (threshrho)
          end do!ngrad_eff
        end do!gridpt
       end if!Gradient calculation
*       ^ end loop over grid points
c         write(6,*)'X Y Z spinDens and grad aft on-top density'
*         do iGrid=1,mGrid
*           write(6,'(3F8.3,3E12.4,4A)')
*     &           (Grid(i,iGrid),i=1,3),
*     &           Rho(1,iGrid), Rho(2,iGrid),(Rho(1,iGrid)+Rho(2,iGrid)),
*     &           '  ciao'
*         end do
*         do iGrid=1,mGrid
*            write(6,'(3F12.6,8E14.6,8A)')
*     &           (Grid(i,iGrid),i=1,3),
*     &           Rho(1,iGrid),Rho(2,iGrid),
*     &           Rho(3,iGrid),Rho(4,iGrid),Rho(5,iGrid),
*     &           Rho(6,iGrid),Rho(7,iGrid),Rho(8,iGrid),
*     &           ' after'
*         end do
      end if
*     ^ End if over GLM stuff
************************************************************************
* FTBLYP,FTPBE,FTREVPBE                                                *
************************************************************************
      if(KSDFA(1:6).eq.'FTBLYP'.or. !GLM
     &   KSDFA(1:5).eq.'FTPBE'.or.
     &   KSDFA(1:6).eq.'FTOPBE'.or.
     &   KSDFA(1:8).eq.'FTREVPBE') then
*  *
        T_Rho=T_X*1.0D-4
        Dens_t1=Dens_t1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,0)
        Dens_a1=Dens_a1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,1)
        Dens_b1=Dens_b1+Comp_d(Weights,mGrid,Rho,nRho,nD,T_Rho,2)

          Call Fzero(P2_ontop,nP2_ontop*mGrid)
       If (.not.Do_Grad) then !regular MO-based run
         Call Do_P2glm(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho,Do_Grad)
       Else !AO-based run for gradients
         nP2_ontop_d = nP2_ontop*mGrid*nGrad_Eff
         P2_ontop_d(:,:,:) = 0
         !Determine number of AOs:
         nAOs = nMOs

      ft=.true.
        Call  Do_P2GLM_grad(TabAO,nTabAO,mAO,mGrid,ipTabAO,
     &          P2_ontop,nP2_ontop,Do_Grad,nGrad_Eff,
     &          list_s,nlist_s,list_bas,Index,nIndex,
     &          P2mo,np2act,D1mo,nd1mo,TabMO,list_g,P2_ontop_d,
     &      Work(ipRhoI),Work(ipRhoA),mRho,nMOs,CMOs,nAOs,nCMO,
     &      TabSO,nsym,ft)
       End If
*
       If(.not.Do_Grad) then
         do iGrid=0,mGrid-1
          dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
          grad_x = Rho(3,iGrid+1)+Rho(6,iGrid+1)
          grad_y = Rho(4,iGrid+1)+Rho(7,iGrid+1)
          grad_z = Rho(5,iGrid+1)+Rho(8,iGrid+1)
          ratio = 0.0d0
          write(LuMT,'(3(F10.6,A),5(F17.10,A))')
     &       Grid(1,iGrid+1),',',
     &       Grid(2,iGrid+1),',',
     &       Grid(3,iGrid+1),',',
     &       Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &       Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &       dTot*Weights(iGrid+1),',',
     &       Weights(iGrid+1),',',
     &       dTot
cGLM      if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
          if(dTot.ge.thrsrho) then
            ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
            if(((1.0d0-ratio).gt.thrsrho2).and.(ratio.lt.thrsrho3)) then
             Zeta  = sqrt(1.0d0-ratio)
             pi_p_x = (1.0d0 - zeta**2.0d0)*dTot*grad_x/2.0d0
             pi_p_y = (1.0d0 - zeta**2.0d0)*dTot*grad_y/2.0d0
             pi_p_z = (1.0d0 - zeta**2.0d0)*dTot*grad_z/2.0d0
* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
* compute gradients
              Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
     &                         + ratio*grad_x/(2.0d0*Zeta)
     &                         - P2_ontop(2,iGrid+1)/(dTot*Zeta)
              Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
     &                         + ratio*grad_y/(2.0d0*Zeta)
     &                         - P2_ontop(3,iGrid+1)/(dTot*Zeta)
              Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
     &                         + ratio*grad_z/(2.0d0*Zeta)
     &                         - P2_ontop(4,iGrid+1)/(dTot*Zeta)
              Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
     &                         - ratio*grad_x/(2.0d0*Zeta)
     &                         + P2_ontop(2,iGrid+1)/(dTot*Zeta)
              Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
     &                         - ratio*grad_y/(2.0d0*Zeta)
     &                         + P2_ontop(3,iGrid+1)/(dTot*Zeta)
              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
     &                         - ratio*grad_z/(2.0d0*Zeta)
     &                         + P2_ontop(4,iGrid+1)/(dTot*Zeta)
*          ^ end if over unwanted ratios
            else if((ratio.ge.thrsrho3).and.(ratio.le.thrsrho4)) then

               Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &        (Bb1*(ratio-1.15d0)**4.0d0) + (Cb1*(ratio-1.15d0)**3.0d0)

             pi_p_x = (1.0d0 - zeta**2.0d0)*dTot*grad_x/2.0d0
             pi_p_y = (1.0d0 - zeta**2.0d0)*dTot*grad_y/2.0d0
             pi_p_z = (1.0d0 - zeta**2.0d0)*dTot*grad_z/2.0d0

* Compute alpha and beta densities

              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0

              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0

* Compute  gradients
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_x))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_x))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_x))

               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_y))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_y))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_y))

               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_z))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_z))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_z))

               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_x))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_x))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_x))

               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_y))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_y))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_y))

              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_z))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_z))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_z))

*         ^ end if over spline
           else if (ratio.gt.thrsrho4) then
             Zeta  = 0.0d0
* Compute alpha and beta densities
              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
* Compute  gradients
              Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
              Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
              Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
              Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
              Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
           end if
*          ^ end if over zeta=0
          end if
*         ^ end if over little density
          write(LuMC,'(3(F10.6,A),7(F17.10,A))')
     &          Grid(1,iGrid+1),',',
     &          Grid(2,iGrid+1),',',
     &          Grid(3,iGrid+1),',',
     &          Rho(1,iGrid+1)*Weights(iGrid+1),',',
     &          Rho(2,iGrid+1)*Weights(iGrid+1),',',
     &          dTot*Weights(iGrid+1),',',
     &          Weights(iGrid+1),',',
     &          dTot,',',
     &          P2_ontop(1,iGrid+1),',',
     &          ratio
        end do
*       ^ end loop over grid points


        Else !GRADIENT CALCULATION
           do iGrid=0,mGrid-1
           dTot=Rho(1,iGrid+1)+Rho(2,iGrid+1)
           grad_x = Rho(3,iGrid+1)+Rho(6,iGrid+1)
           grad_y = Rho(4,iGrid+1)+Rho(7,iGrid+1)
           grad_z = Rho(5,iGrid+1)+Rho(8,iGrid+1)
            do dindex=1,nGrad_Eff
              dTot_d=dRho_dr(1,iGrid+1,dindex)+dRho_dr(2,iGrid+1,dindex)
              dTot_dx=dRho_dr(3,iGrid+1,dindex)
     &               +dRho_dr(6,iGrid+1,dindex)
              dTot_dy=dRho_dr(4,iGrid+1,dindex)
     &               +dRho_dr(7,iGrid+1,dindex)
              dTot_dz=dRho_dr(5,iGrid+1,dindex)
     &               +dRho_dr(8,iGrid+1,dindex)
              ratio = 0.0d0
              ratio_d = 0.0d0
cGLM        if(dTot.ge.thrsrho.and.P2_ontop(1,iGrid+1).ge.thrsrho) then
            if(dTot.ge.thrsrho) then
              ratio = 4.0d0*P2_ontop(1,iGrid+1)/(dTot**2.0d0)
              ratio_d = 4.0d0*P2_ontop_d(1,dindex,iGrid+1)/(dTot**2.0d0)
     &                - 8.0d0*P2_ontop(1,iGrid+1)*dTot_d/(dTot**3.0d0)
              ratio_dx = 4.0d0*P2_ontop(2,iGrid+1)/(dTot**2.0d0)
     &                   -2.0d0*ratio/dTot*grad_x
              ratio_dy = 4.0d0*P2_ontop(3,iGrid+1)/(dTot**2.0d0)
     &                   -2.0d0*ratio/dTot*grad_y
              ratio_dz = 4.0d0*P2_ontop(4,iGrid+1)/(dTot**2.0d0)
     &                   -2.0d0*ratio/dTot*grad_z
             d_ratiox = 4.0d0*P2_ontop_d(2,dindex,iGrid+1)/(dTot**2.0d0)
     &                 -8.0d0*P2_ontop(2,iGrid+1)*dTot_d/(dTot**3.0d0)
     &                 -2.0d0*ratio/dTot*dTot_dx
     &             -8.0d0*grad_x/(dtot**3)*P2_ontop_d(1,dindex,iGrid+1)
     &             +6d0*ratio*grad_x/(dtot**2) * dtot_d
!     &                 -2.0d0*grad_x/dTot*ratio_d
!     &                 +2.0d0*grad_x*ratio/(dTot**2.0d0)*dTot_d
             d_ratioy = 4.0d0*P2_ontop_d(3,dindex,iGrid+1)/(dTot**2.0d0)
     &                 -8.0d0*P2_ontop(3,iGrid+1)*dTot_d/(dTot**3.0d0)
     &                 -2.0d0*ratio/dTot*dTot_dy
     &             -8.0d0*grad_y/(dtot**3)*P2_ontop_d(1,dindex,iGrid+1)
     &             +6d0*ratio*grad_y/(dtot**2) * dtot_d
!     &                 -2.0d0*grad_y/dTot*ratio_d
!     &                 +2.0d0*grad_y*ratio/(dTot**2.0d0)*dTot_d
            d_ratioz = 4.0d0*P2_ontop_d(4,dindex,iGrid+1)/(dTot**2.0d0)
     &                 -8.0d0*P2_ontop(4,iGrid+1)*dTot_d/(dTot**3.0d0)
     &                 -2.0d0*ratio/dTot*dTot_dz
     &             -8.0d0*grad_z/(dtot**3)*P2_ontop_d(1,dindex,iGrid+1)
     &             +6d0*ratio*grad_z/(dtot**2) * dtot_d
!     &                 -2.0d0*grad_z/dTot*ratio_d
!     &                 +2.0d0*grad_z*ratio/(dTot**2.0d0)*dTot_d
           if(((1.0d0-ratio).gt.thrsrho2).and.(ratio.lt.thrsrho3)) then
!           if(.false.) then
               Zeta  = sqrt(1.0d0-ratio)
               Zeta_d = -1.0d0*ratio_d/(2*Zeta)
               Zeta_dx = -1.0d0*ratio_dx/(2*Zeta)
               Zeta_dy = -1.0d0*ratio_dy/(2*Zeta)
               Zeta_dz = -1.0d0*ratio_dz/(2*Zeta)
               !dZeta/dLambda
               d_Zetax = -1.0d0*d_ratiox/(2.0d0*Zeta)
     &                  +ratio_dx/(2.0d0*Zeta**2.0d0)*Zeta_d
               d_Zetay = -1.0d0*d_ratioy/(2.0d0*Zeta)
     &                  +ratio_dy/(2.0d0*Zeta**2.0d0)*Zeta_d
               d_Zetaz = -1.0d0*d_ratioz/(2.0d0*Zeta)
     &                  +ratio_dz/(2.0d0*Zeta**2.0d0)*Zeta_d
! Compute alpha and beta densities
               Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
               Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
     &                          +0.5d0*dTot*Zeta_dx
               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
     &                          +0.5d0*dTot*Zeta_dy
               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
     &                          +0.5d0*dTot*Zeta_dz
               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
     &                          -0.5d0*dTot*Zeta_dx
               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
     &                          -0.5d0*dTot*Zeta_dy
               Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
     &                          -0.5d0*dTot*Zeta_dz

               dRho_dr(1,iGrid+1,dindex) =
     &         Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
               dRho_dr(2,iGrid+1,dindex) =
     &         -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0

               dRho_dr(3,iGrid+1,dindex) =
     &         Zeta_d*grad_x/2.0d0+(1.0d0+Zeta)*dTot_dx/2.0d0
     &         +dTot/2.0d0*d_Zetax+Zeta_dx/2.0d0*dTot_d

               dRho_dr(4,iGrid+1,dindex) =
     &         Zeta_d*grad_y/2.0d0+(1.0d0+Zeta)*dTot_dy/2.0d0
     &         +dTot/2.0d0*d_Zetay+Zeta_dy/2.0d0*dTot_d

               dRho_dr(5,iGrid+1,dindex) =
     &         Zeta_d*grad_z/2.0d0+(1.0d0+Zeta)*dTot_dz/2.0d0
     &         +dTot/2.0d0*d_Zetaz+Zeta_dz/2.0d0*dTot_d

               dRho_dr(6,iGrid+1,dindex) =
     &         -Zeta_d*grad_x/2.0d0+(1.0d0-Zeta)*dTot_dx/2.0d0
     &         -dTot/2.0d0*d_Zetax-Zeta_dx/2.0d0*dTot_d

               dRho_dr(7,iGrid+1,dindex) =
     &         -Zeta_d*grad_y/2.0d0+(1.0d0-Zeta)*dTot_dy/2.0d0
     &         -dTot/2.0d0*d_Zetay-Zeta_dy/2.0d0*dTot_d

               dRho_dr(8,iGrid+1,dindex) =
     &         -Zeta_d*grad_z/2.0d0+(1.0d0-Zeta)*dTot_dz/2.0d0
     &         -dTot/2.0d0*d_Zetaz-Zeta_dz/2.0d0*dTot_d

*
            else if((ratio.ge.thrsrho3).and.(ratio.le.thrsrho4)) then
            !else if(.false.) then

               Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &        (Bb1*(ratio-1.15d0)**4.0d0) + (Cb1*(ratio-1.15d0)**3.0d0)
               Deriv =(5.0d0*Ab1*(ratio-1.15d0)**4.0d0)
     &               +(4.0d0*Bb1*(ratio-1.15d0)**3.0d0)
     &               +(3.0d0*Cb1*(ratio-1.15d0)**2.0d0)
               d_Deriv = (20.0d0*Ab1*(ratio-1.15d0)**3.0d0)
     &               +(12.0d0*Bb1*(ratio-1.15d0)**2.0d0)
     &               +(6.0d0*Cb1*(ratio-1.15d0)**1.0d0)

               Zeta_d = ratio_d*Deriv
               Zeta_dx = ratio_dx*Deriv
               Zeta_dy = ratio_dy*Deriv
               Zeta_dz = ratio_dz*Deriv
               !dZeta/dLambda
               d_Zetax = d_ratiox*Deriv+ratio_dx*d_Deriv*ratio_d
               d_Zetay = d_ratioy*Deriv+ratio_dy*d_Deriv*ratio_d
               d_Zetaz = d_ratioz*Deriv+ratio_dz*d_Deriv*ratio_d


* Compute alpha and beta densities

              Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0

              Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0

* Compute  gradients
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_x))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_x))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_x))
!
               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_y))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_y))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_y))

               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((10.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (5.0d0*ratio*grad_z))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((8.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (4.0d0*ratio*grad_z))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((6.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            - (3.0d0*ratio*grad_z))

               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_x))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_x))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(2,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_x))

               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_y))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_y))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(3,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_y))

              Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
     &           + (Ab1*(ratio-1.15d0)**4.0d0)
     &            * ((-10.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (5.0d0*ratio*grad_z))
     &           + (Bb1*(ratio-1.15d0)**3.0d0)
     &            * ((-8.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (4.0d0*ratio*grad_z))
     &           + (Cb1*(ratio-1.15d0)**2.0d0)
     &            * ((-6.0d0*P2_ontop(4,iGrid+1)/dTot)
     &            + (3.0d0*ratio*grad_z))

               dRho_dr(1,iGrid+1,dindex) =
     &         Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0

               dRho_dr(2,iGrid+1,dindex) =
     &         -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0

               dRho_dr(3,iGrid+1,dindex) =
     &         Zeta_d*grad_x/2.0d0+(1.0d0+Zeta)*dTot_dx/2.0d0
     &         +dTot/2.0D0*d_Zetax+Zeta_dx/2.0d0*dTot_d

               dRho_dr(4,iGrid+1,dindex) =
     &         Zeta_d*grad_y/2.0d0+(1.0d0+Zeta)*dTot_dy/2.0d0
     &         +dTot/2.0d0*d_Zetay+Zeta_dy/2.0d0*dTot_d

               dRho_dr(5,iGrid+1,dindex) =
     &         Zeta_d*grad_z/2.0d0+(1.0d0+Zeta)*dTot_dz/2.0d0
     &         +dTot/2.0d0*d_Zetaz+Zeta_dz/2.0d0*dTot_d

               dRho_dr(6,iGrid+1,dindex) =
     &         -1.0d0*Zeta_d*grad_x/2.0d0+(1.0d0-Zeta)*dTot_dx/2.0d0
     &         -dTot/2.0D0*d_Zetax-Zeta_dx/2.0d0*dTot_d

               dRho_dr(7,iGrid+1,dindex) =
     &         -1.0d0*Zeta_d*grad_y/2.0d0+(1.0d0-Zeta)*dTot_dy/2.0d0
     &         -dTot/2.0D0*d_Zetay-Zeta_dy/2.0d0*dTot_d

               dRho_dr(8,iGrid+1,dindex) =
     &         -1.0d0*Zeta_d*grad_z/2.0d0+(1.0d0-Zeta)*dTot_dz/2.0d0
     &         -dTot/2.0D0*d_Zetaz-Zeta_dz/2.0d0*dTot_d

             else
               Zeta  = 0.0d0
               Zeta_d  = 0.0d0

!               Rho(1,iGrid+1)=0
!               Rho(2,iGrid+1)=0
!               Rho(3,iGrid+1)= 0
!               Rho(4,iGrid+1)= 0
!               Rho(5,iGrid+1)= 0
!               Rho(6,iGrid+1)= 0
!               Rho(7,iGrid+1)= 0
!               Rho(8,iGrid+1)= 0
!               dRho_dr(1,iGrid+1,dindex) =0
!               dRho_dr(2,iGrid+1,dindex) =0
!               dRho_dr(3,iGrid+1,dindex) =0
!               dRho_dr(4,iGrid+1,dindex) =0
!               dRho_dr(5,iGrid+1,dindex) =0
!               dRho_dr(6,iGrid+1,dindex) =0
!               dRho_dr(7,iGrid+1,dindex) =0
!               dRho_dr(8,iGrid+1,dindex) = 0

               Rho(1,iGrid+1)=(1.0d0+Zeta)*dTot/2.0d0
               Rho(2,iGrid+1)=(1.0d0-Zeta)*dTot/2.0d0
               Rho(3,iGrid+1)= (1.0d0+Zeta)*grad_x/2.0d0
               Rho(4,iGrid+1)= (1.0d0+Zeta)*grad_y/2.0d0
               Rho(5,iGrid+1)= (1.0d0+Zeta)*grad_z/2.0d0
               Rho(6,iGrid+1)= (1.0d0-Zeta)*grad_x/2.0d0
               Rho(7,iGrid+1)= (1.0d0-Zeta)*grad_y/2.0d0
               Rho(8,iGrid+1)= (1.0d0-Zeta)*grad_z/2.0d0
!
               dRho_dr(1,iGrid+1,dindex) =
     &         Zeta_d*dTot/2.0d0+(1.0d0+Zeta)*dTot_d/2.0d0
               dRho_dr(2,iGrid+1,dindex) =
     &         -Zeta_d*dTot/2.0d0+(1.0d0-Zeta)*dTot_d/2.0d0
               dRho_dr(3,iGrid+1,dindex) =
     &         Zeta_d*grad_x/2.0d0+(1.0d0+Zeta)*dTot_dx/2.0d0
               dRho_dr(4,iGrid+1,dindex) =
     &         Zeta_d*grad_y/2.0d0+(1.0d0+Zeta)*dTot_dy/2.0d0
               dRho_dr(5,iGrid+1,dindex) =
     &         Zeta_d*grad_z/2.0d0+(1.0d0+Zeta)*dTot_dz/2.0d0
               dRho_dr(6,iGrid+1,dindex) =
     &         -Zeta_d*grad_x/2.0d0+(1.0d0-Zeta)*dTot_dx/2.0d0
               dRho_dr(7,iGrid+1,dindex) =
     &         -Zeta_d*grad_y/2.0d0+(1.0d0-Zeta)*dTot_dy/2.0d0
               dRho_dr(8,iGrid+1,dindex) =
     &         -Zeta_d*grad_z/2.0d0+(1.0d0-Zeta)*dTot_dz/2.0d0
              end if!Threshrho_2
            end if! little density (threshrho)
          end do!ngrad_eff
        end do!gridpt
       end if!
       end if
C        If (Do_Hess)
C    &      Call dRho_dR_GGA(Dens,nDens,nD,dRho_dR,d2Rho_dR2,
C    &                       ndRho_dr,mGrid,list_s,nlist_s,
C    &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
C    &                       nGrad_Eff,list_g,Maps2p,
C    &                       nShell,Grid_Type,Fixed_Grid,
C    &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
C    &                       list_bas,Index,nIndex)
*
      Else If (Functional_type.eq.CASDFT_type) Then
*
         Call Rho_CAS(Dens,nDens,nD,Rho,nRho,mGrid,
     &                list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                Work(ip_Fact),ndc,Work(ipTabAOMax),
     &                list_bas,Index,nIndex)
*
C        If (Do_Grad)
C    &      Call dRho_dR_CAS(Dens,nDens,nD,dRho_dR,ndRho_dr,
C    &                       mGrid,list_s,nlist_s,
C    &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
C    &                       nGrad_Eff,list_g,Maps2p,
C    &                       nShell,Grid_Type,Fixed_Grid,
C    &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
C    &                       list_bas,Index,nIndex)
C        If (Do_Hess)
C    &      Call dRho_dR_CAS(Dens,nDens,nD,dRho_dR,d2Rho_dR2,
C    &                       ndRho_dr,mGrid,list_s,nlist_s,
C    &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
C    &                       nGrad_Eff,list_g,Maps2p,
C    &                       nShell,Grid_Type,Fixed_Grid,
C    &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
C    &                       list_bas,Index,nIndex)
*------- Compute P2_OnTop at the grid
         Call Do_P2new(P2mo,np2act,D1mo,nd1mo,TabMO,mAO,mGrid,
     &                 nMOs,P2_ontop,nP2_ontop,Work(ipRhoI),
     &                 Work(ipRhoA),mRho)
*
      Else If (Functional_type.eq.meta_GGA_type1) Then
*
         Call Rho_meta_GGA1(nD,Rho,nRho,mGrid,
     &                     list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                     Work(ip_Fact),ndc,Work(ipTabAOMax),
     &                     list_bas,Index,nIndex)
*
         If (Do_Grad)
     &      Call dRho_dR_meta_GGA1
     &                      (nD,dRho_dR,ndRho_dr,
     &                       mGrid,list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                       nGrad_Eff,list_g,Maps2p,nShell,
     &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
     &                       list_bas,Index,nIndex)
*
      Else If (Functional_type.eq.meta_GGA_type2) Then
*
         Call Rho_meta_GGA2(nD,Rho,nRho,mGrid,
     &                     list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                     Work(ip_Fact),ndc,Work(ipTabAOMax),
     &                     list_bas,Index,nIndex)
*
         If (Do_Grad)
     &      Call dRho_dR_meta_GGA2
     &                      (nD,dRho_dR,ndRho_dr,
     &                       mGrid,list_s,nlist_s,
     &                       TabAO,ipTabAO,mAO,nTabAO,nSym,
     &                       nGrad_Eff,list_g,Maps2p,nShell,
     &                       Work(ip_Fact),ndc,Work(ipTmp),T_X,
     &                       list_bas,Index,nIndex)
      End If
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
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,Rho,nRho,nD,T_Rho)
      Else If (Functional_type.eq.GGA_type) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,Rho,nRho,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,Rho,nRho,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type1) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,Rho,nRho,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,Rho,nRho,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,Rho,nRho,nD,T_Rho)
      Else If (Functional_type.eq.meta_GGA_type2) Then
         Dens_I=Dens_I+Compute_Rho (Weights,mGrid,Rho,nRho,nD,T_Rho)
         Grad_I=Grad_I+Compute_Grad(Weights,mGrid,Rho,nRho,nD,T_Rho)
         Tau_I =Tau_I +Compute_Tau (Weights,mGrid,Rho,nRho,nD,T_Rho)
      End If
C     Write (*,*) Dens_I,Grad_I,Tau_I
*                                                                      *
************************************************************************
*                                                                      *
*        For New Functionals which use P2_OnTop (like NEWF)
      If (DFTFOCK.eq.'DIFF'.and.nD.eq.2) Then
*        For "Normal" Functionals
         If (Do_MO.and..not.(Functional_type.eq.CASDFT_Type)) Then
            Call Do_RhoIA(nRho,mGrid,Rho,Work(ipRhoI),
     &                    Work(ipRhoA),mRho,TabMO,mAO,
     &                    nMOs,D1mo,nd1mo)
         Else
            Call WarningMessage(2,' We need MO for RhoIA')
            Call Abend()
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QExit('Rho')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-- (A.Ohrn): Here I add the routine which constructs the kernel for
*   the Xhole application. A bit 'cheating' but hey what da hey!
*
      If(l_Xhol) then
        Call Xhole(nRho,mGrid,Rho,Grid,mAO,nMOs,TabMO,ndF_dRho,nD,
     &             dF_dRho,Weights,ip_OrbDip,Func)
        Go To 1979
      Endif
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QEnter('Functional')
#endif
************************************************************************
*                                                                      *
*---- Evaluate the functional on the grid                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(dF_dRho,ndF_dRho*mGrid)
      Call FZero(F_xc,mGrid)
      Call GetMem('F_xca','Allo','Real',ip_F_xca,mGrid)
      Call GetMem('F_xcb','Allo','Real',ip_F_xcb,mGrid)
      Call GetMem('tmpB','Allo','Real',ip_tmpB,mGrid)
      Call FZero(Work(ip_F_xca),mGrid)
      Call FZero(Work(ip_F_xcb),mGrid)
      Call FZero(Work(ip_tmpB),mGrid)
*
*1)   evaluate the energy density, the derivative of the functional with
*     respect to rho and grad rho.
*
      Call Kernel(mGrid,Rho,nRho,P2_ontop,
     &            nP2_ontop,nD,F_xc,dF_dRho,
     &            ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_Rho)
cGLM     &            nP2_ontop,nD,F_xc,F_xca,F_xcb,dF_dRho,
cGLM     &            ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_Rho,Work(ip_tmpB))
*
      If (nD.eq.2.and.DFTFOCK.eq.'DIFF') Then
         Call P2Diff(mGrid,dF_dRho,ndF_dRho,Rho,nRho)
      End If
*
*     Integrate the energy of the functional
*
      Func=Func+DDot_(mGrid,Weights,1,F_xc,1)
cGLM     write(6,*) 'Func in do_batch =', Func
      Funcaa=Funcaa+DDot_(mGrid,Weights,1,Work(ip_F_xca),1)
      Funcbb=Funcbb+DDot_(mGrid,Weights,1,Work(ip_F_xcb),1)
      Funccc=Funccc+DDot_(mGrid,Weights,1,Work(ip_tmpB),1)
         call xflush(6)
      Call GetMem('tmpB','Free','Real',ip_tmpB,mGrid)
      Call GetMem('F_xcb','Free','Real',ip_F_xcb,mGrid)
      Call GetMem('F_xca','Free','Real',ip_F_xca,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
      Call QExit('Functional')
#endif
*                                                                      *
************************************************************************
*                                                                      *
1979  Continue  !Jump here and skip the call to the kernel.

      If (.Not.Do_Grad) Then
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _TIME_
         Call QEnter('Integral')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the DFT contribution to the Fock matrix                  *
*                                                                      *
************************************************************************
*                                                                      *
         If (Functional_type.eq.LDA_Type) Then
*
*2)         form matrix elements for the potential, from derivative of
*           the functional with  respect to rho, respectively.
*
            If (nD.eq.2.and.DFTFOCK.eq.'DIFF') Then
*
               If(Do_TwoEl) Call Calc_PUVX2(TmpPUVX,nTmpPUVX,TabMO,mAO,
     &                                   mGrid,nMOs,dF_dRho,
     &                                   ndF_dRho,nD,Weights)
*
               Call ConvdFdRho(mGrid,dF_dRho,ndF_dRho,
     &                      Work(ipRhoI),Work(ipRhoA),mRho)
               Call Do_Int_CASDFT2(Weights,mGrid,mAO,
     &                             TabSO,nMOs,TabMO,nMOs,
     &                             nD,FckInt,nFckInt,nFckDim,
     &                             Work(ipRhoI),Work(ipRhoA),mRho,
     &                             dF_dRho,ndF_dRho,
     &                             dF_dP2ontop,ndF_dP2ontop,
     &                             TmpPUVX,nTmpPUVX)
*
            Else
!First, calculate some sizes:
             NFINT=nTmpPUVX
             NTOT1=nFckInt

             If(KSDFA(1:5).eq.'TLSDA') then
               If(do_pdftPot) then
!               CALL GETMEM('OE_OT','ALLO','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','ALLO','REAL',LTEG_DB,NFINT)

!               CALL DCOPY_(NTOT1,0.0D0,0,WORK(LOE_DB),1)!NTOT1
!               CALL DCOPY_(NFINT,0.0D0,0,WORK(LTEG_DB),1)

!               Call Get_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Get_dArray('ONTOPT',work(LTEG_DB),NFINT)

               Call Calc_OTPUVX(Work(LTEG_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         D1MO,nD1MO,nsym)

!
               Call Calc_OTOE(Work(LOE_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         nsym)

!               Call Put_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Put_dArray('ONTOPT',work(LTEG_DB),NFINT)

!               CALL GETMEM('OE_OT','FREE','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','FREE','REAL',LTEG_DB,NFINT)
               end if
             else If(KSDFA(1:6).eq.'FTLSDA') then
               If(do_pdftPot) then

!               CALL GETMEM('OE_OT','ALLO','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','ALLO','REAL',LTEG_DB,NFINT)

!               CALL DCOPY_(NTOT1,0.0D0,0,WORK(LOE_DB),1)!NTOT1
!               CALL DCOPY_(NFINT,0.0D0,0,WORK(LTEG_DB),1)

!               Call Get_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Get_dArray('ONTOPT',work(LTEG_DB),NFINT)


               Call Calc_OTPUVX_ftlsda2(Work(LTEG_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         D1MO,nD1MO,nsym)

!               tmpor = Work(loe_db)
!
               Call Calc_OTOEf(Work(LOE_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         nsym)

!               Call Put_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Put_dArray('ONTOPT',work(LTEG_DB),NFINT)

               Call xflush(6)
!               CALL GETMEM('OE_OT','FREE','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','FREE','REAL',LTEG_DB,NFINT)
               end if
             end if

             If(KSDFA(1:5).ne.'TLSDA'.and.KSDFA(1:6).ne.'FTLSDA') then
                 Call DFT_Int(Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                      FckInt,nFckInt,SOTemp,nSOTemp,
     &                      TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                      nSym,nD,Flop,Rho,nRho,Work(ipTmp),nTmp,
     &                      Work(ip_Fact),ndc,mAO,Work(ipTabAOMax),T_Y,
     &                      list_bas,Functional_type,nAOMax)
             End If
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
           If (nD.eq.2.and.DFTFOCK.eq.'DIFF') Then
***************
              If (Do_TwoEl) Call Calc_PUVX2(TmpPUVX,nTmpPUVX,TabMO,
     &                                   mAO,mGrid,nMOs,dF_dRho,
     &                                   ndF_dRho,nD,Weights)
              If (.not.(Functional_Type.eq.CASDFT_Type))
     &        Call ConvdFdRho(mGrid,dF_dRho,ndF_dRho,
     &                             Work(ipRhoI),Work(ipRhoA),mRho)
              Call Do_Int_CASDFT2(Weights,mGrid,mAO,
     &                            TabSO,nMOs,TabMO,nMOs,
     &                            nD,FckInt,nFckInt,nFckDim,
     &                            Work(ipRhoI),Work(ipRhoA),mRho,
     &                            dF_dRho,ndF_dRho,
     &                            dF_dP2ontop,ndF_dP2ontop,
     &                            TmpPUVX,nTmpPUVX)
******************
           Else
             NFINT=nTmpPUVX
             NTOT1=nFckInt

             If(KSDFA(1:4).eq.'TPBE'.or.
     &               KSDFA(1:5).eq.'TOPBE'.or.
     &               KSDFA(1:5).eq.'TBLYP'.or.
     &               KSDFA(1:7).eq.'TREVPBE') then

               If(do_pdftPot) then
!               CALL GETMEM('OE_OT','ALLO','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','ALLO','REAL',LTEG_DB,NFINT)

!               CALL DCOPY_(NTOT1,0.0D0,0,WORK(LOE_DB),1)!NTOT1
!               CALL DCOPY_(NFINT,0.0D0,0,WORK(LTEG_DB),1)

!               Call Get_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Get_dArray('ONTOPT',work(LTEG_DB),NFINT)

               Call Calc_OTPUVXGGA_2(Work(LTEG_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         D1MO,nD1MO,nsym)

               Call Calc_OTOEGGA(Work(LOE_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         nsym)

!               Call Put_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Put_dArray('ONTOPT',work(LTEG_DB),NFINT)

!               CALL GETMEM('OE_OT','FREE','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','FREE','REAL',LTEG_DB,NFINT)
              end if
             Else If(KSDFA(1:5).eq.'FTPBE'.or.
     &               KSDFA(1:6).eq.'FTOPBE'.or.
     &               KSDFA(1:6).eq.'FTBLYP'.or.
     &               KSDFA(1:8).eq.'FTREVPBE') then
               If(do_pdftPot) then
!               CALL GETMEM('OE_OT','ALLO','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','ALLO','REAL',LTEG_DB,NFINT)

!               CALL DCOPY_(NTOT1,0.0D0,0,WORK(LOE_DB),1)!NTOT1
!               CALL DCOPY_(NFINT,0.0D0,0,WORK(LTEG_DB),1)

!               Call Get_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Get_dArray('ONTOPT',work(LTEG_DB),NFINT)

               Call Calc_OTPUVXGGA_ft(Work(LTEG_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         D1MO,nD1MO,nsym)

               Call Calc_OTOEGGA_ft(Work(LOE_DB),TabMO,mAO,mGrid,
     &         nMOs,P2_ontop,nP2_ontop,Rho,nRho,dF_dRho,
     &         ndF_dRho,Work(ipRhoI),Work(ipRhoA),mRho,Weights,
     &         nsym)

!               Call Put_dArray('ONTOPO',work(LOE_DB),NTOT1)
!               Call Put_dArray('ONTOPT',work(LTEG_DB),NFINT)

!               CALL GETMEM('OE_OT','FREE','REAL',LOE_DB,NTOT1)
!               CALL GETMEM('TEG_OT','FREE','REAL',LTEG_DB,NFINT)
              end if
             end if
!FIND NTOT1,NFINT equivalents

!               Call DFT_Int(Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
!     &                      FckInt,nFckInt,SOTemp,nSOTemp,
!     &                      TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
!     &                      nSym,nD,Flop,Rho,nRho,Work(ipTmp),nTmp,
!     &                      Work(ip_Fact),ndc,mAO,Work(ipTabAOMax),T_Y,
!     &                      list_bas,Functional_type,nAOMax)
             If(.not.l_casdft) then
               Call DFT_Int(Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                      FckInt,nFckInt,SOTemp,nSOTemp,
     &                      TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                      nSym,nD,Flop,Rho,nRho,Work(ipTmp),nTmp,
     &                      Work(ip_Fact),ndc,mAO,Work(ipTabAOMax),T_Y,
     &                      list_bas,Functional_type,nAOMax)
             end if
           End If
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
          Call DFT_Int(Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                 FckInt,nFckInt,SOTemp,nSOTemp,
     &                 TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                 nSym,nD,Flop,Rho,nRho,Work(ipTmp),nTmp,
     &                 Work(ip_Fact),ndc,mAO,Work(ipTabAOMax),T_Y,
     &                 list_bas,Functional_type,nAOMax)
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
#ifdef _TIME_
         Call QExit('Integral')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*    Compute the DFT contribution to the gradient                      *
*                                                                      *
************************************************************************
*                                                                      *
      Else
*
         Call DFT_Grad(Grad,nGrad,dF_dRho,ndF_dRho,nD,Grid,mGrid,
     &                 dRho_dR,ndRho_dR,nGrad_Eff,Rho,nRho,IndGrd,
     &                 Weights,iTab,Temp,F_xc,dW_dR,iChBas,MxFnc,iNQ)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If(mRho.ne.-1) Then
        Call GetMem('Rho_I','Free','Real',ipRhoI,mGrid*mRho)
        Call GetMem('Rho_A','Free','Real',ipRhoA,mGrid*mRho)
      End If
#ifdef _DEBUG_
      Debug=Debug_Save
#endif
#ifdef _TIME_
      Call qExit('Do_Batch ')
#endif
      Return
      End
