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
* Copyright (C) 2001, Roland Lindh                                     *
************************************************************************
      Subroutine DrvNQ(Kernel,FckInt,nFckDim,Func,
     &                 Density,nFckInt,nD,
     &                 Do_Grad,Grad,nGrad,
     &                 Do_MO,Do_TwoEl,DFTFOCK)
************************************************************************
*                                                                      *
*     Predriver for numerical integration utility.                     *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             December 2001                                            *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "real.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "setup.fh"
#include "info.fh"
#include "nsd.fh"
#include "debug.fh"
#include "grid_on_disk.fh"
#include "status.fh"
#include "ksdft.fh"
      Real*8 FckInt(nFckInt,nFckDim),Density(nFckInt,nD), Grad(nGrad)
      Logical Do_Grad, Do_MO,Do_TwoEl,PMode
      Logical l_XHol, l_casdft
      Character*4 DFTFOCK
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
#include "nq_structure.fh"
      declare_ip_r_quad
      declare_ip_angular
*                                                                      *
************************************************************************
*                                                                      *
C     Call QEnter('DrvNQ')
*                                                                      *
************************************************************************
*                                                                      *
      RMx=0.0D0
      RNx=0.0D0

      If (Functional_Type.eq.CASDFT_Type) Then
      Do_TwoEl        =.True.
      End If
*
      If (Do_TwoEl) Then
      Do_MO           =.True.
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-- Additions for Xhole. DFTFOCK is here reset to SCF, but we
*   define a new logical variable. This is no DFT-integrations,
*   just dummy-stuff for getting the numerical LoProp-
*   integrations. (A.Ohrn)
*
      If (DFTFOCK(1:4).eq.'XHOL') then
      Write(DFTFOCK,'(A)')'SCF '
      l_Xhol=.true.
      Else
      l_Xhol=.false.
      Endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Allocate enough memory for Maps2p
*
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nShell)
      Call GetMem('s2p','Allo','Inte',ips2p,nShell*nIrrep)
      Call GetMem('R_Min','Allo','Real',ipR_Min,LMax_NQ+1)
*
        NQ_Status=Inactive
      Call Setup_NQ(iWork(ips2p),nShell,nIrrep,nNQ,Do_Grad,Do_MO,nD,
     &              PThr,PMode,Work(ipR_Min),LMax_NQ)
*
      Call GetMem('R_Min','Free','Real',ipR_Min,LMax_NQ)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate scratch memory for the temporary AO matrix
*     elements
*
      nAOInt=0
      Do iShell = 1, nShell
         iBas = iSD(3,iShell)
         iCmp = iSD(2,iShell)
         nAOInt=Max(nAOInt,iBas*iCmp)
      End Do
*
      Call GetMem('AOInt','Allo','Real',ipAOInt,nD*nAOInt**2)
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate scratch memory for the temporary SO matrix
*     elements
*
      iSmLbl=1
      nSOTemp=0
      Do iSkal = 1, nShell
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iShell= iSD(11,iSkal)
         Do jSkal = 1, iSkal
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jShell= iSD(11,jSkal)
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            nSOTemp=Max(nSOTemp,iBas*jBas*nSO)
         End Do
      End Do
      Call GetMem('SO_Temp','Allo','Real',ipSOTemp,nSOTemp*nD)
*                                                                      *
************************************************************************
*                                                                      *
*-----Allocate memory sufficiently large to store all grid points
*     and associated weights.
*
      Call GetMem('Grid','Allo','Real',ip_Grid,nGridMax*3)
*                                                                      *
************************************************************************
*                                                                      *
*     CASDFT stuff:
*
      l_casdft=.false.
      nP2=1
      nCmo=1
      nD1mo=1
      nMos=1
      nTmpTUVX=1
      nTmpPUVX=1
      ipP2mo=ip_Dummy
      ipCmo=ip_Dummy
      ipD1mo=ip_Dummy
      ipDoIt=ip_iDummy
      ipTmpPUVX=ip_Dummy
      ipTmpTUVX=ip_Dummy
*
      NQNAC=0
      If (DFTFOCK.ne.'SCF ') Then
         Do iIrrep = 0, mIrrep - 1
            NQNAC = NQNAC + nAsh(iIrrep)
         End Do
      End If
*
      LuGridFile=31
      LuGridFile=IsFreeUnit(LuGridFile)
      Call Molcas_Open(LuGridFile,'GRIDFILE')
************************************************************************
*
************************************************************************
*                                                                      *
*     Definition of resources needed for the functionals.              *
*                                                                      *
*     mAO: the number of derivatives needed of an basis function.      *
*          Depending of the functional type and if gradients will be   *
*          computed. Numbers will be 1, 4, 10, 20, 35, etc.            *
*     nRho:the number of parameters of the functional. Note that this  *
*          is different for the same functional depending on if it is  *
*          a closed or open-shell case.                                *
*     mdRho_dR: number of derivatives of the parameters with respect   *
*          to the nuclear coordinates. The true number is of course    *
*          three (x,y,z) times this.                                   *
*     nF_drho: the number of derivatives of the functional wrt the     *
*          parameters. Note that grad rho is not a direct parameter    *
*          but that we use gamma.                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
*                                                                      *
************************************************************************
*                                                                      *
*        We need the AOs, for gradients we need the derivatives too.
*
         mAO=1
         If (Do_Grad) mAO=4
*
*        We need rho.
*        For gradients we need derivatives of rho wrt the coordinates
*
         nRho=1*nD
         mdRho_dr=0
         If (Do_Grad) mdRho_dr=nRho
*
*        We need derivatives of the functional with respect to
*        rho(alpha). In case of open-shell calculations we also
*        need rho(beta).
*
         If (nD.eq.1) Then
            ndF_dRho=1
         Else
            ndF_dRho=2
         End If
         nP2_ontop=1
         ndF_dP2ontop=1
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.GGA_type) Then
*                                                                      *
************************************************************************
*                                                                      *
*        We need the AOs and their derivatives, for  gradients we need
*        the second derivatives too.
*
         mAO=4
         If (Do_Grad) mAO=10
*
*        We need rho and grad rho
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=4*nD
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=nRho
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha) and gamma(alpha,beta).
*        In case of open-shell calculations we also
*        need rho(beta) and gamma(beta,beta).
*
         If (nD.eq.1) Then
            ndF_dRho=3
         Else
            ndF_dRho=5
         End If
         nP2_ontop=4
         ndF_dP2ontop=4
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.meta_GGA_type1) Then
*                                                                      *
************************************************************************
*                                                                      *
*        We need the AOs and their derivatives, for  gradients we need
*        the second derivatives too.
*
         mAO=4
         If (Do_Grad) mAO=10
*
*        We need rho, grad rho and tau.
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=5*nD
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=nRho
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha), gamma(alpha,beta) and
*        tau(alpha). In case of open-shell calculations we also
*        need rho(beta), gamma(beta,beta) and tau(beta).
*
         If (nD.eq.1) Then
            ndF_dRho=4
         Else
            ndF_dRho=7
         End If
         nP2_ontop=4
         ndF_dP2ontop=4
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.meta_GGA_type2) Then
*                                                                      *
************************************************************************
*                                                                      *
*        Not debugged yet!
*        We need the AOs and their 1st and 2nd derivatives, for
*        gradients we need the 3rd order derivatives too.
*
         mAO=10
         If (Do_Grad) mAO=20
*
*        We need rho, grad rho, tau, and the Laplacian
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=6*nD
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=nRho
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha), gamma(alpha,beta),
*        tau(alpha) and laplacian(alpha). In case of open-shell
*        calculations we also need rho(beta), gamma(beta,beta),
*        tau(beta) and laplacian(beta).
*
         If (nD.eq.1) Then
            ndF_dRho=5
         Else
            ndF_dRho=9
         End If
         nP2_ontop=4
         ndF_dP2ontop=4
c         Call WarningMessage(2,
c     &        'Meta-GGA functional type 2 not fully DEBUGGED yet!')
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.PAM_type) Then
*                                                                      *
************************************************************************
*                                                                      *
         mAO=5
         nRho=4*nD
         mdRho_dR=0
         If (nD.eq.1) Then
            ndF_dRho=3
         Else
            ndF_dRho=5
         Endif
         nP2_ontop=4
         ndF_dP2ontop=4
         Call WarningMessage(2,'PAM functionals not implemented yet!')
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.CS_type) Then
*                                                                      *
************************************************************************
*                                                                      *
         mAO=5
         nRho=4*nD
         mdRho_dR=0
         If (nD.eq.1) Then
            ndF_dRho=3
         Else
            ndF_dRho=5
         Endif
         nP2_ontop=4
         ndF_dP2ontop=4
         Call WarningMessage(2,'CS functionals not implemented yet!')
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.CASDFT_type) Then
*
*        nD's definition is not consistent with the use here!
*        This needs to be restructured.
*
         mAO=10
         mdRho_dR=0
         If (Do_Grad) Then
             Call WarningMessage(2,'CASDFT: Gradient not available.')
             Call Abend()
         End If
         nRho=4*nD
*
         If (nD.eq.1) Then
            ndF_dRho=3 ! could be 2?
         Else
            ndF_dRho=5
         Endif
         nP2_ontop=6
         ndF_dP2ontop=6
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Functional_type=Other_type
         Call WarningMessage(2,'DrvNQ: Invalid Functional_type!')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('F_xc','Allo','Real',ip_F_xc,nGridMax)
cGLM      Call GetMem('F_xca','Allo','Real',ip_F_xca,nGridMax)
cGLM      Call GetMem('F_xcb','Allo','Real',ip_F_xcb,nGridMax)
      Call GetMem('Rho','Allo','Real',ip_Rho,nRho*nGridMax)
      Call GetMem('dF_dRho','Allo','Real',ip_dFdRho,ndF_dRho*nGridMax)
*
      Call GetMem('Weights','Allo','Real',ip_Weights,nGridMax)
      Call GetMem('list_s','Allo','Inte',iplist_s,2*nIrrep*nShell)
      Call GetMem('list_exp','Allo','Inte',iplist_exp,3*nIrrep*nShell)
      iplist_bas=iplist_exp+nIrrep*nShell
      Call GetMem('list_p','Allo','Inte',iplist_p,nNQ)
      Call GetMem('R2_trail','Allo','Real',ipR2_trail,nNQ)
c      Call GetMem('tmpB','Allo','Real',ip_tmpB,nGridMax)
*                                                                      *
************************************************************************
* Global variable for MCPDFT functionals                               *
      l_casdft = KSDFA(1:5).eq.'TLSDA'   .or.
     &           KSDFA(1:6).eq.'TLSDA5'  .or.
     &           KSDFA(1:5).eq.'TBLYP'   .or.
     &           KSDFA(1:6).eq.'TSSBSW'  .or.
     &           KSDFA(1:5).eq.'TSSBD'   .or.
     &           KSDFA(1:5).eq.'TS12G'   .or.
     &           KSDFA(1:4).eq.'TPBE'    .or.
     &           KSDFA(1:5).eq.'FTPBE'   .or.
     &           KSDFA(1:5).eq.'TOPBE'   .or.
     &           KSDFA(1:6).eq.'FTOPBE'  .or.
     &           KSDFA(1:7).eq.'TREVPBE' .or.
     &           KSDFA(1:8).eq.'FTREVPBE'.or.
     &           KSDFA(1:6).eq.'FTLSDA'  .or.
     &           KSDFA(1:6).eq.'FTBLYP'
      if(Debug) write(6,*) 'l_casdft value at drvnq.f:',l_casdft
      if(Debug.and.l_casdft) write(6,*) 'MCPDFT with functional:', KSDFA
************************************************************************
      If (Do_MO) Then
         If (NQNAC.ne.0) Then
           If(.not.l_casdft) Then
             Call Get_D1MO(ipD1mo,nd1mo)
             Call Get_P2mo(ipP2mo,nP2)
           End If
         End If
         Call Get_CMO(ipCmo,nCmo)
         Call Get_iArray('nAsh',nAsh,mIrrep)
         nMOs=0
         Do iIrrep = 0, mIrrep-1
            nMOs=nMOs+mBas(iIrrep)
         End Do
         Call GetMem('DoIt','Allo','Inte',ipDoIt,nMOs)
         iMO=ipDoIt-1
         Do iIrrep = 0, mIrrep-1
            Do jMO = 1, nISh(iIrrep)+nASh(iIrrep)
               iMO=iMO+1
               iWork(iMO)=1
            End Do
            Do jMO = 1, mBas(iIrrep)-nISh(iIrrep)-nASh(iIrrep)
               iMO=iMO+1
               iWork(iMO)=1
            End Do
         End Do
      End If
***
*     Prepare memory for two-electron integrals:
*     nPUVX, nTUVX
*
      If (Do_TwoEl) Then
         If (.not.Do_MO) Then
            Call WarningMessage(2,
     &              ' Can''t produce 2 el dft integrals without MO')
            Call Abend()
         End If
         NQNACPAR = ( NQNAC**2 + NQNAC )/2
         NQNACPR2 = ( NQNACPAR**2 + NQNACPAR )/2
         nTmpTUVX = NQNACPR2
*
         iStack = 0
         Do iIrrep = 0, mIrrep-1
           iOrb = mBas(iIrrep) - nFro(iIrrep)
           Do jIrrep = 0, mIrrep-1
             jAsh = nAsh(jIrrep)
             ijIrrep=iEor(iIrrep,jIrrep)
             Do kIrrep = 0, mIrrep-1
               kAsh = nAsh(kIrrep)
               ijkIrrep=iEor(ijIrrep,kIrrep)
               If (ijkIrrep.le.kIrrep) Then
                 lAsh = nAsh(ijkIrrep)
                 kl_Orb_pairs = kAsh*lAsh
                 If ( kIrrep.eq.ijkIrrep )
     &                  kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
                 iStack = iStack + iOrb*jAsh*kl_Orb_pairs
               End If
             End Do
           End Do
         End Do
         nTmpPUVX=iStack
*
         Call GetMem('TmpPUVX','Allo','Real',ipTmpPUVX,nTmpPUVX)
         Call dCopy_(nTmpPUVX,[0.0d0],0,Work(ipTmpPUVX),1)
      End If
*
      If (Functional_Type.eq.CASDFT_Type) Then
         Call GetMem('P2_ontop','Allo','Real',ipp2_ontop,
     &               nP2_ontop*nGridMax)
         Call GetMem('dF_dP2ontop','Allo','Real',ipdF_dP2ontop,
     &               ndF_dP2ontop*nGridMax)
      Else
         ipP2_ontop=ip_Dummy
         ipdF_dp2ontop=ip_Dummy
      Endif
*
      If (Do_Grad) Then
         Call GetMem('list_g','Allo','Inte',iplist_g,3*nShell*nIrrep)
         mGrad=3*nAtoms
         Call GetMem('IndGrd','Allo','Inte',ipIndGrd,mGrad)
         Call GetMem('iTab','Allo','Inte',ipiTab,4*mGrad)
         Call GetMem('Temp','Allo','Real',ipTemp,mGrad)
      Else
         iplist_g=ip_iDummy
         ipIndGrd=ip_iDummy
         ipiTab  =ip_iDummy
         ipTemp  =ip_Dummy
      End If

      If (.Not.Do_Grad) Call FZero(FckInt,nFckInt*nFckDim)
*                                                                      *
************************************************************************
*                                                                      *
      Thr=Threshold
      if(Debug) write(6,*) 'l_casdft value at drvnq.f:',l_casdft
      if(Debug.and.l_casdft) write(6,*) 'MCPDFT with functional:', KSDFA
      If(l_casdft) then
        NQNAC=0
        Call Get_iArray('nAsh',nAsh,mIrrep)
        Do iIrrep = 0, mIrrep - 1
          NQNAC = NQNAC + nAsh(iIrrep)
        End Do
        IF(NQNAC.ne.0) then
          NQNACPAR = ( NQNAC**2 + NQNAC )/2
          NQNACPR2 = ( NQNACPAR**2 + NQNACPAR )/2
          Call Get_D1MO(ipD1mo,nd1mo)
cGLM          write(6,*) 'D1MO in drvNQ routine'
cGLM          write(6,*) (Work(ipD1mo+i), i=0,NQNACPAR-1)
          call Get_P2mo(ipP2mo,nP2)
cGLM          write(6,*) 'P2MO in drvNQ routine'
cGLM          write(6,*) (Work(ipP2mo+i), i=0,NQNACPR2-1)
        END IF
         Call GetMem('P2_ontop','Allo','Real',ipp2_ontop,
     &               nP2_ontop*nGridMax)
         Call GetMem('dF_dP2ontop','Allo','Real',ipdF_dP2ontop,
     &               ndF_dP2ontop*nGridMax)
        Call dCopy_(nP2_ontop*nGridMax,[0.0d0],0,Work(ipp2_ontop),1)
        Call dCopy_(ndF_dP2ontop*nGridMax,[0.0d0],0,
     &                                    Work(ipdF_dP2ontop),1)

      end if

      Call DrvNQ_(Kernel,Func,
     &            iWork(ips2p),nIrrep,
     &            iWork(iplist_s),iWork(iplist_exp),iWork(iplist_bas),
     &            nShell,iWork(iplist_p),Work(ipR2_trail),nNQ,
     &            Work(ipAOInt),nAOInt,FckInt,nFckDim,
     &            Density,nFckInt,nD,
     &            Work(ipSOTemp),nSOTemp,
     &            Work(ip_Grid),Work(ip_Weights),Work(ip_Rho),
     &            nGridMax,nRho,
     &            ndF_dRho,nP2_ontop,ndF_dP2ontop,
     &            Do_Mo,Do_TwoEl,l_Xhol,
     &            Work(ipTmpPUVX),nTmpPUVX,
     &            nMOs,Work(ipCMO),nCMO,
     &            iWork(ipDoIt),
     &            Work(ipP2mo),nP2,Work(ipD1mo),nd1mo,Work(ipp2_ontop),
     &            Do_Grad,Grad,nGrad,iWork(iplist_g),
     &            iWork(ipIndGrd),iWork(ipiTab),Work(ipTemp),mGrad,
     &            Work(ip_F_xc),
cGLM     &        Work(ip_F_xca),Work(ip_F_xcb),
     &            Work(ip_dFdRho),work(ipdF_dP2ontop),
     &            DFTFOCK,mAO,mdRho_dR)
*                                                                      *
************************************************************************
*                                                                      *
*-----Deallocate the memory
*
      Call GetMem('O','Free','Real',ip_O,3*3)
      If (Do_Grad) Then
         Call GetMem('Temp','Free','Real',ipTemp,mGrad)
         Call GetMem('iTab','Free','Inte',ipiTab,4*mGrad)
         Call GetMem('IndGrd','Free','Inte',ipIndGrd,mGrad)
         Call GetMem('list_g','Free','Inte',iplist_g,3*nShell*nIrrep)
      End If
      Call GetMem('R2_trail','Free','Real',ipR2_trail,nNQ)
      Call GetMem('list_p','Free','Inte',iplist_p,nNQ)
      Call GetMem('list_exp','Free','Inte',iplist_exp,3*nIrrep*nShell)
      Call GetMem('list_s','Free','Inte',iplist_s,2*nIrrep*nShell)
      Call GetMem('Weights','Free','Real',ip_Weights,nGridMax)
      Call GetMem('dF_dRho','Free','Real',ip_dFdRho,ndF_dRho*nGridMax)
*Do_TwoEl
      If(ipP2mo.ne.ip_Dummy) Call Free_Work(ipP2mo)
      If(ipD1MO.ne.ip_Dummy) Call Free_Work(ipD1MO)
      If(ipCMO.ne.ip_Dummy)  Call Free_Work(ipCMO)
      If(ipDoIt.ne.ip_iDummy) Call GetMem('DoIt','Free','Inte',
     &                                    ipDoIt,nMOs)
      If(ipTmpPUVX.ne.ip_Dummy) Then
         Call Put_dArray('DFT_TwoEl',Work(ipTmpPUVX),nTmpPUVX)
         Call GetMem('TmpPUVX','Free','Real',ipTmpPUVX,nTmpPUVX)
      End If
*
      Call GetMem('Rho','Free','Real',ip_Rho,nRho*nGridMax)
      Call GetMem('F_xc','Free','Real',ip_F_xc,nGridMax)
cGLM      Call GetMem('F_xca','Free','Real',ip_F_xca,nGridMax)
cGLM      Call GetMem('F_xcb','Free','Real',ip_F_xcb,nGridMax)
      Call GetMem('Grid','Free','Real',ip_Grid,nGridMax*3)
c      Call GetMem('tmpB','Free','Real',ip_tmpB,nGridMax)

      if(Debug) write(6,*) 'l_casdft value at drvnq.f:',l_casdft
      if(Debug.and.l_casdft) write(6,*) 'MCPDFT with functional:', KSDFA
      If (Functional_type.eq.CASDFT_Type.or.l_casdft) Then
         Call GetMem('P2_ontop','Free','Real',ipP2_ontop,
     &               nP2_ontop*nGridMax)
         Call GetMem('dF_dP2ontop','Free','Real',ipdF_dP2ontop,
     &               ndF_dP2ontop*nGridMax)
      End If
*
      Do iNQ = 1, nNQ
         ip_iRx=ip_of_iWork_d(Work(ip_R_Quad(iNQ)))
         ip_Rx=iWork(ip_iRx)
         Call GetMem('Radial','Free','Real',ip_Rx,iDum)
         ip_iA=ip_of_iWork_d(Work(ip_Angular(iNQ)))
         ip_A=iWork(ip_iA)
         Call GetMem('ip_Angular','Free','Inte',ip_A,iDum)
      End Do
      Call GetMem('NumRadEff','Free','Inte',ip_nR_eff,nNQ)
      Call GetMem('Coor','FREE','REAL',ipCoor,3*8*nAtoms)

      Call GetMem('SO_Temp','Free','Real',ipSOTemp,nSOTemp)
      Call GetMem('AOInt','Free','Real',ipAOInt,nD*nAOInt**2)
      Call GetMem('nq_centers','Free','Real',ipNQ,nShell*l_NQ)
      Call GetMem('nMem','Free','Real',ipMem,nMem)
      Call GetMem('Tmp','Free','Real',ipTmp,nTmp)
      Call GetMem('Dijs','Free','Real',ipDijs,MxDij)
      Call Free_Work(ip_Fact)
      Call GetMem('s2p','Free','Inte',ips2p,nshell)
      NQ_Status=Inactive
*                                                                      *
************************************************************************
*                                                                      *
*---- Write the status flag and TOC.
*
      If (iGrid_Set.eq.Intermediate .and.
     &     Grid_Status.eq.Regenerate) iDisk_Set(Final)=iDisk_Grid
      G_S(iGrid_Set)=Use_Old
*
      iDisk_Grid=0
      Call iDaFile(Lu_Grid,1,G_S,5,iDisk_Grid)
      iDisk_Grid=iDisk_Set(iGrid_Set)
      Call iDaFile(Lu_Grid,1,iWork(ip_GridInfo),
     &             2*number_of_subblocks,iDisk_Grid)
*
      Call DaClos(Lu_Grid)
*
      Call Free_iWork(ip_GridInfo)
*                                                                      *
************************************************************************
*                                                                      *
      Call IniPkR8(PThr,PMode)
*
      Call xFlush(LuGridFile)
      Close(LuGridFile)
C     Call QExit('DrvNQ')
      Return
      End
