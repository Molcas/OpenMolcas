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
      Subroutine DrvNQ(Kernel,FckInt,nFckDim,Funct,
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
      use Symmetry_Info, only: nIrrep
      use KSDFT_Info, only: KSDFA
      use nq_Grid, only: Rho, GradRho, Sigma, Tau, Lapl
      use nq_Grid, only: vRho, vSigma, vTau, vLapl
      use nq_Grid, only: Grid, Weights
      use nq_Grid, only: nRho, nGradRho, nTau, nSigma, nLapl, nGridMax
      use nq_Grid, only: l_CASDFT, kAO
      use nq_Grid, only: F_xc, F_xca, F_xcb
      use nq_Grid, only: Coor, R2_trial, Pax
      use nq_pdft, only: lft, lGGA
      use nq_MO, only: DoIt, CMO, D1MO, P2MO, P2_ontop
      use libxc
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "setup.fh"
#include "nsd.fh"
#include "debug.fh"
#include "grid_on_disk.fh"
#include "status.fh"
#include "ksdft.fh"
      Real*8 FckInt(nFckInt,nFckDim),Density(nFckInt,nD), Grad(nGrad)
      Logical Do_Grad, Do_MO,Do_TwoEl,PMode
      Character*4 DFTFOCK
      Integer nBas(8), nDel(8)
      Integer, Allocatable:: Maps2p(:,:), List_s(:,:), List_Exp(:),
     &                       List_Bas(:,:), List_P(:), List_G(:,:)
      Integer, Allocatable:: IndGrd(:), iTab(:,:)
      Real*8, Allocatable:: R_Min(:), Temp(:)
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
*                                                                      *
************************************************************************
*                                                                      *
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
*-----Allocate enough memory for Maps2p
*
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nShell)
      Call mma_allocate(Maps2p,nShell,nIrrep,Label='Maps2p')
      Call mma_allocate(R_Min,LMax_NQ+1,Label='R_Min')
*
      NQ_Status=Inactive
      Call Setup_NQ(Maps2p,nShell,nIrrep,nNQ,Do_Grad,Do_MO,nD,
     &              PThr,PMode,R_Min,LMax_NQ)
*
      Call mma_deallocate(R_Min)
*                                                                      *
************************************************************************
*                                                                      *
*-----Allocate memory sufficiently large to store all grid points
*     and associated weights.
*
      Call mma_Allocate(Grid,3,nGridMax,Label='Grid')
      Call mma_Allocate(Weights,nGridMax,Label='Weights')
*                                                                      *
************************************************************************
*                                                                      *
*     CASDFT stuff:
*
      nTmpPUVX=1
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
*                                                                      *
************************************************************************
* Global variable for MCPDFT functionals                               *

      l_casdft = KSDFA.eq.'TLSDA'   .or.
     &           KSDFA.eq.'TLSDA5'  .or.
     &           KSDFA.eq.'TBLYP'   .or.
     &           KSDFA.eq.'TSSBSW'  .or.
     &           KSDFA.eq.'TSSBD'   .or.
     &           KSDFA.eq.'TS12G'   .or.
     &           KSDFA.eq.'TPBE'    .or.
     &           KSDFA.eq.'FTPBE'   .or.
     &           KSDFA.eq.'TOPBE'   .or.
     &           KSDFA.eq.'FTOPBE'  .or.
     &           KSDFA.eq.'TREVPBE' .or.
     &           KSDFA.eq.'FTREVPBE'.or.
     &           KSDFA.eq.'FTLSDA'  .or.
     &           KSDFA.eq.'FTBLYP'

      lft      = KSDFA.eq.'FTPBE'   .or.
     &           KSDFA.eq.'FTOPBE'  .or.
     &           KSDFA.eq.'FTREVPBE'.or.
     &           KSDFA.eq.'FTLSDA'  .or.
     &           KSDFA.eq.'FTBLYP'

      if(Debug) write(6,*) 'l_casdft value at drvnq.f:',l_casdft
      if(Debug.and.l_casdft) write(6,*) 'MCPDFT with functional:', KSDFA
************************************************************************
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
         kAO=mAO
         If (Do_Grad) mAO=4
*
*        We need rho.
*        For gradients we need derivatives of rho wrt the coordinates
*
         nRho=nD
         mdRho_dr=0
         If (Do_Grad) mdRho_dr=nD
         nSigma=0
         nGradRho=0
         nLapl=0
         nTau=0
*
*        We need derivatives of the functional with respect to
*        rho(alpha). In case of open-shell calculations we also
*        need rho(beta).
*
         nP2_ontop=1
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
         kAO=mAO
         If (Do_Grad) mAO=10
*
*        We need rho and grad rho
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=nD
         nSigma=nD*(nD+1)/2
         nGradRho=nD*3
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=4*nD
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha) and gamma(alpha,beta).
*        In case of open-shell calculations we also
*        need rho(beta) and gamma(beta,beta).
*
         nP2_ontop=4
         lGGA=.True.
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
         kAO=mAO
         If (Do_Grad) mAO=10
*
*        We need rho, grad rho and tau.
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=nD
         nSigma=nD*(nD+1)/2
         nGradRho=nD*3
*        nLapl=0
         nLapl=nD
         nTau=nD
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=5*nD
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha), gamma(alpha,beta) and
*        tau(alpha). In case of open-shell calculations we also
*        need rho(beta), gamma(beta,beta) and tau(beta).
*
         nP2_ontop=4
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
         kAO=mAO
         If (Do_Grad) mAO=20
*
*        We need rho, grad rho, tau, and the Laplacian
*        For gradients we need the derrivatives wrt the coordinates
*
         nRho=nD
         nSigma=nD*(nD+1)/2
         nGradRho=nD*3
         nTau=nD
         nLapl=nD
         mdRho_dR=0
         If (Do_Grad) mdRho_dR=6*nD
*
*        We need derivatives of the functional with respect to
*        rho(alpha), gamma(alpha,alpha), gamma(alpha,beta),
*        tau(alpha) and laplacian(alpha). In case of open-shell
*        calculations we also need rho(beta), gamma(beta,beta),
*        tau(beta) and laplacian(beta).
*
         nP2_ontop=4
*                                                                      *
************************************************************************
*                                                                      *
      Else If ( Functional_type.eq.CASDFT_type) Then
*                                                                      *
************************************************************************
*                                                                      *
*        nD  definition is not consistent with the use here!
*        This needs to be restructured.
*
         mAO=10
         kAO=mAO

         nRho=nD
         nSigma=nD*(nD+1)/2
         nGradRho=nD*3
         nTau=nD
         nLapl=0
         If (Do_Grad) mdRho_dR=4*nD
         If (Do_Grad) Then
             Call WarningMessage(2,'CASDFT: Gradient not available.')
             Call Abend()
         End If
*
         nP2_ontop=6
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Functional_type=Other_type
         Call WarningMessage(2,'DrvNQ: Invalid Functional_type!')
         Call Abend()
         nRho=0
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Rho,nRho,nGridMax,Label='Rho')
      Call mma_allocate(vRho,nRho,nGridMax,Label='vRho')
      Call mma_allocate(dfunc_drho,nRho,nGridMax,Label='dfunc_drho')
      If (nSigma.ne.0) Then
         Call mma_Allocate(Sigma,nSigma,nGridMax,Label='Sigma')
         Call mma_Allocate(vSigma,nSigma,nGridMax,Label='vSigma')
         Call mma_Allocate(dfunc_dSigma,nSigma,nGridMax,
     &                     Label='dfunc_dSigma')
      End If
      If (nGradRho.ne.0) Then
         Call mma_Allocate(GradRho,nGradRho,nGridMax,Label='GradRho')
      End If
      If (nTau.ne.0) Then
         Call mma_allocate(Tau,nTau,nGridMax,Label='Tau')
         Call mma_allocate(vTau,nTau,nGridMax,Label='vTau')
         Call mma_allocate(dfunc_dTau,nTau,nGridMax,Label='dfunc_dTau')
         Tau(:,:)=Zero
      End If
      If (nLapl.ne.0) Then
         Call mma_allocate(Lapl,nLapl,nGridMax,Label='Lapl')
         Call mma_allocate(vLapl,nLapl,nGridMax,Label='vLapl')
         Call mma_allocate(dfunc_dLapl,nLapl,nGridMax,
     &                     Label='dfunc_dLapl')
         Lapl(:,:)=Zero
      End If

      Call mma_allocate(F_xc,nGridMax,Label='F_xc')
      Call mma_allocate(func,nGridMax,Label='func')
      If (l_casdft) Then
         Call mma_allocate(F_xca,nGridMax,Label='F_xca')
         Call mma_allocate(F_xcb,nGridMax,Label='F_xcb')
      End If
*
      Call mma_allocate(List_S,2,nIrrep*nShell,Label='List_S')
      Call mma_allocate(List_Exp,nIrrep*nShell,Label='List_Exp')
      Call mma_allocate(List_Bas,2,nIrrep*nShell,Label='List_Bas')
      Call mma_allocate(List_P,nNQ,Label='List_P')
      Call mma_allocate(R2_trial,nNQ,Label='R2_trial')

      If (Do_MO) Then
         If (NQNAC.ne.0) Then
           If(.not.l_casdft) Then
             NQNACPAR = ( NQNAC**2 + NQNAC )/2
             nd1mo=NQNACPAR
             Call mma_allocate(D1MO,nd1mo,Label='D1MO')
             Call Get_D1MO(D1MO,nD1MO)
             NQNACPR2 = ( NQNACPAR**2 + NQNACPAR )/2
             Call mma_Allocate(P2MO,nP2,Label='P2MO')
             Call Get_P2mo(P2MO,nP2)

           End If
         End If
         Call Get_iArray('nBas',nBas,mIrrep)
         Call Get_iArray('nDel',nDel,mIrrep)
         nCMO=0
         Do i = 1, mIrrep
            nCMO = nCMO + nBas(i)*(nBas(i)-nDel(i))
         End Do
         Call mma_allocate(CMO,nCMO,Label='CMO')
         Call Get_CMO(CMO,nCMO)
         Call Get_iArray('nAsh',nAsh,mIrrep)
         nMOs=0
         Do iIrrep = 0, mIrrep-1
            nMOs=nMOs+mBas(iIrrep)
         End Do
         Call mma_Allocate(DoIt,nMOs,Label='DoIt')
         iMO=0
         Do iIrrep = 0, mIrrep-1
            Do jMO = 1, nISh(iIrrep)+nASh(iIrrep)
               iMO=iMO+1
               DoIt(iMO)=1
            End Do
            Do jMO = 1, mBas(iIrrep)-nISh(iIrrep)-nASh(iIrrep)
               iMO=iMO+1
               DoIt(iMO)=1
            End Do
         End Do
      End If
***
*     Prepare memory for two-electron integrals:
*     nPUVX
*
      If (Do_TwoEl) Then
         If (.not.Do_MO) Then
            Call WarningMessage(2,
     &              ' Can''t produce 2 el dft integrals without MO')
            Call Abend()
         End If
         NQNACPAR = ( NQNAC**2 + NQNAC )/2
         NQNACPR2 = ( NQNACPAR**2 + NQNACPAR )/2
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
      End If
*
      If (Functional_Type.eq.CASDFT_Type) Then
         Call mma_allocate(P2_ontop,nP2_ontop,nGridMax,Label='P2_ontop')
      Endif
*
      If (Do_Grad) Then
         Call mma_allocate(List_g,3,nShell*nIrrep,Label='List_G')
         mGrad=3*nAtoms
         Call mma_allocate(IndGrd,mGrad,Label='IndGrd')
         Call mma_allocate(iTab,4,mGrad,Label='iTab')
         Call mma_allocate(Temp,mGrad,Label='Temp')
      End If

      If (.Not.Do_Grad) Call FZero(FckInt,nFckInt*nFckDim)
*                                                                      *
************************************************************************
*                                                                      *
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
          nD1MO = NQNACPAR
          Call mma_allocate(D1MO,nD1MO,Label='D1MO')
          Call Get_D1MO(D1MO,nD1MO)
          nP2 = NQNACPR2
          Call mma_Allocate(P2MO,nP2,Label='P2MO')
          call Get_P2mo(P2MO,nP2)
        END IF
        Call mma_allocate(P2_ontop,nP2_ontop,nGridMax,Label='P2_ontop')
        P2_ontop(:,:)=Zero
      end if

      Call DrvNQ_Inner(Kernel,Funct,Maps2p,nIrrep,List_S,List_Exp,
     &                 List_bas,nShell,List_P,nNQ,
     &                 FckInt,nFckDim,Density,nFckInt,nD,
     &                 nGridMax,nP2_ontop,Do_Mo,nTmpPUVX,
     &                 Do_Grad,Grad,nGrad,List_G,IndGrd,iTab,
     &                 Temp,mGrad,mAO,mdRho_dR)
*                                                                      *
************************************************************************
*                                                                      *
*-----Deallocate the memory
*
      Call mma_deallocate(Pax)
      If (Do_Grad) Then
         Call mma_deallocate(Temp)
         Call mma_deallocate(iTab)
         Call mma_deallocate(IndGrd)
         Call mma_deallocate(List_G)
      End If
      Call mma_deallocate(R2_trial)
      Call mma_deallocate(List_P)
      Call mma_deallocate(List_Bas)
      Call mma_deallocate(List_Exp)
      Call mma_deallocate(List_S)
*Do_TwoEl
      If (Allocated(D1MO)) Call mma_deallocate(D1MO)
      If (Allocated(P2MO)) Call mma_deallocate(P2MO)
      If (Allocated(CMO)) Call mma_deallocate(CMO)
      If (Allocated(DoIt)) Call mma_deallocate(DoIt)
      If (l_casdft) Then
         Call mma_deallocate(F_xcb)
         Call mma_deallocate(F_xca)
      End If
      Call mma_deallocate(func)
      Call mma_deallocate(F_xc)
*
      If (Allocated(Lapl)) Then
         Call mma_deallocate(dfunc_dLapl)
         Call mma_deallocate(vLapl)
         Call mma_deallocate(Lapl)
      End If
      If (Allocated(Tau)) Then
         Call mma_deallocate(dfunc_dTau)
         Call mma_deallocate(vTau)
         Call mma_deallocate(Tau)
      End If
      If (Allocated(GradRho)) Call mma_deallocate(GradRho)
      If (Allocated(Sigma)) Then
         Call mma_deallocate(dfunc_dSigma)
         Call mma_deallocate(vSigma)
         Call mma_deallocate(Sigma)
      End If
      Call mma_deallocate(dfunc_dRho)
      Call mma_deallocate(vRho)
      Call mma_deallocate(Rho)

      Call mma_deallocate(Weights)
      Call mma_deallocate(Grid)

      if(Debug) write(6,*) 'l_casdft value at drvnq.f:',l_casdft
      if(Debug.and.l_casdft) write(6,*) 'MCPDFT with functional:', KSDFA
      If (Functional_type.eq.CASDFT_Type.or.l_casdft) Then
         Call mma_deallocate(P2_ontop)
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
      Call mma_deallocate(Coor)

      Call GetMem('nq_centers','Free','Real',ipNQ,nShell*l_NQ)
      Call GetMem('nMem','Free','Real',ipMem,nMem)
      Call GetMem('Tmp','Free','Real',ipTmp,nTmp)
      Call Free_Work(ip_Fact)
      Call mma_deallocate(Maps2p)
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
      Return
      End
