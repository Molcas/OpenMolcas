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
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine DrvNQ_(Kernel,Func,
     &                  Maps2p,nSym,list_s,list_exp,list_bas,
     &                  nShell,list_p,R2_trial,nNQ,
     &                  AOInt,nAOInt,FckInt,nFckDim,
     &                  Density,nFckInt,nD,
     &                  SOTemp,nSOTemp,
     &                  Grid,Weights,Rho,mGrid,nRho,
     &                  ndF_dRho,nP2_ontop,ndF_dP2ontop,
     &                  Do_Mo,Do_TwoEl,l_Xhol,
     &                  TmpPUVX,nTmpPUVX,
     &                  nMOs,
     &                  CMOs,nCMO,DoIt,P2mo,np2act,D1mo,nd1mo,P2_ontop,
     &                  Do_Grad,Grad,nGrad,list_g,IndGrd,iTab,Temp,
     &                  mGrad,F_xc,dF_dRho,dF_dP2ontop,
     &                  DFTFOCK,mAO,mdRho_dR)
************************************************************************
*                                                                      *
* Object: Driver for numerical quadrature.                             *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
      External Kernel, Rsv_Tsk
#include "real.fh"
#include "WrkSpc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "itmax.fh"
#include "info.fh"
#include "nq_info.fh"
#include "grid_on_disk.fh"
#include "k2.fh"
#include "debug.fh"
#include "ksdft.fh"
      Integer Maps2p(nShell,0:nSym-1),
     &        list_s(nSym*nShell), list_exp(nSym*nShell),
     &        list_p(nNQ), DoIt(nMOs), List_g(3,nSym*nShell),
     &        IndGrd(mGrad), iTab(4,mGrad), list_bas(2,nSym*nShell)
      Real*8 AOInt(nAOInt,nAOInt,nD), FckInt(nFckInt,nFckDim),
     &       Density(nFckInt,nD),
     &       SOTemp(nSOTemp,nD), Grid(3,mGrid), Weights(mGrid),
     &       Rho(nRho,mGrid), R2_trial(nNQ),
     &       CMOs(nCMO),P2mo(np2act),D1mo(nd1mo), Temp(mGrad),
     &       P2_ontop(nP2_ontop,mGrid), Grad(nGrad),
     &       F_xc(mGrid),dF_dRho(ndF_dRho,mGrid),
     &       dF_dP2ontop(ndF_dp2ontop,mGrid)
      Real*8 TmpPUVX(nTmpPUVX)
      Logical Check, Do_Grad, Rsv_Tsk
      Logical Do_Mo,Do_TwoEl,l_Xhol,l_casdft,Exist
      Character*4 DFTFOCK
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
      iGridInfo(i,iNQ)=iWork(ip_GridInfo+(iNQ-1)*2+i-1)
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
* Initializations for MC-PDFT                                          *
************************************************************************
      l_casdft = .false.
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
************************************************************************
* Open file for MC-PDFT to store density, pair density and ratio:      *
*                   ratio = 4pi/rho^2                                  *
************************************************************************
      IF(l_casdft) then
!
        PUVX_Time= 0d0
        FA_Time = 0d0
        sp_time = 0d0
        FI_time = 0d0
!
        LuMC=37
        LuMT=37
        call OpnFl('MCPDFT',LuMC,Exist)
        Call append_file(LuMC)
        call OpnFl('MCTRUD',LuMT,Exist)
        write(LuMC,'(A)') ' Here densities are MCPDFT modified ones'
        write(LuMC,'(A)') '       X         Y        Z'//
     &   '            d_alpha     d_beta       dTot         P2'//
     &   '         ratio'
        write(LuMT,'(A)') '     X    ,     Y    ,     Z    ,'//
     &                    '       d_a*W     ,       d_b*W     ,'//
     &                    '       dTot*W    ,       Weights   ,'//
     &                    '       dTot '
      END IF
************************************************************************
*
*----- Desymmetrize the 1-particle density matrix
*
      Call Allok2_Funi(nD)
      Call DeDe_Funi(Density,nFckInt,nD,mDens,ipDq)
*
*----- Setup for differential Rho
*
#ifdef _RDIFF_
      iDisk_Now=iDisk_Grid
      If (Grid_Status.eq.Use_Old.and..Not.Do_Grad.and.
     &    Functional_Type.eq.Old_Functional_Type) Then
*
*---------Read desymmetrized densities from previous iteration
*
         Call Allocate_Work(ipDOld,nDeDe_DFT)
         Call dDaFile(Lu_Grid,2,Work(ipDOld),nDeDe_DFT,iDisk_Now)
*
      End If
*
      Call dDaFile(Lu_Grid,1,Work(ipDeDe),nDeDe_DFT,iDisk_Grid)
*
      If (Grid_Status.eq.Use_Old.and..Not.Do_Grad.and.
     &    Functional_Type.eq.Old_Functional_Type) Then
*
*------- Form the differential desymmetrized density
*
         Call DaXpY_(nDeDe_DFT,-One,Work(ipDOld),1,Work(ipDeDe),1)
         Call Free_Work(ipDOld)
*
      End If
#endif

      If(l_casdft.and.do_pdftPot) then
        CALL GETMEM('OE_OT','ALLO','REAL',LOE_DB,nFckInt)
        CALL GETMEM('TEG_OT','ALLO','REAL',LTEG_DB,nTmpPUVX)
        Call GETMEM('FI_V','ALLO','REAL',ifiv,nFckInt)
        Call GETMEM('FI_A','ALLO','REAL',ifav,nFckInt)
!        Call GETMEM('FI_V','ALLO','REAL',ifiv_n,nFckInt)
!        Call GETMEM('FI_A','ALLO','REAL',ifav_n,nFckInt)

        CALL DCOPY_(nFckInt,[0.0D0],0,WORK(LOE_DB),1)!NTOT1
        CALL DCOPY_(nTmpPUVX,[0.0D0],0,WORK(LTEG_DB),1)
        CALL DCOPY_(nFckInt,[0.0D0],0,WORK(ifiv),1)
        CALL DCOPY_(nFckInt,[0.0D0],0,WORK(ifav),1)
!        CALL DCOPY_(nFckInt,[0.0D0],0,WORK(ifiv_n),1)
!        CALL DCOPY_(nFckInt,[0.0D0],0,WORK(ifav_n),1)
      Else
        LOE_DB = ip_dummy
        LTEG_DB = ip_dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     For a parallel implementation the iterations over
*     subblocks are parallelized.
*
      Call Init_Tsk(id,number_of_subblocks) ! Initialize parallelization
*
*-----Loop over subblocks
*
      Flop=Zero
      iSB = 0
C     Do iSB = 1, number_of_subblocks
*
*------- Start of parallelized loop here!
*
 100     Continue
         If (Grid_Status.eq.Regenerate) Then
*           Try to get an iSB to execute. If fail: done and branch out!
            If (.Not.Rsv_Tsk(id,iSB)) Go To 200
         Else
*           Try to find a subblock which was generated by this processor.
            iSB = iSB + 1
            If (iSB.gt.number_of_subblocks) Go To 200
            If (iGridInfo(2,iSB).eq.0) Go To 100
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Eliminate redundant subblocks in case of symmetry.
*        This is only done for the Lebedev grids!
*
         If (nIrrep.ne.1.and.Check(iOpt_Angular,3)) Then
*
*---------- Resolve triplet index
*
            iyz = 1 + (iSB-1)/nx
            ix  = iSB - (iyz-1)*nx
            iz  = 1 + (iyz-1)/ny
            iy  = iyz - (iz-1)*ny
*
*---------- Do symmetry by procastination.
*
            Do iIrrep = 1, nIrrep-1
               jx=ix
               If (iAnd(iOper(iIrrep),1).ne.0) jx=nx-jx+1
               jy=iy
               If (iAnd(iOper(iIrrep),2).ne.0) jy=ny-jy+1
               jz=iz
               If (iAnd(iOper(iIrrep),4).ne.0) jz=nz-jz+1
*
               jyz = (jz -1)*ny    + jy
               jxyz= (jyz-1)*nx    + jx
C              If (jxyz.gt.iSB) Go To 777
               If (jxyz.gt.iSB) Go To 100 ! go for the next task.
            End Do
*
         End If
         Debug=.False.
C        If (iSB.eq.58) Debug=.True.
C        Debug=.True.
         If (Debug) Write (6,*) 'DrvNQ_: iSB=',iSB
*                                                                      *
************************************************************************
*                                                                      *
*-----Here the List_S is the list of all the complete shells for
*     the whole system.
*
         Call Get_Subblock(Kernel,Func,iSB,
     &                     Maps2p,list_s,list_exp,list_bas,nShell,nSym,
     &                     list_p,R2_trial,nNQ,
     &                     AOInt,nAOInt,FckInt,nFckDim,nFckInt,
     &                     SOTemp,nSOTemp,
     &                     Work(ipDq),mDens,nD,
     &                     Grid,Weights,Rho,mGrid,nRho,
     &                     ndF_dRho,nP2_ontop,ndF_dP2ontop,
     &                     Do_Mo,Do_TwoEl,l_Xhol,
     &                     TmpPUVX,nTmpPUVX,nMOs,CMOs,nCMO,DoIt,
     &                     P2mo,np2act,D1mo,nd1mo,P2_ontop,
     &                     Do_Grad,Grad,nGrad,List_G,IndGrd,iTab,Temp,
cGLM     &                     mGrad,F_xc,F_xca,F_xcb,dF_dRho,dF_dP2ontop,
     &                     mGrad,F_xc,dF_dRho,dF_dP2ontop,
     &                     DFTFOCK,mAO,mdRho_dR,
     &                     LOE_DB,LTEG_DB)
*                                                                      *
************************************************************************
*                                                                      *
         Go To 100  ! go back and try to do another task
C777     Continue
C     End Do ! number_of_subblocks
 200  Continue ! Done!
      Call Free_Tsk(id)
      Flop=Flop/DBLE(nFckInt)


*                                                                      *
************************************************************************
*                                                                      *
*---- Scale result with respect to the degeneracy of the grid points
*
      If (nIrrep.ne.1.and.Check(iOpt_Angular,3)) Then
*
         Func   = DBLE(nIrrep)*Func
         Funcaa  = DBLE(nIrrep)*Funcaa
         Funcbb  = DBLE(nIrrep)*Funcbb
         Funccc  = DBLE(nIrrep)*Funccc
         Dens_I = DBLE(nIrrep)*Dens_I
         Dens_a1 = DBLE(nIrrep)*Dens_a1
         Dens_b1 = DBLE(nIrrep)*Dens_b1
         Dens_a2 = DBLE(nIrrep)*Dens_a2
         Dens_b2 = DBLE(nIrrep)*Dens_b2
         Dens_t1 = DBLE(nIrrep)*Dens_t1
         Dens_t2 = DBLE(nIrrep)*Dens_t2
         Grad_I = DBLE(nIrrep)*Grad_I
         Tau_I  = DBLE(nIrrep)*Tau_I
*
         Call DScal_(nFckInt*nFckDim,DBLE(nIrrep),FckInt,1)
         If (Do_TwoEl)
     &   Call DScal_(nTmpPUVX,DBLE(nIrrep),TmpPUVX,1)
*
      End If

*                                                                      *
************************************************************************
*                                                                      *
*---- Free memory associated with the density
*
      Call Free_DeDe_Funi(Density,nFckInt,ipDq)
*
*---- Free memory for angular grids
*
      Do iSet = 1, nAngularGrids
         Call GetMem('AngRW','Free','Real',Info_Ang(3,iSet),nDum)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Debug=.True.
      If (Debug.and..Not.Do_Grad) Then
         Write (6,*) 'Func=',Func
         iOff=1
         Do iIrrep=0,nIrrep-1
            nB=nBas(iIrrep)
            If (nB.gt.0) Then
               Call TriPrt('Final FckInt(Alpha)',' ',FckInt(iOff,1),nB)
               lB=nB*(nB+1)/2
               iOff = iOff + lB
            End If
         End Do
         If (nD.eq.1) Go To 98
         iOff=1
         Do iIrrep=0,nIrrep-1
            nB=nBas(iIrrep)
            If (nB.gt.0) Then
               Call TriPrt('Final FckInt(Beta)',' ',FckInt(iOff,2),nB)
               lB=nB*(nB+1)/2
               iOff = iOff + lB
            End If
         End Do
 98      Continue
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     For parallel implementation syncronize here!
*
*     Data to be syncronized: FckInt, Func, Dens,  and Grad.
*
*
      If (Do_Grad) Then
         Call GADSum(Grad,nGrad)
      Else
         Call GADSum_SCAL(Func)
         Call GADSum_SCAL(Funcaa)
         Call GADSum_SCAL(Funcbb)
         Call GADSum_SCAL(Funccc)
         Call GADSum_SCAL(Dens_I)
         Call GADSum_SCAL(Dens_t1)
         Call GADSum_SCAL(Dens_t2)
         Call GADSum_SCAL(Dens_a1)
         Call GADSum_SCAL(Dens_a2)
         Call GADSum_SCAL(Dens_b1)
         Call GADSum_SCAL(Dens_b2)
         Call GADSum_SCAL(Grad_I)
         Call GADSum_SCAL(Tau_I)
         Call GADSum(FckInt,nFckInt*nD)
        If(l_casdft.and.do_pdftPot) then
          Call GADSum(Work(LOE_DB),nFckInt)
          Call GADSum(Work(LTEG_DB),nTmpPUVX)
          Call GADSum(Work(ifiv),nFckInt)
          Call GADSum(Work(ifav),nFckInt)
        End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If(l_casdft.and.do_pdftPot) then


!        CALL DCOPY_(nFckInt,Work(ifav_n),1,WORK(ifav),1)
!        CALL DCOPY_(nFckInt,Work(ifiv_n),1,WORK(ifiv),1)

        Call Put_dArray('ONTOPO',work(LOE_DB),nFckInt)
        Call Put_dArray('ONTOPT',work(LTEG_DB),nTmpPUVX)
        Call Put_dArray('FI_V',Work(ifiv),nFckInt)
        Call Put_dArray('FA_V',Work(ifav),nFckInt)


        CALL GETMEM('OE_OT','Free','REAL',LOE_DB,nFckInt)
        CALL GETMEM('TEG_OT','Free','REAL',LTEG_DB,nTmpPUVX)
        CALL GETMEM('FI_V','FREE','REAL',ifiv,nFckInt)
!        CALL GETMEM('FI_V','FREE','REAL',ifiv_n,nFckInt)
        CALL GETMEM('FA_V','FREE','REAL',ifav,nFckInt)
!        CALL GETMEM('FA_V','FREE','REAL',ifav_n,nFckInt)

!      write(*,*) 'Potential timings:'
!      write(*,*) 'PUVX time: ',PUVX_Time
!      write(*,*) 'FA time: ',FA_Time
!      write(*,*) 'FI time: ',FI_Time
!      write(*,*) 'SP time: ',SP_Time
!        PUVX_Time= 0d0
!        FA_Time = 0d0
!        FI_time = 0d0
!        sp_time = 0d0
      End If

      IF(debug. and. l_casdft) THEN
        write(6,*) 'Dens_I in drvnq_ :', Dens_I
        write(6,*) 'Dens_a1 in drvnq_ :', Dens_a1
        write(6,*) 'Dens_b1 in drvnq_ :', Dens_b1
        write(6,*) 'Dens_a2 in drvnq_ :', Dens_a2
        write(6,*) 'Dens_b2 in drvnq_ :', Dens_b2
        write(6,*) 'Dens_t1 in drvnq_ :', Dens_t1
        write(6,*) 'Dens_t2 in drvnq_ :', Dens_t2
        write(6,*) 'Func in drvnq_ :', Func
        write(6,*) 'Funcaa in drvnq_ :', Funcaa
        write(6,*) 'Funcbb in drvnq_ :', Funcbb
        write(6,*) 'Funccc in drvnq_ :', Funccc
      END IF
*
* Close these files...
      If(l_casdft) then
        Close(LuMC)
        Close(LuMT)
      End if

      Return
      End
