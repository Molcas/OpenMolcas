***********************************************************************
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
      Subroutine DrvNQ_Inner(Kernel,Func,
     &                  Maps2p,nSym,list_s,list_exp,list_bas,
     &                  nShell,list_p,nNQ,
     &                  FckInt,nFckDim,Density,nFckInt,nD,
     &                  mGrid,nP2_ontop,Do_Mo,nTmpPUVX,
     &                  Do_Grad,Grad,nGrad,mAO,mdRho_dR)
************************************************************************
*                                                                      *
* Object: Driver for numerical quadrature.                             *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
#ifdef _DEBUGPRINT_
      use Basis_Info, only: nBas
#endif
      use Real_Spherical
      use Symmetry_Info, only: nIrrep, iOper
      use KSDFT_Info, only: KSDFA, LuMC, LuMT, Funcaa, Funcbb, Funccc
      use nq_Grid, only: l_casdft, D1UnZip, P2UnZip
      use nq_MO, only: D1MO, P2MO
      use Grid_On_Disk
      Implicit Real*8 (A-H,O-Z)
      External Kernel, Rsv_Tsk
#include "real.fh"
#include "WrkSpc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "nq_info.fh"
#include "debug.fh"
#include "ksdft.fh"
#include "stdalloc.fh"
      Integer Maps2p(nShell,0:nSym-1),
     &        list_s(nSym*nShell), list_exp(nSym*nShell),
     &        list_p(nNQ), list_bas(2,nSym*nShell)
      Real*8 FckInt(nFckInt,nFckDim),Density(nFckInt,nD), Grad(nGrad)
      Logical Check, Do_Grad, Rsv_Tsk
      Logical Do_Mo,Exist,l_tgga
      REAL*8,DIMENSION(:),Allocatable:: PDFTPot1,PDFTFocI,PDFTFocA
      Real*8, Allocatable:: OE_OT(:), EG_OT(:)
      Real*8, Allocatable:: FI_V(:), FA_V(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
* Initializations for MC-PDFT                                          *
************************************************************************
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
       IF(Debug) THEN
        LuMC=37
        call OpnFl('MCPDFT',LuMC,Exist)
c        Call append_file(LuMC)
        write(LuMC,'(A)') ' Here densities are MCPDFT modified ones.'
        write(LuMC,*)     ' Used by translated functional: ', KSDFA(1:8)
        write(LuMC,'(A)') '     X    ,     Y    ,     Z    ,'//
     &                    '       d_a*W     ,       d_b*W     ,'//
     &                    '       dTot*W    ,       Weights   ,'//
     &                    '          dTot   ,       P2        ,   ratio'
        LuMT=37
        call OpnFl('MCTRUD',LuMT,Exist)
c        Call append_file(LuMT)
        write(LuMT,'(A)') ' Here densities are original ones.'
        write(LuMT,*)     ' Used by translated functional: ', KSDFA(1:8)
        write(LuMT,'(A)') '     X    ,     Y    ,     Z    ,'//
     &                    '       d_a*W     ,       d_b*W     ,'//
     &                    '       dTot*W    ,       Weights   ,'//
     &                    '       dTot '
       END IF

      CALL CalcOrbOff()
      NASHT4=NASHT**4
      CALL mma_allocate(P2Unzip,NASHT4)
      CALL mma_allocate(D1Unzip,NASHT**2)
      CALL UnzipD1(D1Unzip,D1MO,SIZE(D1MO))
      CALL UnzipP2(P2Unzip,P2MO,SIZE(P2MO))
      END IF
************************************************************************
*
*----- Desymmetrize the 1-particle density matrix
*
      Call Allok2_Funi(nD)
      Call DeDe_Funi(Density,nFckInt,nD)
*
      If(l_casdft.and.do_pdftPot) then
        CALL mma_allocate(PDFTPot1,nPot1)
        CALL mma_allocate(PDFTFocI,nPot1)
        CALL mma_allocate(PDFTFocA,nPot1)
        CALL mma_allocate(OE_OT,nFckInt,Label='OE_OT')
        CALL mma_allocate(EG_OT,nTmpPUVX,Label='EG_OT')
        Call mma_allocate(FI_V,nFckInt,Label='FI_V')
        Call mma_allocate(FA_V,nFckInt,Label='FA_V')

        OE_OT(:)=Zero
        EG_OT(:)=Zero
        FI_V(:)=Zero
        FA_V(:)=Zero
        CALL FZero(PDFTPot1,nPot1)
        CALL FZero(PDFTFocI,nPot1)
        CALL FZero(PDFTFocA,nPot1)
        CALL CalcPUVXOff()
      Else
        nPot1=1
        CALL mma_allocate(OE_OT,nPot1,Label='OE_OT')
        CALL mma_allocate(EG_OT,nPot1,Label='EG_OT')
        Call mma_allocate(FI_V,nPot1,Label='FI_V')
        Call mma_allocate(FA_V,nPot1,Label='FA_V')
        CALL mma_allocate(PDFTPot1,nPot1)
        CALL mma_allocate(PDFTFocI,nPot1)
        CALL mma_allocate(PDFTFocA,nPot1)
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
            If (GridInfo(2,iSB).eq.0) Go To 100
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
     &                     list_p,nNQ,FckInt,nFckDim,nFckInt,nD,
     &                     mGrid,nP2_ontop,Do_Mo,
     &                     Do_Grad,Grad,nGrad,
     &                     mAO,mdRho_dR,
     &                     EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
*                                                                      *
************************************************************************
*                                                                      *
         Go To 100  ! go back and try to do another task
C777     Continue
C     End Do ! number_of_subblocks
 200  Continue ! Done!
      Call Free_Tsk(id)


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
*
      End If

*                                                                      *
************************************************************************
*                                                                      *
*---- Free memory associated with the density
*
      Call Free_DeDe_Funi()
*
*---- Free memory for angular grids
*
      Do iSet = 1, nAngularGrids
         Call GetMem('AngRW','Free','Real',Info_Ang(3,iSet),nDum)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
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

      IF(l_casdft) THEN
        CALL mma_deallocate(D1Unzip)
        CALL mma_deallocate(P2Unzip)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
*     For parallel implementation syncronize here!
*
*     Data to be syncronized: FckInt, Func, Dens,  and Grad.
*
*
        l_tgga = .true.
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
          Call GADSum(OE_OT,nFckInt)
          Call GADSum(EG_OT,nTmpPUVX)
          Call GADSum(FI_V,nFckInt)
          Call GADSum(FA_V,nFckInt)
          if(l_tgga) then
           CALL GADSum(PDFTPot1,nPot1)
           CALL GADSum(PDFTFocI,nPot1)
           CALL GADSum(PDFTFocA,nPot1)
          end if
        End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If(l_casdft.and.do_pdftPot) then


        If(l_tgga) Then
         CALL PackPot1(OE_OT,PDFTPot1,nFckInt,dble(nIrrep)*0.5d0)
         CALL DScal_(nPot2,dble(nIrrep),EG_OT,1)
         CALL PackPot1(FI_V,PDFTFocI,nFckInt,dble(nIrrep)*0.25d0)
         CALL PackPot1(FA_V,PDFTFocA,nFckInt,dble(nIrrep)*0.5d0)
        End If
        Call Put_dArray('ONTOPO',OE_OT,nFckInt)
        Call Put_dArray('ONTOPT',EG_OT,nTmpPUVX)
        Call Put_dArray('FI_V',FI_V,nFckInt)
        Call Put_dArray('FA_V',FA_V,nFckInt)

      End If
      Call mma_deallocate(OE_OT)
      Call mma_deallocate(EG_OT)
      Call mma_deallocate(FA_V)
      Call mma_deallocate(FI_V)
      CALL mma_deallocate(PDFTPot1)
      CALL mma_deallocate(PDFTFocI)
      CALL mma_deallocate(PDFTFocA)

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
*
* Close these files...
        Close(LuMC)
        Close(LuMT)
      END IF

      Return
      End
