!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************

subroutine DrvNQ_Inner(Kernel,Func,Maps2p,nSym,list_s,list_exp,list_bas,nShell,list_p,nNQ,FckInt,nFckDim,Density,nFckInt,nD,mGrid, &
                       nP2_ontop,Do_Mo,nTmpPUVX,Do_Grad,Grad,nGrad,mAO,mdRho_dR)
!***********************************************************************
!                                                                      *
! Object: Driver for numerical quadrature.                             *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept of Chemical Physics,                                *
!             University of Lund, Sweden                               *
!             August 1999                                              *
!***********************************************************************

#ifdef _DEBUGPRINT_
use Basis_Info, only: nBas
#endif
use Real_Spherical
use Symmetry_Info, only: nIrrep, iOper
use KSDFT_Info, only: do_pdftpot, FA_time, FI_time, Funcaa, Funcbb, Funccc, KSDFA, LuMC, LuMT, PUVX_time, sp_time
use nq_Grid, only: l_casdft, D1UnZip, P2UnZip
use nq_MO, only: D1MO, P2MO
use nq_Structure, only: Close_Info_Ang
use Grid_On_Disk
use nq_Info

implicit real*8(A-H,O-Z)
external Kernel, Rsv_Tsk
#include "real.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "debug.fh"
#include "stdalloc.fh"
integer Maps2p(nShell,0:nSym-1), list_s(nSym*nShell), list_exp(nSym*nShell), list_p(nNQ), list_bas(2,nSym*nShell)
real*8 FckInt(nFckInt,nFckDim), Density(nFckInt,nD), Grad(nGrad)
logical Check, Do_Grad, Rsv_Tsk
logical Do_Mo, Exist, l_tgga
real*8, dimension(:), allocatable :: PDFTPot1, PDFTFocI, PDFTFocA
real*8, allocatable :: OE_OT(:), EG_OT(:)
real*8, allocatable :: FI_V(:), FA_V(:)
! Statement functions
Check(i,j) = iand(i,2**(j-1)) /= 0

!***********************************************************************
! Initializations for MC-PDFT                                          *
!***********************************************************************
!***********************************************************************
! Open file for MC-PDFT to store density, pair density and ratio:      *
!                   ratio = 4pi/rho^2                                  *
!***********************************************************************
if (l_casdft) then

  PUVX_Time = 0d0
  FA_Time = 0d0
  sp_time = 0d0
  FI_time = 0d0

  if (Debug) then
    LuMC = 37
    call OpnFl('MCPDFT',LuMC,Exist)
    !call append_file(LuMC)
    write(LuMC,'(A)') ' Here densities are MCPDFT modified ones.'
    write(LuMC,*) ' Used by translated functional: ',KSDFA(1:8)
    write(LuMC,'(A)') '     X    ,     Y    ,     Z    ,       d_a*W     ,       d_b*W     ,       dTot*W    ,'// &
                      '       Weights   ,          dTot   ,       P2        ,   ratio'
    LuMT = 37
    call OpnFl('MCTRUD',LuMT,Exist)
    !call append_file(LuMT)
    write(LuMT,'(A)') ' Here densities are original ones.'
    write(LuMT,*) ' Used by translated functional: ',KSDFA(1:8)
    write(LuMT,'(A)') '     X    ,     Y    ,     Z    ,       d_a*W     ,       d_b*W     ,       dTot*W    ,'// &
                      '       Weights   ,       dTot '
  end if

  call CalcOrbOff()
  NASHT4 = NASHT**4
  call mma_allocate(P2Unzip,NASHT4)
  call mma_allocate(D1Unzip,NASHT**2)
  call UnzipD1(D1Unzip,D1MO,size(D1MO))
  call UnzipP2(P2Unzip,P2MO,size(P2MO))
end if
!***********************************************************************
!
! Desymmetrize the 1-particle density matrix

call Allok2_Funi(nD)
call DeDe_Funi(Density,nFckInt,nD)

if (l_casdft .and. do_pdftPot) then
  call mma_allocate(PDFTPot1,nPot1)
  call mma_allocate(PDFTFocI,nPot1)
  call mma_allocate(PDFTFocA,nPot1)
  call mma_allocate(OE_OT,nFckInt,Label='OE_OT')
  call mma_allocate(EG_OT,nTmpPUVX,Label='EG_OT')
  call mma_allocate(FI_V,nFckInt,Label='FI_V')
  call mma_allocate(FA_V,nFckInt,Label='FA_V')

  OE_OT(:) = Zero
  EG_OT(:) = Zero
  FI_V(:) = Zero
  FA_V(:) = Zero
  call FZero(PDFTPot1,nPot1)
  call FZero(PDFTFocI,nPot1)
  call FZero(PDFTFocA,nPot1)
  call CalcPUVXOff()
else
  nPot1 = 1
  call mma_allocate(OE_OT,nPot1,Label='OE_OT')
  call mma_allocate(EG_OT,nPot1,Label='EG_OT')
  call mma_allocate(FI_V,nPot1,Label='FI_V')
  call mma_allocate(FA_V,nPot1,Label='FA_V')
  call mma_allocate(PDFTPot1,nPot1)
  call mma_allocate(PDFTFocI,nPot1)
  call mma_allocate(PDFTFocA,nPot1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For a parallel implementation the iterations over
! subblocks are parallelized.

call Init_Tsk(id,number_of_subblocks) ! Initialize parallelization

! Loop over subblocks

iSB = 0
!do iSB=1,number_of_subblocks

! Start of parallelized loop here!

100 continue
if (Grid_Status == Regenerate) then
  ! Try to get an iSB to execute. If fail: done and branch out!
  if (.not. Rsv_Tsk(id,iSB)) Go To 200
else
  ! Try to find a subblock which was generated by this processor.
  iSB = iSB+1
  if (iSB > number_of_subblocks) Go To 200
  if (GridInfo(2,iSB) == 0) Go To 100
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Eliminate redundant subblocks in case of symmetry.
! This is only done for the Lebedev grids!

if ((nIrrep /= 1) .and. Check(iOpt_Angular,3)) then

  ! Resolve triplet index

  iyz = 1+(iSB-1)/nx
  ix = iSB-(iyz-1)*nx
  iz = 1+(iyz-1)/ny
  iy = iyz-(iz-1)*ny

  ! Do symmetry by procastination.

  do iIrrep=1,nIrrep-1
    jx = ix
    if (iand(iOper(iIrrep),1) /= 0) jx = nx-jx+1
    jy = iy
    if (iand(iOper(iIrrep),2) /= 0) jy = ny-jy+1
    jz = iz
    if (iand(iOper(iIrrep),4) /= 0) jz = nz-jz+1

    jyz = (jz-1)*ny+jy
    jxyz = (jyz-1)*nx+jx
    !if (jxyz > iSB) Go To 777
    if (jxyz > iSB) Go To 100 ! go for the next task.
  end do

end if
Debug = .false.
!if (iSB == 58) Debug = .true.
!Debug = .true.
if (Debug) write(6,*) 'DrvNQ_: iSB=',iSB
!                                                                      *
!***********************************************************************
!                                                                      *
! Here the List_S is the list of all the complete shells for
! the whole system.

call Get_Subblock(Kernel,Func,iSB,Maps2p,list_s,list_exp,list_bas,nShell,nSym,list_p,nNQ,FckInt,nFckDim,nFckInt,nD,mGrid, &
                  nP2_ontop,Do_Mo,Do_Grad,Grad,nGrad,mAO,mdRho_dR,EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
!                                                                      *
!***********************************************************************
!                                                                      *
Go To 100 ! go back and try to do another task
!777 continue
!end do ! number_of_subblocks
200 continue ! Done!
call Free_Tsk(id)

!                                                                      *
!***********************************************************************
!                                                                      *
! Scale result with respect to the degeneracy of the grid points

if ((nIrrep /= 1) .and. Check(iOpt_Angular,3)) then

  Func = dble(nIrrep)*Func
  Funcaa = dble(nIrrep)*Funcaa
  Funcbb = dble(nIrrep)*Funcbb
  Funccc = dble(nIrrep)*Funccc
  Dens_I = dble(nIrrep)*Dens_I
  Dens_a1 = dble(nIrrep)*Dens_a1
  Dens_b1 = dble(nIrrep)*Dens_b1
  Dens_a2 = dble(nIrrep)*Dens_a2
  Dens_b2 = dble(nIrrep)*Dens_b2
  Dens_t1 = dble(nIrrep)*Dens_t1
  Dens_t2 = dble(nIrrep)*Dens_t2
  Grad_I = dble(nIrrep)*Grad_I
  Tau_I = dble(nIrrep)*Tau_I

  call DScal_(nFckInt*nFckDim,dble(nIrrep),FckInt,1)

end if

!                                                                      *
!***********************************************************************
!                                                                      *
! Free memory associated with the density

call Free_DeDe_Funi()

! Free memory for angular grids

call Close_Info_Ang()
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Debug = .true.
if (Debug .and. (.not. Do_Grad)) then
  write(6,*) 'Func=',Func
  iOff = 1
  do iIrrep=0,nIrrep-1
    nB = nBas(iIrrep)
    if (nB > 0) then
      call TriPrt('Final FckInt(Alpha)',' ',FckInt(iOff,1),nB)
      lB = nB*(nB+1)/2
      iOff = iOff+lB
    end if
  end do
  if (nD == 1) Go To 98
  iOff = 1
  do iIrrep=0,nIrrep-1
    nB = nBas(iIrrep)
    if (nB > 0) then
      call TriPrt('Final FckInt(Beta)',' ',FckInt(iOff,2),nB)
      lB = nB*(nB+1)/2
      iOff = iOff+lB
    end if
  end do
98 continue
end if
#endif

if (l_casdft) then
  call mma_deallocate(D1Unzip)
  call mma_deallocate(P2Unzip)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For parallel implementation synchronize here!
!
! Data to be synchronized: FckInt, Func, Dens,  and Grad.

l_tgga = .true.
if (Do_Grad) then
  call GADSum(Grad,nGrad)
else
  call GADSum_SCAL(Func)
  call GADSum_SCAL(Funcaa)
  call GADSum_SCAL(Funcbb)
  call GADSum_SCAL(Funccc)
  call GADSum_SCAL(Dens_I)
  call GADSum_SCAL(Dens_t1)
  call GADSum_SCAL(Dens_t2)
  call GADSum_SCAL(Dens_a1)
  call GADSum_SCAL(Dens_a2)
  call GADSum_SCAL(Dens_b1)
  call GADSum_SCAL(Dens_b2)
  call GADSum_SCAL(Grad_I)
  call GADSum_SCAL(Tau_I)
  call GADSum(FckInt,nFckInt*nD)
  if (l_casdft .and. do_pdftPot) then
    call GADSum(OE_OT,nFckInt)
    call GADSum(EG_OT,nTmpPUVX)
    call GADSum(FI_V,nFckInt)
    call GADSum(FA_V,nFckInt)
    if (l_tgga) then
      call GADSum(PDFTPot1,nPot1)
      call GADSum(PDFTFocI,nPot1)
      call GADSum(PDFTFocA,nPot1)
    end if
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (l_casdft .and. do_pdftPot) then

  if (l_tgga) then
    call PackPot1(OE_OT,PDFTPot1,nFckInt,dble(nIrrep)*0.5d0)
    call DScal_(nPot2,dble(nIrrep),EG_OT,1)
    call PackPot1(FI_V,PDFTFocI,nFckInt,dble(nIrrep)*0.25d0)
    call PackPot1(FA_V,PDFTFocA,nFckInt,dble(nIrrep)*0.5d0)
  end if
  call Put_dArray('ONTOPO',OE_OT,nFckInt)
  call Put_dArray('ONTOPT',EG_OT,nTmpPUVX)
  call Put_dArray('FI_V',FI_V,nFckInt)
  call Put_dArray('FA_V',FA_V,nFckInt)

end if
call mma_deallocate(OE_OT)
call mma_deallocate(EG_OT)
call mma_deallocate(FA_V)
call mma_deallocate(FI_V)
call mma_deallocate(PDFTPot1)
call mma_deallocate(PDFTFocI)
call mma_deallocate(PDFTFocA)

if (debug .and. l_casdft) then
  write(6,*) 'Dens_I in drvnq_ :',Dens_I
  write(6,*) 'Dens_a1 in drvnq_ :',Dens_a1
  write(6,*) 'Dens_b1 in drvnq_ :',Dens_b1
  write(6,*) 'Dens_a2 in drvnq_ :',Dens_a2
  write(6,*) 'Dens_b2 in drvnq_ :',Dens_b2
  write(6,*) 'Dens_t1 in drvnq_ :',Dens_t1
  write(6,*) 'Dens_t2 in drvnq_ :',Dens_t2
  write(6,*) 'Func in drvnq_ :',Func
  write(6,*) 'Funcaa in drvnq_ :',Funcaa
  write(6,*) 'Funcbb in drvnq_ :',Funcbb
  write(6,*) 'Funccc in drvnq_ :',Funccc

  ! Close these files...
  close(LuMC)
  close(LuMT)
end if

return

end subroutine DrvNQ_Inner
