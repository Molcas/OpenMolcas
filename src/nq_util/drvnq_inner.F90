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

use Symmetry_Info, only: nIrrep, iOper
use KSDFT_Info, only: do_pdftpot, FA_time, FI_time, Funcaa, Funcbb, Funccc, PUVX_time, sp_time
use nq_Grid, only: l_casdft, D1UnZip, P2UnZip
use nq_MO, only: D1MO, P2MO
use nq_Structure, only: Close_Info_Ang
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Grad_I, iOpt_Angular, NASHT, nPot1, nPot2, &
                   number_of_subblocks, nx, ny, nz, Tau_I
use Grid_On_Disk, only: Grid_Status, GridInfo, Regenerate
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, Quart
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Basis_Info, only: nBas
use KSDFT_Info, only: KSDFA, LuMC, LuMT
use Definitions, only: u6
#endif

implicit none
external :: Kernel
integer(kind=iwp), intent(in) :: nShell, nSym, Maps2p(nShell,0:nSym-1), nNQ, nFckDim, nFckInt, nD, mGrid, nP2_ontop, nTmpPUVX, &
                                 nGrad, mAO, mdRho_dR
integer(kind=iwp), intent(out) :: list_s(nSym*nShell), list_exp(nSym*nShell), list_bas(2,nSym*nShell), list_p(nNQ)
real(kind=wp), intent(inout) :: Func, FckInt(nFckInt,nFckDim), Grad(nGrad)
real(kind=wp), intent(in) :: Density(nFckInt,nD)
logical(kind=iwp), intent(in) :: Do_Mo, Do_Grad
integer(kind=iwp) :: id, iIrrep, iSB, ix, iy, iyz, iz, jx, jxyz, jy, jyz, jz
logical(kind=iwp) :: l_tgga
real(kind=wp), allocatable :: EG_OT(:), FA_V(:), FI_V(:), OE_OT(:), PDFTFocA(:), PDFTFocI(:), PDFTPot1(:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iOff, lB, nB
logical(kind=iwp) :: Exists
#endif
logical(kind=iwp), external :: Rsv_Tsk

!***********************************************************************
! Initializations for MC-PDFT                                          *
!***********************************************************************
!***********************************************************************
! Open file for MC-PDFT to store density, pair density and ratio:      *
!                   ratio = 4pi/rho^2                                  *
!***********************************************************************
if (l_casdft) then

  PUVX_Time = Zero
  FA_Time = Zero
  sp_time = Zero
  FI_time = Zero

# ifdef _DEBUGPRINT_
  LuMC = 37
  call OpnFl('MCPDFT',LuMC,Exists)
  !call append_file(LuMC)
  write(LuMC,'(A)') ' Here densities are MCPDFT modified ones.'
  write(LuMC,*) ' Used by translated functional: ',KSDFA(1:8)
  write(LuMC,'(A)') '     X    ,     Y    ,     Z    ,       d_a*W     ,       d_b*W     ,       dTot*W    ,       Weights   ,'// &
                    '          dTot   ,       P2        ,   ratio'
  LuMT = 37
  call OpnFl('MCTRUD',LuMT,Exists)
  !call append_file(LuMT)
  write(LuMT,'(A)') ' Here densities are original ones.'
  write(LuMT,*) ' Used by translated functional: ',KSDFA(1:8)
  write(LuMT,'(A)') '     X    ,     Y    ,     Z    ,       d_a*W     ,       d_b*W     ,       dTot*W    ,       Weights   ,'// &
                    '          dTot '
# endif

  call CalcOrbOff()
  call mma_allocate(P2Unzip,NASHT,NASHT,NASHT,NASHT)
  call mma_allocate(D1Unzip,NASHT,NASHT)
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
  PDFTPot1(:) = Zero
  PDFTFocI(:) = Zero
  PDFTFocA(:) = Zero
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

outer: do
  if (Grid_Status == Regenerate) then
    ! Try to get an iSB to execute. If fail: done and branch out!
    if (.not. Rsv_Tsk(id,iSB)) exit outer
  else
    ! Try to find a subblock which was generated by this processor.
    iSB = iSB+1
    if (iSB > number_of_subblocks) exit outer
    if (GridInfo(2,iSB) == 0) cycle outer
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Eliminate redundant subblocks in case of symmetry.
  ! This is only done for the Lebedev grids!

  if ((nIrrep /= 1) .and. btest(iOpt_Angular,2)) then

    ! Resolve triplet index

    iyz = 1+(iSB-1)/nx
    ix = iSB-(iyz-1)*nx
    iz = 1+(iyz-1)/ny
    iy = iyz-(iz-1)*ny

    ! Do symmetry by procrastination.

    do iIrrep=1,nIrrep-1
      jx = ix
      if (btest(iOper(iIrrep),0)) jx = nx-jx+1
      jy = iy
      if (btest(iOper(iIrrep),1)) jy = ny-jy+1
      jz = iz
      if (btest(iOper(iIrrep),2)) jz = nz-jz+1

      jyz = (jz-1)*ny+jy
      jxyz = (jyz-1)*nx+jx
      !if (jxyz > iSB) exit outer
      if (jxyz > iSB) cycle outer ! go for the next task.
    end do

  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'DrvNQ_Inner: iSB=',iSB
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Here the List_S is the list of all the complete shells for
  ! the whole system.

  call Get_Subblock(Kernel,Func,iSB,Maps2p,list_s,list_exp,list_bas,nShell,nSym,list_p,nNQ,FckInt,nFckDim,nFckInt,nD,mGrid, &
                    nP2_ontop,Do_Mo,Do_Grad,Grad,nGrad,mAO,mdRho_dR,EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
  ! go back and try to do another task
end do outer
!end do ! number_of_subblocks
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Tsk(id)

!                                                                      *
!***********************************************************************
!                                                                      *
! Scale result with respect to the degeneracy of the grid points

if ((nIrrep /= 1) .and. btest(iOpt_Angular,2)) then

  Func = real(nIrrep,kind=wp)*Func
  Funcaa = real(nIrrep,kind=wp)*Funcaa
  Funcbb = real(nIrrep,kind=wp)*Funcbb
  Funccc = real(nIrrep,kind=wp)*Funccc
  Dens_I = real(nIrrep,kind=wp)*Dens_I
  Dens_a1 = real(nIrrep,kind=wp)*Dens_a1
  Dens_b1 = real(nIrrep,kind=wp)*Dens_b1
  Dens_a2 = real(nIrrep,kind=wp)*Dens_a2
  Dens_b2 = real(nIrrep,kind=wp)*Dens_b2
  Dens_t1 = real(nIrrep,kind=wp)*Dens_t1
  Dens_t2 = real(nIrrep,kind=wp)*Dens_t2
  Grad_I = real(nIrrep,kind=wp)*Grad_I
  Tau_I = real(nIrrep,kind=wp)*Tau_I

  FckInt(:,:) = real(nIrrep,kind=wp)*FckInt(:,:)

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
#ifdef _DEBUGPRINT_
if (.not. Do_Grad) then
  write(u6,*) 'Func=',Func
  iOff = 1
  do iIrrep=0,nIrrep-1
    nB = nBas(iIrrep)
    if (nB > 0) then
      call TriPrt('Final FckInt(Alpha)',' ',FckInt(iOff,1),nB)
      lB = nB*(nB+1)/2
      iOff = iOff+lB
    end if
  end do
  if (nD /= 1) then
    iOff = 1
    do iIrrep=0,nIrrep-1
      nB = nBas(iIrrep)
      if (nB > 0) then
        call TriPrt('Final FckInt(Beta)',' ',FckInt(iOff,2),nB)
        lB = nB*(nB+1)/2
        iOff = iOff+lB
      end if
    end do
  end if
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
    call PackPot1(OE_OT,PDFTPot1,nFckInt,real(nIrrep,kind=wp)*Half)
    EG_OT(1:nPot2) = real(nIrrep,kind=wp)*EG_OT(1:nPot2)
    call PackPot1(FI_V,PDFTFocI,nFckInt,real(nIrrep,kind=wp)*Quart)
    call PackPot1(FA_V,PDFTFocA,nFckInt,real(nIrrep,kind=wp)*Half)
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

#ifdef _DEBUGPRINT_
if (l_casdft) then
  write(u6,*) 'Dens_I in DrvNQ_Inner :',Dens_I
  write(u6,*) 'Dens_a1 in DrvNQ_Inner :',Dens_a1
  write(u6,*) 'Dens_b1 in DrvNQ_Inner :',Dens_b1
  write(u6,*) 'Dens_a2 in DrvNQ_Inner :',Dens_a2
  write(u6,*) 'Dens_b2 in DrvNQ_Inner :',Dens_b2
  write(u6,*) 'Dens_t1 in DrvNQ_Inner :',Dens_t1
  write(u6,*) 'Dens_t2 in DrvNQ_Inner :',Dens_t2
  write(u6,*) 'Func in DrvNQ_Inner :',Func
  write(u6,*) 'Funcaa in DrvNQ_Inner :',Funcaa
  write(u6,*) 'Funcbb in DrvNQ_Inner :',Funcbb
  write(u6,*) 'Funccc in DrvNQ_Inner :',Funccc

  ! Close these files...
  close(LuMC)
  close(LuMT)
end if
#endif

return

end subroutine DrvNQ_Inner
