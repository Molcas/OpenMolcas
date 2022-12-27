!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001, Roland Lindh                                     *
!***********************************************************************

subroutine DrvNQ(Kernel,FckInt,nFckDim,Funct,Density,nFckInt,nD,Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK,IsFT)
!***********************************************************************
!                                                                      *
!     Predriver for numerical integration utility.                     *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept of Chemical Physics,                                *
!             University of Lund, Sweden                               *
!             December 2001                                            *
!***********************************************************************

use Symmetry_Info, only: nIrrep
use nq_Grid, only: Angular, Coor, F_xc, F_xca, F_xcb, Fact, GradRho, Grid, IndGrd, iTab, kAO, l_CASDFT, Lapl, List_G, Mem, &
                   nGridMax, nR_Eff, nRho, Pax, R2_trial, Rho, Sigma, Tau, Temp, vLapl, vRho, vSigma, vTau, Weights
use nq_pdft, only: lft, lGGA
use nq_MO, only: nMOs, CMO, D1MO, P2MO, P2_ontop
use nq_Structure, only: Close_NQ_Data
use nq_Info, only: Functional_type, GGA_type, LDA_type, LMax_NQ, mBas, meta_GGA_type1, meta_GGA_type2, mIrrep, nAsh, nAtoms, nFro, &
                   number_of_subblocks, Other_type
use Grid_On_Disk, only: Final_Grid, G_S, Grid_Status, GridInfo, iDisk_Grid, iDisk_Set, iGrid_Set, Intermediate, Lu_Grid, &
                        LuGridFile, Old_Functional_Type, Regenerate, Use_Old
use libxc, only: dfunc_dLapl, dfunc_drho, dfunc_dsigma, dfunc_dTau, func
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
external :: Kernel
integer(kind=iwp), intent(in) :: nFckDim, nFckInt, nD, nGrad
real(kind=wp), intent(inout) :: FckInt(nFckInt,nFckDim), Funct, Grad(nGrad)
real(kind=wp), intent(in) :: Density(nFckInt,nD)
logical(kind=iwp), intent(in) :: Do_Grad, Do_TwoEl, IsFT
logical(kind=iwp), intent(inout) :: Do_MO
character(len=4), intent(in) :: DFTFOCK
#include "status.fh"
integer(kind=iwp) :: i, iDum(1), iIrrep, ijIrrep, ijkIrrep, iOrb, iStack, jAsh, jIrrep, kAsh, kIrrep, kl_Orb_pairs, lAsh, mAO, &
                     mdRho_dr, mGrad, nBas(8), nCMO, nD1MO, nDel(8), nGradRho, nLapl, nNQ, nP2, nP2_ontop, nSigma, nTau, NQNAC, &
                     NQNACPAR, NQNACPR2, nShell, nTmpPUVX
real(kind=wp) :: PThr
logical(kind=iwp) :: PMode
integer(kind=iwp), allocatable :: List_Bas(:,:), List_Exp(:), List_P(:), List_s(:,:), Maps2p(:,:)
real(kind=wp), allocatable :: R_Min(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
lft = IsFT
if (Do_TwoEl) Do_MO = .true.
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate enough memory for Maps2p

call Set_Basis_Mode('Valence')
call Nr_Shells(nShell)
call mma_allocate(Maps2p,nShell,nIrrep,Label='Maps2p')
call mma_allocate(R_Min,LMax_NQ+1,Label='R_Min')

NQ_Status = Inactive
call Setup_NQ(Maps2p,nShell,nIrrep,nNQ,Do_Grad,Do_MO,PThr,PMode,R_Min,LMax_NQ)

call mma_deallocate(R_Min)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory sufficiently large to store all grid points
! and associated weights.

call mma_Allocate(Grid,3,nGridMax,Label='Grid')
call mma_Allocate(Weights,nGridMax,Label='Weights')
!                                                                      *
!***********************************************************************
!                                                                      *
! CASDFT stuff:
!
! Note, l_CASDFT=.True. implies that both Do_MO and Do_Twoel are true.

nTmpPUVX = 1

NQNAC = 0
if (l_casdft) call Get_iArray('nAsh',nAsh,mIrrep)
if ((DFTFOCK /= 'SCF ') .or. l_CASDFT) then
  do iIrrep=0,mIrrep-1
    NQNAC = NQNAC+nAsh(iIrrep)
  end do
end if
NQNACPAR = (NQNAC**2+NQNAC)/2
NQNACPR2 = (NQNACPAR**2+NQNACPAR)/2

LuGridFile = 31
LuGridFile = IsFreeUnit(LuGridFile)
call Molcas_Open(LuGridFile,'GRIDFILE')

#ifdef _DEBUGPRINT_
write(u6,*) 'l_casdft value at drvnq:',l_casdft
if (l_casdft) write(u6,*) 'MCPDFT with functional:',KSDFA
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     Definition of resources needed for the functionals.              *
!                                                                      *
!     mAO: the number of derivatives needed of an basis function.      *
!          Depending of the functional type and if gradients will be   *
!          computed. Numbers will be 1, 4, 10, 20, 35, etc.            *
!     nRho:the number of parameters of the functional. Note that this  *
!          is different for the same functional depending on if it is  *
!          a closed or open-shell case.                                *
!     mdRho_dR: number of derivatives of the parameters with respect   *
!          to the nuclear coordinates. The true number is of course    *
!          three (x,y,z) times this.                                   *
!     nF_drho: the number of derivatives of the functional wrt the     *
!          parameters. Note that grad rho is not a direct parameter    *
!          but that we use gamma.                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
select case (Functional_type)

  case (LDA_type)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! We need the AOs, for gradients we need the derivatives too.

    mAO = 1
    kAO = mAO
    if (Do_Grad) mAO = 4

    ! We need rho.
    ! For gradients we need derivatives of rho wrt the coordinates

    nRho = nD
    mdRho_dr = 0
    if (Do_Grad) mdRho_dr = nD
    nSigma = 0
    nGradRho = 0
    nLapl = 0
    nTau = 0

    ! We need derivatives of the functional with respect to
    ! rho(alpha). In case of open-shell calculations we also
    ! need rho(beta).

    nP2_ontop = 1
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (GGA_type)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! We need the AOs and their derivatives, for  gradients we need
    ! the second derivatives too.

    mAO = 4
    kAO = mAO
    if (Do_Grad) mAO = 10

    ! We need rho and grad rho
    ! For gradients we need the derrivatives wrt the coordinates

    nRho = nD
    nSigma = nD*(nD+1)/2
    nGradRho = nD*3
    nLapl = 0
    nTau = 0
    mdRho_dR = 0
    if (Do_Grad) mdRho_dR = 4*nD

    ! We need derivatives of the functional with respect to
    ! rho(alpha), gamma(alpha,alpha) and gamma(alpha,beta).
    ! In case of open-shell calculations we also
    ! need rho(beta) and gamma(beta,beta).

    nP2_ontop = 4
    lGGA = .true.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! We need the AOs and their derivatives, for  gradients we need
    ! the second derivatives too.

    mAO = 4
    kAO = mAO
    if (Do_Grad) mAO = 10

    ! We need rho, grad rho and tau.
    ! For gradients we need the derrivatives wrt the coordinates

    nRho = nD
    nSigma = nD*(nD+1)/2
    nGradRho = nD*3
    !nLapl = 0
    nLapl = nD
    nTau = nD
    mdRho_dR = 0
    if (Do_Grad) mdRho_dR = 5*nD

    ! We need derivatives of the functional with respect to
    ! rho(alpha), gamma(alpha,alpha), gamma(alpha,beta) and
    ! tau(alpha). In case of open-shell calculations we also
    ! need rho(beta), gamma(beta,beta) and tau(beta).

    nP2_ontop = 4
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type2)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! We need the AOs and their 1st and 2nd derivatives, for
    ! gradients we need the 3rd order derivatives too.

    mAO = 10
    kAO = mAO
    if (Do_Grad) mAO = 20

    ! We need rho, grad rho, tau, and the Laplacian
    ! For gradients we need the derrivatives wrt the coordinates

    nRho = nD
    nSigma = nD*(nD+1)/2
    nGradRho = nD*3
    nLapl = nD
    nTau = nD
    mdRho_dR = 0
    if (Do_Grad) mdRho_dR = 6*nD

    ! We need derivatives of the functional with respect to
    ! rho(alpha), gamma(alpha,alpha), gamma(alpha,beta),
    ! tau(alpha) and laplacian(alpha). In case of open-shell
    ! calculations we also need rho(beta), gamma(beta,beta),
    ! tau(beta) and laplacian(beta).

    nP2_ontop = 4
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case default
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    Functional_type = Other_type
    call WarningMessage(2,'DrvNQ: Invalid Functional_type!')
    call Abend()
    nRho = 0
    nSigma = 0
    nGradRho = 0
    nLapl = 0
    nTau = 0
    !                                                                  *
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Rho,nRho,nGridMax,Label='Rho')
call mma_allocate(vRho,nRho,nGridMax,Label='vRho')
call mma_allocate(dfunc_drho,nRho,nGridMax,Label='dfunc_drho')
if (nSigma /= 0) then
  call mma_Allocate(Sigma,nSigma,nGridMax,Label='Sigma')
  call mma_Allocate(vSigma,nSigma,nGridMax,Label='vSigma')
  call mma_Allocate(dfunc_dSigma,nSigma,nGridMax,Label='dfunc_dSigma')
end if
if (nGradRho /= 0) then
  call mma_Allocate(GradRho,nGradRho,nGridMax,Label='GradRho')
end if
if (nTau /= 0) then
  call mma_allocate(Tau,nTau,nGridMax,Label='Tau')
  call mma_allocate(vTau,nTau,nGridMax,Label='vTau')
  call mma_allocate(dfunc_dTau,nTau,nGridMax,Label='dfunc_dTau')
  Tau(:,:) = Zero
end if
if (nLapl /= 0) then
  call mma_allocate(Lapl,nLapl,nGridMax,Label='Lapl')
  call mma_allocate(vLapl,nLapl,nGridMax,Label='vLapl')
  call mma_allocate(dfunc_dLapl,nLapl,nGridMax,Label='dfunc_dLapl')
  Lapl(:,:) = Zero
end if

call mma_allocate(F_xc,nGridMax,Label='F_xc')
call mma_allocate(func,nGridMax,Label='func')
if (l_casdft) then
  call mma_allocate(F_xca,nGridMax,Label='F_xca')
  call mma_allocate(F_xcb,nGridMax,Label='F_xcb')
end if

call mma_allocate(List_S,2,nIrrep*nShell,Label='List_S')
call mma_allocate(List_Exp,nIrrep*nShell,Label='List_Exp')
call mma_allocate(List_Bas,2,nIrrep*nShell,Label='List_Bas')
call mma_allocate(List_P,nNQ,Label='List_P')
call mma_allocate(R2_trial,nNQ,Label='R2_trial')

if (Do_MO) then
  if (NQNAC /= 0) then
    nD1MO = NQNACPAR
    call mma_allocate(D1MO,nD1MO,Label='D1MO')
    call Get_dArray_chk('D1mo',D1MO,nD1MO)
    nP2 = NQNACPR2
    call mma_Allocate(P2MO,nP2,Label='P2MO')
    call Get_dArray_chk('P2mo',P2MO,nP2)
  end if
  call Get_iArray('nBas',nBas,mIrrep)
  call Get_iArray('nDel',nDel,mIrrep)
  nCMO = 0
  do i=1,mIrrep
    nCMO = nCMO+nBas(i)*(nBas(i)-nDel(i))
  end do
  call mma_allocate(CMO,nCMO,Label='CMO')
  call Get_dArray_chk('Last orbitals',CMO,nCMO)
  call Get_iArray('nAsh',nAsh,mIrrep)
  nMOs = 0
  do iIrrep=0,mIrrep-1
    nMOs = nMOs+mBas(iIrrep)
  end do
end if

! Prepare memory for two-electron integrals:
! nPUVX

if (Do_TwoEl) then
  if (.not. Do_MO) then
    call WarningMessage(2,' Can''t produce 2 el dft integrals without MO')
    call Abend()
  end if

  iStack = 0
  do iIrrep=0,mIrrep-1
    iOrb = mBas(iIrrep)-nFro(iIrrep)
    do jIrrep=0,mIrrep-1
      jAsh = nAsh(jIrrep)
      ijIrrep = ieor(iIrrep,jIrrep)
      do kIrrep=0,mIrrep-1
        kAsh = nAsh(kIrrep)
        ijkIrrep = ieor(ijIrrep,kIrrep)
        if (ijkIrrep <= kIrrep) then
          lAsh = nAsh(ijkIrrep)
          kl_Orb_pairs = kAsh*lAsh
          if (kIrrep == ijkIrrep) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
          iStack = iStack+iOrb*jAsh*kl_Orb_pairs
        end if
      end do
    end do
  end do
  nTmpPUVX = iStack

end if

if (Do_Grad) then
  call mma_allocate(List_g,3,nShell*nIrrep,Label='List_G')
  mGrad = 3*nAtoms
  call mma_allocate(IndGrd,mGrad,Label='IndGrd')
  call mma_allocate(iTab,4,mGrad,Label='iTab')
  call mma_allocate(Temp,mGrad,Label='Temp')
end if

if (.not. Do_Grad) FckInt(:,:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'l_casdft value at drvnq:',l_casdft
if (l_casdft) write(u6,*) 'MCPDFT with functional:',KSDFA
#endif

if (l_casdft) then
  call mma_allocate(P2_ontop,nP2_ontop,nGridMax,Label='P2_ontop')
  P2_ontop(:,:) = Zero
end if

call DrvNQ_Inner(Kernel,Funct,Maps2p,nIrrep,List_S,List_Exp,List_bas,nShell,List_P,nNQ,FckInt,nFckDim,Density,nFckInt,nD,nGridMax, &
                 nP2_ontop,Do_Mo,nTmpPUVX,Do_Grad,Grad,nGrad,mAO,mdRho_dR)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate the memory

call mma_deallocate(Pax)
if (Do_Grad) then
  call mma_deallocate(Temp)
  call mma_deallocate(iTab)
  call mma_deallocate(IndGrd)
  call mma_deallocate(List_G)
end if
call mma_deallocate(R2_trial)
call mma_deallocate(List_P)
call mma_deallocate(List_Bas)
call mma_deallocate(List_Exp)
call mma_deallocate(List_S)
! Do_TwoEl
if (allocated(D1MO)) call mma_deallocate(D1MO)
if (allocated(P2MO)) call mma_deallocate(P2MO)
if (allocated(CMO)) call mma_deallocate(CMO)
if (l_casdft) then
  call mma_deallocate(F_xcb)
  call mma_deallocate(F_xca)
end if
call mma_deallocate(func)
call mma_deallocate(F_xc)
!
if (allocated(Lapl)) then
  call mma_deallocate(dfunc_dLapl)
  call mma_deallocate(vLapl)
  call mma_deallocate(Lapl)
end if
if (allocated(Tau)) then
  call mma_deallocate(dfunc_dTau)
  call mma_deallocate(vTau)
  call mma_deallocate(Tau)
end if
if (allocated(GradRho)) call mma_deallocate(GradRho)
if (allocated(Sigma)) then
  call mma_deallocate(dfunc_dSigma)
  call mma_deallocate(vSigma)
  call mma_deallocate(Sigma)
end if
call mma_deallocate(dfunc_dRho)
call mma_deallocate(vRho)
call mma_deallocate(Rho)

call mma_deallocate(Weights)
call mma_deallocate(Grid)

#ifdef _DEBUGPRINT_
write(u6,*) 'l_casdft value at drvnq:',l_casdft
if (l_casdft) write(u6,*) 'MCPDFT with functional:',KSDFA
#endif
if (allocated(P2_ontop)) call mma_deallocate(P2_ontop)

call mma_deallocate(nR_Eff)
call mma_deallocate(Coor)

call Close_NQ_Data()
call mma_deallocate(Mem)
call mma_deallocate(Angular)
call mma_deallocate(Fact)
call mma_deallocate(Maps2p)
NQ_Status = Inactive
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the status flag and TOC.

if ((iGrid_Set == Intermediate) .and. (Grid_Status == Regenerate)) iDisk_Set(Final_Grid) = iDisk_Grid
if (Do_Grad) then
  G_S(iGrid_Set) = Regenerate
else
  G_S(iGrid_Set) = Use_Old
end if

iDisk_Grid = 0
call iDaFile(Lu_Grid,1,G_S,2,iDisk_Grid)
call iDaFile(Lu_Grid,1,iDisk_Set,2,iDisk_Grid)
iDum(1) = Old_Functional_Type
call iDaFile(Lu_Grid,1,iDum,1,iDisk_Grid)
iDisk_Grid = iDisk_Set(iGrid_Set)
call iDaFile(Lu_Grid,1,GridInfo,2*number_of_subblocks,iDisk_Grid)

call DaClos(Lu_Grid)

call mma_deallocate(GridInfo)
!                                                                      *
!***********************************************************************
!                                                                      *
call IniPkR8(PThr,PMode)

call xFlush(LuGridFile)
close(LuGridFile)

return

end subroutine DrvNQ
