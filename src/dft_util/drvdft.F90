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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

subroutine DrvDFT(h1,nh1,KSDFT,ExFac,Do_Grad,Grad,nGrad,iSpin,DFTFOCK)

use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc, KSDFA
use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Energy_integrated, Grad_I, mBas, mIrrep, nFro, &
                   nIsh, Tau_I
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad, iSpin
real(kind=wp), intent(inout) :: h1(nh1), Grad(nGrad)
character(len=*), intent(in) :: KSDFT
real(kind=wp), intent(out) :: ExFac
logical(kind=iwp), intent(in) :: Do_Grad
character(len=4), intent(in) :: DFTFOCK
#include "debug.fh"
integer(kind=iwp) :: i, nD, nFckDim
real(kind=wp) :: d_Alpha, d_Beta, DSpn, DTot, Fact, Func, PDFT_Ratio, Vxc_ref(2), WF_Ratio
logical(kind=iwp) :: Do_HPDFT, Do_MO, Do_TwoEl
real(kind=wp), allocatable :: D_DS(:,:), F_DFT(:,:)
real(kind=wp), external :: DDot_, Get_ExFac

!                                                                      *
!***********************************************************************
!                                                                      *
KSDFA = KSDFT  ! Store a local copy
Debug = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call Put_iScalar('Multiplicity',iSpin)
call Get_iScalar('nSym',mIrrep)
call Get_iArray('nBas',mBas,mIrrep)

call Set_Basis_Mode('Valence')
call Setup_iSD()

call Get_dScalar('DFT exch coeff',CoefX)
call Get_dScalar('DFT corr coeff',CoefR)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_Grad) Grad(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
if (iSpin == 1) then
  nD = 1
else
  nD = 2
end if

! What is this?

if (DFTFOCK == 'ROKS') nD = 2
call mma_allocate(D_DS,nh1,nD,Label='D_DS')

! Get the total density

call Get_dArray_chk('D1ao',D_DS,nh1)
!call RecPrt('D1ao',' ',D_DS(:,1),nh1,1)

! Get the spin density

if (nD /= 1) then
  call Get_dArray_chk('D1sao',D_DS(:,2),nh1)
  !call RecPrt('D1Sao',' ',D_DS(:,2),nh1,1)
end if

! Compute alpha and beta densities

!call RecPrt('DTot',' ',D_DS(:,1),nh1,1)
!call RecPrt('DSpn',' ',D_DS(:,2),nh1,1)
if (nD == 1) then
  D_DS(:,1) = Half*D_DS(:,1)
else
  do i=1,nh1
    DTot = D_DS(i,1)
    DSpn = D_DS(i,2)
    d_Alpha = Half*(DTot+DSpn)
    d_Beta = Half*(DTot-DSpn)
    D_DS(i,1) = d_Alpha
    D_DS(i,2) = d_Beta
  end do
end if
!call RecPrt('Da',' ',D_DS(:,1),nh1,1)
!call RecPrt('Db',' ',D_DS(:,2),nh1,1)

if (KSDFT(1:3) /= 'SCF') then
  call Get_iArray('nIsh',nIsh,mIrrep)
  call Get_iArray('nFro',nFro,mIrrep)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! DFT functionals, compute integrals over the potential

Func = Zero
Funcaa = Zero
Funcbb = Zero
Funccc = Zero
Dens_I = Zero
Dens_a1 = Zero
Dens_b1 = Zero
Dens_a2 = Zero
Dens_b2 = Zero
Dens_t1 = Zero
Dens_t2 = Zero
Grad_I = Zero
Tau_I = Zero
Do_MO = .false.
Do_TwoEl = .false.

! nFckDim: number of different types of Fock matrices. Normally for
! conventional functionals we have one Fock matrix for closed shell
! calculations and two (F_alpha and F_beta) for open shell systems.
! For CASDFT we have always two (F_inactive and F_active)

nFckDim = nD
call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
F_DFT(:,:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call Driver(KSDFA,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nD,DFTFOCK)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_Grad) then
  Do_HPDFT = .false.
  call qpg_DScalar('R_WF_HMC',Do_HPDFT)
  if (Do_HPDFT) then
    write(u6,*) 'DFT gradient is scaled in a hybrid formalism.'
    call Get_DScalar('R_WF_HMC',WF_Ratio)
    PDFT_Ratio = One-WF_Ratio
    Grad(:) = PDFT_Ratio*Grad
  end if
end if

!                                                                      *
!***********************************************************************
!                                                                      *
ExFac = Get_ExFac(KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
Energy_integrated = Func
!                                                                      *
!***********************************************************************
!                                                                      *
if ((KSDFT == 'Overlap') .or. (KSDFT == 'NucAtt')) then
  h1(:) = F_DFT(:,1)
  if (KSDFT == 'NucAtt') Energy_integrated = Func
else

  ! Put out the integrated DFT energy and the DFT Fock matrices on the RUNFILE

  !call Put_DFT_Energy(Energy_integrated)
  call Poke_dScalar('KSDFT energy',Energy_integrated)
  call Put_dScalar('CASDFT energy',Energy_integrated)
  call Put_dArray('dExcdRa',F_DFT,nFckDim*nh1)
  !write(u6,'(a,f22.16)') ' Energy in drvdft ',Energy_integrated
# ifdef _DEBUGPRINT_
  write(u6,'(a,f22.16)') ' Energy ',Energy_integrated
  if (nFckDim == 1) then
    do i=1,nh1
      write(u6,'(i4,f22.16)') i,F_DFT(i,1)
    end do
  else
    do i=1,nh1
      write(u6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2),Half*(F_DFT(i,1)+F_DFT(i,2))
    end do
  end if
# endif

  ! In the SCF program (traclc.f) the program computes the trace
  ! of the one-electron hamiltonian over a set of densities. The
  ! DFT contribution is not linear with respect to variations of
  ! the density. However, with the following term we can include
  ! the linear component in that code.

  Fact = Two
  if (nD /= 1) Fact = One
  Vxc_ref(1) = Fact*DDot_(nh1,F_DFT(:,1),1,D_DS,1)
  if (nD /= 1) then
    Vxc_ref(2) = DDot_(nh1,F_DFT(:,2),1,D_DS(:,2),1)
  else
    Vxc_ref(2) = Zero
  end if
  call Put_Temp('Vxc_ref ',Vxc_ref,2)
end if

call mma_deallocate(F_DFT)
call mma_deallocate(D_DS)
call Free_iSD()

return

end subroutine DrvDFT
