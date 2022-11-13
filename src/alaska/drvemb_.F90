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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine DrvEMB_(nh1,KSDFT,Do_Grad,Grad,nGrad,DFTFOCK)
!***********************************************************************
!***********************************************************************
!** Orbital-Free Embedding calculation (gradients)                   ***
!**                                                                  ***
!**                                                                  ***
!** Author: F. Aquilante, Geneva Nov  2010                           ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************

use OFembed, only: Xsigma, dFMD
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad
character(len=*), intent(inout) :: KSDFT
logical(kind=iwp), intent(in) :: Do_Grad
real(kind=wp), intent(out) :: Grad(nGrad)
character(len=4), intent(in) :: DFTFOCK
#include "debug.fh"
integer(kind=iwp) :: i, iSpin, kSpin, nD, nFckDim
real(kind=wp) :: d_Alpha, d_Beta, DSpn, DTot, Energy_NAD, Fakt_, Func_A, Func_AB, Func_B, Func_X
real(kind=wp), allocatable :: Grad_A(:), F_DFT(:,:), D_DS(:,:), Fcorr(:,:)
real(kind=wp), external :: Xlambda

Debug = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Do_Grad) then
  call WarningMessage(2,'DrvEMB_: Do_Grad must be .true.')
  call Abend()
end if
Grad(:) = Zero
call mma_allocate(Grad_A,nGrad,Label='Grad_A')
Grad_A(:) = Zero
!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys B (environment)             *
!                                                                      *
!***********************************************************************
call NameRun('AUXRFIL') ! switch RUNFILE name
!                                                                      *
!***********************************************************************
!                                                                      *
nD = 4
call mma_allocate(F_DFT,nh1,nD,Label='F_DFT')
call mma_allocate(D_DS,nh1,nD,Label='D_DS')

! Get the density matrix of the environment (rho_B)

call Get_iScalar('Multiplicity',kSpin)
call Get_dArray_chk('D1ao',D_DS(1,1),nh1)
!call RecPrt('D_DS(1,1)',' ',D_DS(1,1),nh1,1)

! Get the spin density matrix of the environment

if (kSpin /= 1) then
  call Get_dArray_chk('D1sao',D_DS(1,2),nh1)
  !call RecPrt('D1Sao',' ',D_DS(1,2),nh1,1)
end if

! Compute alpha and beta density matrices of the environment

nFckDim = 2
if (kSpin == 1) then
  D_DS(:,1) = Half*D_DS(:,1)
  D_DS(:,2) = D_DS(:,1)
  nFckDim = 1
else
  do i=1,nh1
    DTot = D_DS(i,1)
    DSpn = D_DS(i,2)
    d_Alpha = Half*(DTot+DSpn)
    d_Beta = Half*(DTot-DSpn)
    D_DS(i,1) = d_Alpha
    D_DS(i,2) = d_Beta
  end do
  !call RecPrt('Da',' ',D_DS(1,1),nh1,1)
  !call RecPrt('Db',' ',D_DS(1,2),nh1,1)
end if

if (KSDFT(1:4) == 'NDSD') then

  call wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func_B,D_DS,nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

  KSDFT(1:4) = 'LDTF' !set to Thomas-Fermi for subsequent calls
end if

!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys A                           *
!                                                                      *
!***********************************************************************
call NameRun('#Pop')     ! switch back RUNFILE name

! Get the density matrix for rho_A

call Get_dArray_chk('D1ao',D_DS(1,3),nh1)
!call RecPrt('D_DS(1,3)',' ',D_DS(1,3),nh1,1)

call Get_iScalar('Multiplicity',iSpin)
if ((iSpin == 1) .and. (kSpin /= 1)) then
  call WarningMessage(0,' Non-singlet environment perturbation on singlet state!'// &
                      ' Spin-components of the OFE potential will be averaged.')
end if

! Get the spin density matrix of A

if (iSpin /= 1) then
  call Get_dArray_chk('D1sao',D_DS(1,4),nh1)
  !call RecPrt('D1Sao',' ',D_DS(1,4),nh1,1)
end if

! Compute alpha and beta density matrices of subsystem A

nFckDim = 2
if (iSpin == 1) then
  D_DS(:,3) = Half*D_DS(:,3)
  D_DS(:,4) = D_DS(:,3)
  if (kSpin == 1) nFckDim = 1
else
  do i=1,nh1
    DTot = D_DS(i,3)
    DSpn = D_DS(i,4)
    d_Alpha = Half*(DTot+DSpn)
    d_Beta = Half*(DTot-DSpn)
    D_DS(i,3) = d_Alpha
    D_DS(i,4) = d_Beta
  end do
  !call RecPrt('Da',' ',D_DS(1,3),nh1,1)
  !call RecPrt('Db',' ',D_DS(1,4),nh1,1)
end if

call wrap_DrvNQ(KSDFT,F_DFT(1,3),nFckDim,Func_A,D_DS(1,3),nh1,nFckDim,Do_Grad,Grad_A,nGrad,DFTFOCK)

call daxpy_(nGrad,-One,Grad_A,1,Grad,1)

! Fraction of correlation potential from A (cases: HF or Trunc. CI)
if (dFMD > Zero) then

  Grad_A(:) = Zero
  call mma_allocate(Fcorr,nh1,nFckDim,Label='Fcorr')

  call cwrap_DrvNQ(KSDFT,nFckDim,Func_A,D_DS(1,3),nh1,nFckDim,Do_Grad,Grad_A,nGrad,DFTFOCK,Fcorr)

  call get_dScalar('NAD dft energy',Energy_NAD)
  Fakt_ = Xlambda(abs(Energy_NAD),Xsigma)
  call daxpy_(nGrad,Fakt_,Grad_A,1,Grad,1)

  call mma_deallocate(Fcorr)
end if

call mma_deallocate(Grad_A)

call NameRun('AUXRFIL') ! switch RUNFILE name

call wrap_DrvNQ('NUCATT_EMB',F_DFT,nFckDim,Func_X,D_DS(1,3),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

call NameRun('#Pop')    ! switch back RUNFILE name

!***********************************************************************
!                                                                      *
!     Calculation on the supermolecule                                 *
!                                                                      *
!***********************************************************************
nFckDim = 2
if ((iSpin == 1) .and. (kSpin == 1)) then
  nFckDim = 1
  call daxpy_(nh1,One,D_DS(1,3),1,D_DS(1,1),1)
else
  call daxpy_(nh1,One,D_DS(1,3),1,D_DS(1,1),1)
  call daxpy_(nh1,One,D_DS(1,4),1,D_DS(1,2),1)
end if

call wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func_AB,D_DS,nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

call mma_deallocate(F_DFT)
call mma_deallocate(D_DS)

return

end subroutine DrvEMB_
