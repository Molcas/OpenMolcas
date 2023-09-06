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
! Copyright (C) 1991, Roland Lindh                                     *
!               1991, Anders Bernhardsson                              *
!***********************************************************************

subroutine Drvh2(Hess,Temp,nHess,show)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradient with respect to the one-  *
!         electron hamiltonian and the overlap matrix. The former will *
!         be contracted with the "variational" first order density     *
!         matrix and the latter will be contracted with the generalized*
!         Fock matrix.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!             Anders Bernhardsson Dept. of Theoretical Chemistry,      *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nCnttp, nBas
use Symmetry_Info, only: nIrrep
use McK_interface, only: hss_kernel, mck_mem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nHess
real(kind=wp), intent(inout) :: Hess(nHess)
real(kind=wp), intent(out) :: Temp(nHess)
logical(kind=iwp), intent(in) :: show
#include "rctfld.fh"
integer(kind=iwp) :: i, iIrrep, nComp, nDens, nFock
real(kind=wp) :: TCpu1, TCpu2, TWall1, TWall2
character(len=80) :: Label
logical(kind=iwp) :: DiffOp, lECP
integer(kind=iwp), allocatable :: lOper(:)
real(kind=wp), allocatable :: Coor(:,:), D0(:), Fock(:)
procedure(hss_kernel) :: KneHss, M1Hss, NAHss, OvrHss, PCMHss, PrjHss, SROHss
procedure(mck_mem) :: KneMmH, M1MmH, NAMmH, OvrMmH, PCMMmH, PrjMmH, SROMmH

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)
call StatusLine(' McKinley:',' Computing 1-electron 2rd order derivatives')
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the variational density matrix and the occupied Fock matrix.

nFock = 0
nDens = 0
do iIrrep=0,nIrrep-1
  nFock = nFock+nTri_Elem(nBas(iIrrep))
  nDens = nDens+nTri_Elem(nBas(iIrrep))
end do

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
call mma_allocate(D0,nDens,Label='D0')
call Get_D1ao_Var(D0,nDens)
! Read the generalized Fock matrix
! Fock matrix in AO/SO basis
call mma_allocate(Fock,nFock,Label='Fock')
call Get_dArray_chk('FockOcc',Fock,nFock)
!                                                                      *
!***********************************************************************
!                                                                      *
! Prologue
nComp = 1
call mma_allocate(Coor,3,nComp,Label='Coor')
Coor(:,:) = Zero
call mma_allocate(lOper,nComp,Label='lOper')
lOper(:) = 1
!***********************************************************************
!1)                                                                    *
!     Trace the generalized Fock matrix with the gradient of the       *
!     overlap matrix.                                                  *
!                                                                      *
!***********************************************************************

DiffOp = .false.
Temp(:) = Zero
Label = ' The Renormalization Contribution'
call Dot1El(OvrHss,OvrMmH,Temp,nHess,DiffOp,Coor,Fock,nFock,lOper,nComp)
if (show) write(u6,*) label
if (show) call HssPrt(Temp,nHess)
Hess(:) = Hess(:)-Temp(:)

!***********************************************************************
!2)                                                                    *
!     Trace the "variational" zero order density matrix with the       *
!     gradient of the kinetic energy integrals.                        *
!                                                                      *
!***********************************************************************

DiffOp = .false.
Temp(:) = Zero
Label = ' The Kinetic Energy Contribution'
call Dot1El(KneHss,KneMmH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
if (show) write(u6,*) label
if (show) call HssPrt(Temp,nHess)
Hess(:) = Hess(:)-Temp(:)

!***********************************************************************
!3)                                                                    *
!     Trace the "variational" zero order density matrix with the       *
!     gradient of the nuclear attraction integrals.                    *
!                                                                      *
!***********************************************************************

DiffOp = .true.
Label = ' The Nuclear Attraction Contribution'
Temp(:) = Zero
call Dot1El(NAHss,NAMmH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
if (show) write(u6,*) label
if (show) call HssPrt(Temp,nHess)
Hess(:) = Hess(:)+Temp(:)

!***********************************************************************
!3)                                                                    *
!     Trace the "variational" zero order density matrix with the       *
!     gradient of the ECP integrals.                                   *
!                                                                      *
!***********************************************************************

lECP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
end do
if (lECP) then
  DiffOp = .true.
  Label = ' The Projection (ECP) Contribution'
  Temp(:) = Zero
  call Dot1El(PrjHss,PRJMMH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
  if (show) write(u6,*) label
  if (show) call HssPrt(Temp,nHess)
  Hess(:) = Hess(:)+Temp(:)

  DiffOp = .true.
  Label = ' The Spec. Res. (ECP) Contribution'
  Temp(:) = Zero
  call Dot1El(SROHss,SROMMH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
  if (show) write(u6,*) Label,'first part '
  if (show) call HssPrt(Temp,nHess)
  Hess(:) = Hess(:)+Temp(:)

  DiffOp = .true.
  Label = ' The M1 (ECP) Contribution'
  Temp(:) = Zero
  call Dot1El(m1Hss,m1MMH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
  if (show) write(u6,*) Label,'second part '
  if (show) call HssPrt(Temp,nHess)
  Hess(:) = Hess(:)+Temp(:)
end if

!***********************************************************************
!4)                                                                    *
!     Trace the "variational" zero order density matrix with the       *
!     gradient of the PCM integrals.                                   *
!                                                                      *
!***********************************************************************

if (PCM) then
  DiffOp = .true.
  Label = ' The PCM Contribution'
  Temp(:) = Zero
  call Dot1El(PCMHss,PCMMMH,Temp,nHess,DiffOp,Coor,D0,nFock,lOper,nComp)
  if (show) write(u6,*) label
  if (show) call HssPrt(Temp,nHess)
  Hess(:) = Hess(:)+Temp(:)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(lOper)
call mma_deallocate(Coor)
call mma_deallocate(Fock)
call mma_deallocate(D0)
if (Show) call HssPrt(Hess,nHess)
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)

return

end subroutine Drvh2
