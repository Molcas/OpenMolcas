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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine get_vsoc()
!***********************************************************************
! Purpose:  transformation of V(SOC) matrix from the basis
!           of spin-free states as it is obtained from RASSI
!           to the basis of CSFs
!***********************************************************************
! V_SO    : SO-Hamiltonian matrix
! REV_SO  : Real part of V_SO
! IMV_SO  : Imaginary part of V_SO
! REV_CSF : U_CI*REV_SO*U_CI**T
! IMV_CSF : U_CI*IMV_SO*U_CI**T

use rhodyn_data, only: ipglob, lrootstot, nconftot, prep_vcsfi, prep_vcsfr, U_CI, V_CSF, V_SO, threshold
use rhodyn_utils, only: dashes, transform, check_hermicity
use mh5, only: mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: REV_SO(:,:), REV_CSF(:,:), IMV_SO(:,:), IMV_CSF(:,:)

call mma_allocate(REV_SO,lrootstot,lrootstot)
call mma_allocate(IMV_SO,lrootstot,lrootstot)
call mma_allocate(REV_CSF,nconftot,nconftot)
call mma_allocate(IMV_CSF,nconftot,nconftot)

REV_SO(:,:) = real(V_SO)
IMV_SO(:,:) = aimag(V_SO)

call check_hermicity(V_SO,lrootstot,'V_SO in SF basis',threshold)

if (ipglob > 3) then
  call dashes()
  write(u6,*) 'Printout the Spin-orbit Hamiltonian in SF basis'
  call dashes()
  do i=1,6
    write(u6,*) (V_SO(i,j),j=1,6)
  end do
end if

if (ipglob > 2) write(u6,*) 'Begin transform the SO-Hamiltonian'
! Transform the SO-Hamiltonian from SF states to CSFs
call transform(REV_SO,U_CI,REV_CSF,.false.)
call transform(IMV_SO,U_CI,IMV_CSF,.false.)

V_CSF(:,:) = cmplx(REV_CSF,IMV_CSF,kind=wp)

call check_hermicity(V_CSF,nconftot,'V_SO in CSF basis',threshold)

! Store pure SOC matrix in CSF basis to PREP file
call mh5_put_dset(prep_vcsfr,REV_CSF)
call mh5_put_dset(prep_vcsfi,IMV_CSF)

if (allocated(REV_SO)) call mma_deallocate(REV_SO)
if (allocated(IMV_SO)) call mma_deallocate(IMV_SO)
if (allocated(REV_CSF)) call mma_deallocate(REV_CSF)
if (allocated(IMV_CSF)) call mma_deallocate(IMV_CSF)

end subroutine get_vsoc
