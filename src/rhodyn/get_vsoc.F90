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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************

subroutine get_vsoc()
!***********************************************************************
! Purpose:  transformation of V(SOC) matrix from the basis
!           of spin-free states to the basis of CSFs
!***********************************************************************
! V_SO    : SO-Hamiltonian matrix
! REV_SO  : Real part of V_SO
! IMV_SO  : Imaginary part of V_SO
! REV_CSF : U_CI*REV_SO*U_CI**T
! IMV_CSF : U_CI*IMV_SO*U_CI**T

use rhodyn_data, only: lrootstot, nconftot, ipglob, threshold, V_SO, V_CSF, U_CI, prep_vcsfr, prep_vcsfi, int2real
use rhodyn_utils, only: transform, dashes
use definitions, only: wp, iwp, u6
use stdalloc, only: mma_allocate, mma_deallocate
use mh5, only: mh5_put_dset

implicit none
real(kind=wp), dimension(:,:), allocatable :: REV_SO, REV_CSF, IMV_SO, IMV_CSF
integer(kind=iwp) :: i, j

if (ipglob > 2) write(u6,*) 'Begin of get_vsoc'
if (ipglob > 2) call dashes()

call mma_allocate(REV_SO,lrootstot,lrootstot)
call mma_allocate(IMV_SO,lrootstot,lrootstot)
call mma_allocate(REV_CSF,nconftot,nconftot)
call mma_allocate(IMV_CSF,nconftot,nconftot)

REV_SO(:,:) = dble(V_SO)
IMV_SO(:,:) = aimag(V_SO)

! check whether V_SO in SF basis is hermitian.
if (ipglob > 3) then
  call dashes()
  write(u6,*) 'Check if Spin-orbit coupling V_SO is hermitain'
  call dashes()
  do i=1,lrootstot
    do j=1,i
      if ((abs(REV_SO(i,j)-REV_SO(j,i)) >= threshold) .or. (abs(IMV_SO(i,j)+IMV_SO(j,i)) >= threshold)) then
        write(u6,*) 'ERROR: V_SO is not Hermitian; check element',i,j,REV_SO(i,j),REV_SO(j,i),IMV_SO(i,j),IMV_SO(j,i), &
                    (REV_SO(i,j)-REV_SO(j,i)),(IMV_SO(i,j)+IMV_SO(j,i))
      end if
    end do
  end do
  call dashes()
  write(u6,*) 'If there is no error printout, V_SO is hermitain!'
  call dashes()
end if

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

! Check whether V_CSF is hermitian
if (ipglob > 3) then
  call dashes()
  write(u6,*) 'Check whether SO-Hamiltonian in CSF is hermitian'
  call dashes()
  do i=1,nconftot
    do j=i+1,nconftot
      if (abs(REV_CSF(i,j)-REV_CSF(j,i)) >= threshold) then
        write(u6,int2real) 'WARNING!!!: REV_CSF is not hermitian:',i,j,REV_CSF(i,j),REV_CSF(j,i)
      else if (abs(IMV_CSF(i,j)+IMV_CSF(j,i)) >= threshold) then
        write(u6,int2real) 'WARNING!!!: IMV_CSF is not hermitian:',i,j,IMV_CSF(i,j),IMV_CSF(j,i)
      end if
    end do
  end do
  call dashes()
  write(u6,*) 'If there is no WARNING!!! printout'
  write(u6,*) 'SO-Hamiltonian in CSF basis is hermitian'
  call dashes()
end if

V_CSF(:,:) = dcmplx(REV_CSF,IMV_CSF)

! Store pure SOC matrix in CSF basis to PREP file
call mh5_put_dset(prep_vcsfr,REV_CSF)
call mh5_put_dset(prep_vcsfi,IMV_CSF)

if (allocated(REV_SO)) call mma_deallocate(REV_SO)
if (allocated(IMV_SO)) call mma_deallocate(IMV_SO)
if (allocated(REV_CSF)) call mma_deallocate(REV_CSF)
if (allocated(IMV_CSF)) call mma_deallocate(IMV_CSF)

if (ipglob > 2) write(u6,*) 'End of get_vsoc'

end subroutine get_vsoc
