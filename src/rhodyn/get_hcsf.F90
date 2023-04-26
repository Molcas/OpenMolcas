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

subroutine get_hcsf()
!***********************************************************************
!
! Purpose: Construct the H(RASSCF) matrix in the total CSF basis of
!          different spin manifolds. H(RASSCF) has been already read
!          in the matrix H_CSF(:,:,N) for different spin manifolds
!
!***********************************************************************

use rhodyn_data, only: flag_so, H_CSF, HTOT_CSF, HTOTRE_CSF, ipglob, ispin, N, nconf, nconftot, prep_fhi, prep_fhr, threshold, &
                       V_CSF
use rhodyn_utils, only: dashes, check_hermicity
use mh5, only: mh5_put_dset
use Constants, only: Zero, cZero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, l, ii, jj, kk

! Construct the Hamiltonian matrix H(RASSCF) in CSF basis include all
! spin manifolds to HTOTRE_CSF
HTOTRE_CSF(:,:) = Zero
ii = 0
jj = 0
kk = 0
if (flag_so) then
  do i=1,N
    if (i /= 1) ii = ii+nconf(i-1)*ispin(i-1)
    do j=1,nconf(i)
      do k=1,nconf(i)
        do l=1,ispin(i)
          jj = ii+(j-1)*ispin(i)+l
          kk = ii+(k-1)*ispin(i)+l
          HTOTRE_CSF(jj,kk) = H_CSF(j,k,i)
        end do
      end do
    end do
  end do
else
  ! without accounting for spin-degeneracy
  do i=1,N
    if (i /= 1) ii = ii+nconf(i-1)
    do j=1,nconf(i)
      do k=1,nconf(i)
        jj = ii+j
        kk = ii+k
        HTOTRE_CSF(jj,kk) = H_CSF(j,k,i)
      end do
    end do
  end do
end if

if (ipglob > 2) then
  write(u6,*) 'Construct the full Hamiltonian with (possibly) SOC!'
  call dashes()
  !write(u6,sint) 'number of NCSF:',nconftot
  !call dashes()
end if
HTOT_CSF(:,:) = cZero

! If consider the spin-orbit coupling
! Computing the total Hamiltonian H(RASSCF)+H(SO) in CSF basis
! to HTOT_CSF
if (.not. flag_so) then
  HTOT_CSF(:,:) = HTOTRE_CSF
else
  HTOT_CSF(:,:) = HTOTRE_CSF+V_CSF
end if
if (ipglob > 2) write(u6,*) 'end constructing full Hamiltonian'

! dimensionality of this check should be verified (is it always of nconftot size)
if (ipglob > 3) call check_hermicity(HTOT_CSF,nconftot,'Hamiltonian in CSF basis',threshold)

call mh5_put_dset(prep_fhr,real(HTOT_CSF))
call mh5_put_dset(prep_fhi,aimag(HTOT_CSF))

end subroutine get_hcsf
