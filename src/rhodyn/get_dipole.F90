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

subroutine get_dipole()
! Purpose :  Read in the dipole matrix from the MOLCAS output (SO)

use rhodyn_data, only: a_einstein, dipole, dysamp, dysamp_bas, E_SO, emiss, flag_dyson, flag_emiss, ispin, lroots, &
                       lrootstot, N, prep_dipolei, prep_dipoler, prep_do, runmode, SO_CI
use rhodyn_utils, only: transform
use mh5, only: mh5_put_dset
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, j, k, l, ii, jj

! To put the imaginary part of diagonal elements to 0 just in case
do j=1,lrootstot
  dipole(j,j,:) = real(dipole(j,j,:))
end do
! process Dyson amplitudes
if (flag_dyson .and. (N > 2)) then
  ! Transformation of Dyson amplitudes matrix to SF basis
  call transform(cmplx(dysamp,kind=wp),SO_CI,dysamp_bas,.false.)
  ! nullify non-neighbouring SF blocks, ii,jj - block indices
  ii = 0
  jj = 0
  do k=1,N
    do l=1,k
      if ((k-l) > 1) then
        do i=ii,(ii+lroots(k)*ispin(k))
          do j=jj,(jj+lroots(l)*ispin(l))
            dysamp_bas(i,j) = cZero
            dysamp_bas(j,i) = cZero
          end do
        end do
      end if
      jj = jj+lroots(l)*ispin(l)
    end do
    ii = ii+lroots(k)*ispin(k)
    jj = 0
  end do
else if (flag_dyson) then
  dysamp_bas(:,:) = dysamp
end if
if (flag_dyson) then
  dysamp_bas(:,:) = abs(dysamp_bas**2)
end if

! calculate matrix of Einstein coefficient A if emission spectrum needed
!if ((DM_basis /= 'SO') .and. (DM_basis /= 'CSF_SO') .and. (DM_basis /= 'SF_SO')) then
!  flag_emiss=.False.
!  write(u6,*) 'Emission spectra can be calculated only if SOC  density matrix available. '// &
!              'Set DMBasis keyword to SO, CSF_SO, or SF_SO'
!end if
if (flag_emiss) then
  a_einstein = Zero
  emiss = Zero
  ii = 1
  do j=1,(lrootstot-1)
    do i=(j+1),lrootstot
      do k=1,3
        a_einstein(i,j) = a_einstein(i,j)+abs(dipole(i,j,k))**2
      end do
      ! write down frequencies
      emiss(ii) = abs(E_SO(i)-E_SO(j))
      a_einstein(i,j) = a_einstein(i,j)*emiss(ii)**3
      ii = ii+1
    end do
  end do
end if

if (runmode /= 4) then
  ! not CM case
  call mh5_put_dset(prep_dipoler,real(dipole))
  call mh5_put_dset(prep_dipolei,aimag(dipole))
  if (flag_dyson) then
    call mh5_put_dset(prep_do,real(dysamp_bas))
  end if
end if

end subroutine get_dipole
