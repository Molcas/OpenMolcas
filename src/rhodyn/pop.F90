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

subroutine pop(time,popcount)
!***********************************************************************
! prints diagonal of the density matrix densityt in reqired basis
! at the current time
!***********************************************************************

use rhodyn_data, only: basis, DM_basis, density0, densityt, d, dgl, out_fmt_csf, out_fmt, out_tout, out_dm_csf, out_dm_sf, &
                       out_dm_so, out_emiss, emiss, n_freq, a_einstein, lu_csf, lu_sf, lu_so, lu_dip, tmp, pulse_vec, flag_emiss, &
                       flag_dipole, dipole_basis, SO_CI, CSF2SO, U_CI_compl, Nstate, nconftot
use rhodyn_utils, only: mult, transform
use definitions, only: wp, iwp
use constants, only: auToFs
use mh5, only: mh5_put_dset

implicit none
integer(kind=iwp), intent(in) :: popcount
real(kind=wp), intent(in) :: time
real(kind=wp) :: norm
real(kind=wp), dimension(nconftot) :: dgl_csf
complex(kind=wp), dimension(nconftot,nconftot) :: density_csf
character(len=64) :: sline
integer(kind=iwp) :: i, j, l

!   here notation d is dimension of all the basis matrices
!!! density0 (can't be) used as a temporary storage for dm in required basis

write(sline,'(f10.3)') time*auToFs
call StatusLine('SpinDyn current time: ',trim(sline))

call mh5_put_dset(out_tout,[time*auToFs],[1],[popcount-1])

if (basis == 'CSF') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in CSF basis
    dgl_csf = dble((/(densityt(i,i),i=1,d)/))
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt) time*auToFs,(dgl_csf(i),i=1,d),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from CSF to SF states basis
    call transform(densityt,U_CI_compl,density0)
    dgl(:) = dble((/(density0(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from CSF to SO states basis
    call transform(densityt,CSF2SO,density0)
    dgl(:) = dble((/(density0(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if

else if (basis == 'SF') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SF to CSF states basis
    call transform(densityt,U_CI_compl,density_csf,.false.)
    dgl_csf = dble((/(density_csf(i,i),i=1,nconftot)/))
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt_csf) time*auToFs,(dgl_csf(i),i=1,nconftot),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,nconftot],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in SF basis
    dgl(:) = dble((/(densityt(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SF to SO states basis
    call transform(densityt,SO_CI,density0)
    dgl(:) = dble((/(density0(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if

else if (basis == 'SO') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SO to CSF states basis
    call transform(densityt,CSF2SO,density_csf,.false.)
    dgl_csf = dble((/(density_csf(i,i),i=1,nconftot)/))
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt_csf) time*auToFs,(dgl_csf(i),i=1,nconftot),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,nconftot],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SO to SF states basis
    call transform(densityt,SO_CI,density0,.false.)
    dgl(:) = dble((/(density0(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in SO basis
    dgl(:) = dble((/(densityt(i,i),i=1,d)/))
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if
end if

! time-dependent dipole moment
! pulse_vec is used as storage for tr(rho d)
if (flag_dipole) then
  do i=1,3
    call mult(densityt,dipole_basis(:,:,i),tmp)
    dgl(:) = dble((/(tmp(i,i),i=1,d)/))
    pulse_vec(i) = sum(dgl)
  end do
  write(lu_dip,'(7(g25.15e3,2x))') time*auToFs,(dble(pulse_vec(i)),aimag(pulse_vec(i)),i=1,3)
end if

! emission spectra
if (flag_emiss) then
  ! in dgl the SOC populations are left
  l = 1
  emiss = 0d0
  do j=1,(Nstate-1)
    do i=(j+1),Nstate
      emiss(l) = a_einstein(i,j)*dgl(i)
      l = l+1
    end do
  end do
  call mh5_put_dset(out_emiss,emiss,[1,n_freq],[popcount-1,0])
end if

end subroutine pop
